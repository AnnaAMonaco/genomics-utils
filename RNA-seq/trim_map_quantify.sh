#!/usr/bin/env bash
# Trim reads, align with STAR, generate BigWigs, and quantify with Salmon (optional)
# Requires: trim_galore, samtools >= 1.13, GNU parallel, STAR, salmon, deepTools
#   - <sample_list.txt> should list one sample ID per line (no extensions)
#   - Each sample must have paired FASTQs: ${sample}_R1.fastq.gz and ${sample}_R2.fastq.gz

set -euo pipefail

# Deaful settings
THREADS=4
RUN_STAR=true
RUN_SALMON=true
STARIDX=""
SLMIDX=""
SAMPLES=""
FQDIR=""
OUTDIR=""

# Argument parsing
while [[ $# -gt 0 ]]; do
    case "$1" in
        -s)
            SAMPLES="$2"
            shift 2
            ;;
        -f)
            FQDIR="$2"
            shift 2
            ;;
        -o)
            OUTDIR="$2"
            shift 2
            ;;
        --star-index)
            STARIDX="$2"
            shift 2
            ;;
        --salmon-index)
            SLMIDX="$2"
            shift 2
            ;;
        -t)
            THREADS="$2"
            shift 2
            ;;
        --no-salmon)
            RUN_SALMON=false
            shift
            ;;
        --no-star)
            RUN_STAR=false
            shift
            ;;
        -h|--help)
            echo "Trim reads, align with STAR, generate BigWigs, and quantify with Salmon (optional)\n\
			Usage: $0 -s <sample_list.txt> -f <fastq_dir> -o <output_dir> --star-index [star_index] --salmon-index [salmon_index] -t [num_threads] [--no-salmon] [--no-star]\n"
            echo "  -s        		sample list - should list one sample ID per line (no extensions)"
            echo "  -f      		FASTQ directory"
            echo "  -o         		output directory"
            echo "  --star-index	STAR index path (skip with --no-star)"
            echo "  --salmon-index	Salmon index path (skip with --no-salmon)"
            echo "  -t        		number of threads (default: 4)"
            echo "  --no-star       skip STAR alignment"
            echo "  --no-salmon     skip Salmon quantification"
            exit 0
            ;;
        *)
            echo "[ERROR] Unknown option: $1"
            exit 1
            ;;
    esac
done

# Validate required arguments
if [[ -z "$SAMPLES" || -z "$FQDIR" || -z "$OUTDIR" ]]; then
    echo "[ERROR] --samples, --fastq-dir, and --outdir are required."
    exit 1
fi

if ! $RUN_STAR && ! $RUN_SALMON; then
    echo "[ERROR] Cannot disable both STAR and Salmon."
    exit 1
fi

if $RUN_STAR && [[ -z "$STARIDX" ]]; then
    echo "[ERROR] STAR enabled but no STAR index provided (--star-index)."
    exit 1
fi

if $RUN_SALMON && [[ -z "$SLMIDX" ]]; then
    echo "[ERROR] Salmon enabled but no Salmon index provided (--salmon-index)."
    exit 1
fi


if ! $RUN_SALMON && ! $RUN_STAR; then
    echo "[ERROR] --no-salmon and --no-star cannot be used together."
    exit 1
fi

# Show pipeline configuration
echo "[INFO] Pipeline configuration:"
echo "[INFO]   Sample list:       $SAMPLES"
echo "[INFO]   FASTQ directory:   $FQDIR"
echo "[INFO]   Output directory:  $OUTDIR"
echo "[INFO]   Threads:           $THREADS"
echo "[INFO]   STAR alignment:    $RUN_STAR"
echo "[INFO]   Salmon quant:      $RUN_SALMON"
$RUN_STAR   && echo "[INFO]   STAR index:        $STARIDX"
$RUN_SALMON && echo "[INFO]   Salmon index:      $SLMIDX"

# ensure all downstream directories are present
mkdir -p "$OUTDIR"/{data/{bamFiles,bigwig,counts},tmp/trim,logs}

# Export for GNU parallel subprocesses
export OUTDIR FQDIR STARIDX SLMIDX THREADS RUN_SALMON RUN_STAR

# Check dependencies
REQUIRED_CMDS=(trim_galore parallel)
$RUN_SALMON && REQUIRED_CMDS+=(salmon)
$RUN_STAR && REQUIRED_CMDS+=(samtools STAR bamCoverage)
for cmd in "${REQUIRED_CMDS[@]}"; do
    if ! command -v "$cmd" &>/dev/null; then
        echo "[ERROR] Required tool '$cmd' not found in PATH."
        exit 1
    fi
done

# Read sample IDs into an array
mapfile -t samples < <(grep -vE '^#|^[[:space:]]*$' "$SAMPLES")
echo "[INFO] Found ${#samples[@]} samples to process."


## Trim reads
echo "[INFO] Starting trimming with trim_galore..."
# Run trim_galore only for listed samples
trim_samples=()
for sample in "${samples[@]}"; do
    if [[ ! -f "$OUTDIR/tmp/trim/${sample}_val_1.fq.gz" || ! -f "$OUTDIR/tmp/trim/${sample}_val_2.fq.gz" ]]; then
        trim_samples+=("$sample")
    fi
done

if (( ${#trim_samples[@]} > 0 )); then
    printf "%s\n" "${trim_samples[@]}" | parallel -j 4 '
        trim_galore --illumina --paired --basename {1} \
            "$FQDIR"/{1}_R1.fastq.gz "$FQDIR"/{1}_R2.fastq.gz \
            -o "$OUTDIR"/tmp/trim > "$OUTDIR"/logs/{1}_trim.log 2>&1
    '
fi
# Only progress once both reads of each samples has been trimmed
echo "[INFO] Trimming complete for all samples."


## Define processing function
process_sample() {
    local sample=$1
    echo "[INFO] Processing ${sample}..."

    local R1="$OUTDIR/tmp/trim/${sample}_val_1.fq.gz"
    local R2="$OUTDIR/tmp/trim/${sample}_val_2.fq.gz"
    local STARprefix="$OUTDIR/data/bamFiles/${sample}"

    # salmon quant
    if $RUN_SALMON; then
        if [[ ! -f "$OUTDIR/data/counts/${sample}_quant/quant.sf" ]]; then
            salmon quant -i "$SLMIDX" -l A -1 "$R1" -2 "$R2" \
                -o "$OUTDIR/data/counts/${sample}_quant" > "$OUTDIR/logs/${sample}_salmon.log" 2>&1
        fi
    else
        echo "[INFO] Skipping Salmon quantification for ${sample}."
    fi

    # STAR alignment
    if $RUN_STAR; then
        if [[ ! -f "$OUTDIR/data/bamFiles/${sample}.filtered.sorted.bam" || ! -f "$OUTDIR/data/bamFiles/${sample}.filtered.sorted.bam.bai" ]]; then
            STAR --runThreadN "$THREADS" --genomeDir "$STARIDX" --readFilesIn "$R1" "$R2" \
                --outFileNamePrefix "$STARprefix" --outSAMunmapped Within --outSAMtype BAM Unsorted \
                --readFilesCommand zcat --alignEndsType Local > "$OUTDIR/logs/${sample}_STAR.log" 2>&1
    # BAM processing
            samtools view -@ "$THREADS" -b -F 2308 -q 20 "${STARprefix}Aligned.out.bam" | \
                samtools sort -@ "$THREADS" -n -T "$OUTDIR/tmp/${sample}" -o "$OUTDIR/tmp/${sample}.nameSorted.bam" -
            samtools fixmate -@ "$THREADS" "$OUTDIR/tmp/${sample}.nameSorted.bam" "$OUTDIR/tmp/${sample}.matefixed.bam"
            samtools sort -@ "$THREADS" -T "$OUTDIR/tmp/${sample}" \
                -o "$OUTDIR/data/bamFiles/${sample}.filtered.sorted.bam" "$OUTDIR/tmp/${sample}.matefixed.bam"
            samtools index "$OUTDIR/data/bamFiles/${sample}.filtered.sorted.bam"
        fi
    # BigWig generation
        if [[ ! -f "$OUTDIR/data/bigwig/${sample}.bw" ]]; then
            bamCoverage -b "$OUTDIR/data/bamFiles/${sample}.filtered.sorted.bam" \
                -o "$OUTDIR/data/bigwig/${sample}.bw" \
                --normalizeUsing RPKM > "$OUTDIR/logs/${sample}_bw.log" 2>&1
        fi
    else
        echo "[INFO] Skipping STAR alignment for ${sample}."
    fi

    echo "[INFO] Finished ${sample}."
}

export -f process_sample


## Run all samples in parallel
echo "[INFO] Running $( $RUN_STAR && echo 'STAR MEM + bamCoverage' || echo '(no STAR)' ) $( $RUN_SALMON && echo '+ Salmon' || echo '(no Salmon)' )..."
printf "%s\n" "${samples[@]}" | parallel -j 4 process_sample {}

echo "[DONE] Pipeline completed successfully. All outputs in $OUTDIR"
