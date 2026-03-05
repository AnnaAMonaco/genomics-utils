#!/usr/bin/env bash
# Trim reads, align with STAR, generate BigWigs, and quantify with Salmon (optional)
# Requires: trim_galore, samtools >= 1.13, GNU parallel, STAR, salmon, deepTools
#   - <sample_list.txt> should list one sample ID per line (no extensions)
#   - Each sample must have paired FASTQs: ${sample}_R1.fastq.gz and ${sample}_R2.fastq.gz

set -euo pipefail

if [ "$#" -lt 4 ]; then
    echo -e "Trim reads, align with STAR, generate BigWigs, and quantify with Salmon (optional)\n\
    Usage: $0 <sample_list.txt> <fastq_dir> <output_dir> [star_index] [salmon_index] [num_threads] [--no-salmon]\n\
    Note: <sample_list.txt> should list one sample ID per line (no extensions)\n\
	To skip Salmon quantification, use --no-salmon.\n\
    To skip STAR and bigWig generation, use --no-star."
    exit 1
fi

SAMPLES=$1
FQDIR=$2
OUTDIR=$3
STARIDX=${4:-""}
SLMIDX=${5:-""}
THREADS=${6:-4}
RUN_SALMON=true
RUN_STAR=true

# detect if salmon or STAR should be skipped
shift 6 2>/dev/null || shift "$#" 2>/dev/null || true
for arg in "$@"; do
    case "$arg" in
        --no-salmon) RUN_SALMON=false ;;
        --no-star) RUN_STAR=false ;;
        *) echo "[ERROR] Unknown option: $arg"; exit 1 ;;
    esac
done

if ! $RUN_SALMON && ! $RUN_STAR; then
    echo "[ERROR] --no-salmon and --no-star cannot be used together."
    exit 1
fi

echo "[INFO] Pipeline configuration:"
echo "[INFO]   STAR alignment:   $RUN_STAR"
echo "[INFO]   Salmon quant:     $RUN_SALMON"
echo "[INFO]   Threads:          $THREADS"

# ensure all downstream directories are present
mkdir -p "$OUTDIR"/{data/{bamFiles,bigwig,counts},tmp/trim,logs}

# Export for GNU parallel subprocesses
export OUTDIR FQDIR STARIDX SLMIDX THREADS RUN_SALMON

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
            bamCoverage -b "$OUTDIR/data/bamFiles/${sample}.filtered.sorted.rmDup.bam" \
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
