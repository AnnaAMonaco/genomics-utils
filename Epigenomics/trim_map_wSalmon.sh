#!/usr/bin/env bash
# Trim reads, align with BWA, generate BigWigs, and quantify with Salmon
# Requires: trim_galore, samtools >= 1.13, GNU parallel, BWA, salmon, deepTools
#   - <sample_list.txt> should list one sample ID per line (no extensions)
#   - Each sample must have paired FASTQs: ${sample}_R1.fastq.gz and ${sample}_R2.fastq.gz

set -euo pipefail

# Deaful settings
THREADS=4
RUN_BWA=true
RUN_SALMON=true
BWAIDX=""
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
        --bwa-index)
            BWAIDX="$2"
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
        --no-bwa)
            RUN_BWA=false
            shift
            ;;
        -h|--help)
            echo "Trim reads, align with bwa mem, generate BigWigs, and quantify with Salmon (optional)\n\
			Usage: $0 -s <sample_list.txt> -f <fastq_dir> -o <output_dir> --bwa-index [bwa_index] --salmon-index [salmon_index] -t [num_threads] [--no-salmon] [--no-bwa]\n"
            echo "  -s        		sample list - should list one sample ID per line (no extensions)"
            echo "  -f      		FASTQ directory"
            echo "  -o         		output directory"
            echo "  --bwa-index		bwa index path (skip with --no-bwa)"
            echo "  --salmon-index	Salmon index path (skip with --no-salmon)"
            echo "  -t        		number of threads (default: 4)"
            echo "  --no-bwa       	skip bwa mem alignment"
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

if ! $RUN_BWA && ! $RUN_SALMON; then
    echo "[ERROR] Cannot disable both bwa and Salmon."
    exit 1
fi

if $RUN_BWA && [[ -z "$BWAIDX" ]]; then
    echo "[ERROR] bwa enabled but no bwa index provided (--bwa-index)."
    exit 1
fi

if $RUN_SALMON && [[ -z "$SLMIDX" ]]; then
    echo "[ERROR] Salmon enabled but no Salmon index provided (--salmon-index)."
    exit 1
fi

if ! $RUN_SALMON && ! $RUN_BWA; then
    echo "[ERROR] --no-salmon and --no-bwa cannot be used together."
    exit 1
fi

# Show pipeline configuration
echo "[INFO] Pipeline configuration:"
echo "[INFO]   Sample list:       $SAMPLES"
echo "[INFO]   FASTQ directory:   $FQDIR"
echo "[INFO]   Output directory:  $OUTDIR"
echo "[INFO]   Threads:           $THREADS"
echo "[INFO]   bwa alignment:     $RUN_BWA"
echo "[INFO]   Salmon quant:      $RUN_SALMON"
$RUN_BWA   && echo "[INFO]   bwa index:        $BWAIDX"
$RUN_SALMON && echo "[INFO]   Salmon index:      $SLMIDX"

# ensure all downstream directories are present
mkdir -p "$OUTDIR"/{data/{bamFiles,bigwig,counts},tmp/trim,logs}

# Export for GNU parallel subprocesses
export OUTDIR FQDIR BWAIDX SLMIDX THREADS RUN_SALMON RUN_BWA

# Check dependencies
REQUIRED_CMDS=(trim_galore parallel)
$RUN_SALMON && REQUIRED_CMDS+=(salmon)
$RUN_BWA && REQUIRED_CMDS+=(bwa samtools gatk bamCoverage)

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
printf "%s\n" "${samples[@]}" | parallel -j "$THREADS" '
	trim_galore --illumina --paired --basename {1} \
		"$FQDIR"/{1}_R1.fastq.gz "$FQDIR"/{1}_R2.fastq.gz \
		-o "$OUTDIR"/tmp/trim > "$OUTDIR"/logs/{1}_trim.log 2>&1
'

# Only progress once both reads of each samples has been trimmed
echo "[INFO] Trimming complete for all samples."


## Define processing function
process_sample() {
    local sample=$1

    echo "[INFO] Processing ${sample}..."

    local R1="$OUTDIR/tmp/trim/${sample}_val_1.fq.gz"
    local R2="$OUTDIR/tmp/trim/${sample}_val_2.fq.gz"

    # BWA MEM alignment
    bwa mem -t "$THREADS" "$BWAIDX" "$R1" "$R2" \
        -R "@RG\tID:${sample}\tSM:${sample}\tPL:ILLUMINA\tLB:${sample}\tPU:unit1" \
        2> "$OUTDIR/logs/${sample}_bwa.log" | \
        samtools view -@ "$THREADS" -b -F 3332 -q 20 - \
        > "$OUTDIR/data/bamFiles/${sample}.bwa.mem.bam"

    # BAM filtering and sorting
    gatk SortSam -I "$OUTDIR/data/bamFiles/${sample}.bwa.mem.bam" \
        -O "$OUTDIR/tmp/${sample}.nameSorted.bam" --SORT_ORDER queryname
    gatk FixMateInformation -I "$OUTDIR/tmp/${sample}.nameSorted.bam" \
        -O "$OUTDIR/tmp/${sample}.nameSorted.matefixed.bam"
    gatk SortSam -I "$OUTDIR/tmp/${sample}.nameSorted.matefixed.bam" \
        -O "$OUTDIR/tmp/${sample}.filtered.sorted.bam" --SORT_ORDER coordinate
    gatk MarkDuplicates -I "$OUTDIR/tmp/${sample}.filtered.sorted.bam" \
        -O "$OUTDIR/data/bamFiles/${sample}.filtered.sorted.rmDup.bam" \
        -M "$OUTDIR/data/bamFiles/${sample}.dupMetrics.txt" \
        --REMOVE_DUPLICATES true 
    samtools index "$OUTDIR/data/bamFiles/${sample}.filtered.sorted.rmDup.bam"

    # BigWig generation
    bamCoverage -b "$OUTDIR/data/bamFiles/${sample}.filtered.sorted.rmDup.bam" \
        -o "$OUTDIR/data/bigwig/${sample}.bw" \
        --normalizeUsing RPKM --numberOfProcessors "$THREADS" \
        > "$OUTDIR/logs/${sample}_bw.log" 2>&1

    # Optional Salmon quantification
    if $RUN_SALMON; then
        salmon quant -i "$SLMIDX" -l A -1 "$R1" -2 "$R2" -p "$THREADS" \
            -o "$OUTDIR/data/counts/${sample}_quant" > "$OUTDIR/logs/${sample}_salmon.log" 2>&1
    else
        echo "[INFO] Skipping Salmon quantification for ${sample}."
    fi

    echo "[INFO] Finished ${sample}."
}

export -f process_sample

## Run all samples in parallel
echo "[INFO] Running $( $RUN_BWA && echo 'bwa mem + bamCoverage' || echo '(no BWA)' ) $( $RUN_SALMON && echo '+ Salmon' || echo '(no Salmon)' )..."
printf "%s\n" "${samples[@]}" | parallel -j 4 process_sample {}

echo "[DONE] Pipeline completed successfully. All outputs in $OUTDIR"
