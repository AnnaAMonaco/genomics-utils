#!/usr/bin/env bash
# Trim reads, align with STAR, generate BigWigs, and quantify with Salmon (optional)
# Requires: trim_galore, samtools >= 1.13, GNU parallel, STAR, salmon, deepTools
#   - <sample_list.txt> should list one sample ID per line (no extensions)
#   - Each sample must have paired FASTQs: ${sample}_R1.fastq.gz and ${sample}_R2.fastq.gz

set -euo pipefail

if [ "$#" -lt 5 ]; then
    echo -e "Trim reads, align with STAR, generate BigWigs, and quantify with Salmon (optional)\n\
    Usage: $0 <sample_list.txt> <fastq_dir> <output_dir> <star_index> <salmon_index> [num_threads] [--no-salmon]\n\
    Note: <sample_list.txt> should list one sample ID per line (no extensions)\n\
	To skip Salmon quantification, use --no-salmon."
    exit 1
fi

SAMPLES=$1
FQDIR=$2
OUTDIR=$3
STARIDX=$4
SLMIDX=$5
THREADS=${6:-4}
RUN_SALMON=true

# detect if salmon quant should be skipped
if [[ "${7:-}" == "--no-salmon" ]]; then
    RUN_SALMON=false
    echo "[INFO] Salmon quantification disabled (--no-salmon flag set)."
fi

# ensure all downstream directories are present
mkdir -p "$OUTDIR"/{data/{bamFiles,bigwig,counts},tmp/trim,logs}
# salmon=/home/amonaco/miniconda3/bin/salmon

# Export for GNU parallel subprocesses
export OUTDIR FQDIR STARIDX SLMIDX THREADS RUN_SALMON

# Check dependencies
REQUIRED_CMDS=(trim_galore samtools STAR salmon parallel bamCoverage)
$RUN_SALMON && REQUIRED_CMDS+=(salmon)
# $RUN_SALMON && REQUIRED_CMDS+=($salmon)

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


## 2. Define processing function
process_sample() {
    local sample=$1

    echo "[INFO] Processing ${sample}..."

    local R1="$OUTDIR/tmp/trim/${sample}_val_1.fq.gz"
    local R2="$OUTDIR/tmp/trim/${sample}_val_2.fq.gz"
    local STARprefix="$OUTDIR/data/bamFiles/${sample}"

    # STAR alignment
    STAR --runThreadN "$THREADS" --genomeDir "$STARIDX" --readFilesIn "$R1" "$R2" \
        --outFileNamePrefix "$STARprefix" --outSAMunmapped Within --outSAMtype BAM Unsorted \
        --readFilesCommand zcat --alignEndsType Local > "$OUTDIR/logs/${sample}_STAR.log" 2>&1

    # BAM filtering and sorting
    samtools view -@ "$THREADS" -b -F 2308 -q 20 "${STARprefix}Aligned.out.bam" | \
        samtools sort -@ "$THREADS" -n -T "$OUTDIR/tmp/${sample}" -o "$OUTDIR/tmp/${sample}.nameSorted.bam" -
    samtools fixmate -@ "$THREADS" "$OUTDIR/tmp/${sample}.nameSorted.bam" "$OUTDIR/tmp/${sample}.matefixed.bam"
    samtools sort -@ "$THREADS" -T "$OUTDIR/tmp/${sample}" \
    	-o "$OUTDIR/data/bamFiles/${sample}.filtered.sorted.bam" "$OUTDIR/tmp/${sample}.matefixed.bam"
    samtools index "$OUTDIR/data/bamFiles/${sample}.filtered.sorted.bam"

    # BigWig generation
    bamCoverage -b "$OUTDIR/data/bamFiles/${sample}.filtered.sorted.bam" -o "$OUTDIR/data/bigwig/${sample}.bw" \
        --normalizeUsing RPKM --numberOfProcessors "$THREADS" > "$OUTDIR/logs/${sample}_bw.log" 2>&1

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
echo "[INFO] Running STAR + bamCoverage $( $RUN_SALMON && echo '+ Salmon' || echo '(no Salmon)' )..."
printf "%s\n" "${samples[@]}" | parallel -j "$THREADS" process_sample {}

echo "[DONE] Pipeline completed successfully. All outputs in $OUTDIR"
