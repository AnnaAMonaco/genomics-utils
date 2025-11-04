#!/usr/bin/env bash
# Trim reads, align with STAR, generate BigWigs, and quantify with Salmon
# Requires: trim_galore, samtools >= 1.13, GNU parallel, STAR, salmon, deepTools
#   - <sample_list.txt> should list one sample ID per line (no extensions)
#   - Each sample must have paired FASTQs: ${sample}_R1.fastq.gz and ${sample}_R2.fastq.gz
#   - Each step runs in parallel per sample; trimming must finish before alignment/quantification

set -euo pipefail

if [ "$#" -lt 5 ]; then
    echo -e "Trim reads, align with STAR, generate BigWigs, and quantify with Salmon \n\
    Usage: $0 <sample_list.txt> <fastq_dir> <output_dir> <star_index> <salmon_index> [num_threads] \n\
    Note: <sample_list.txt> should list one sample ID per line (no extensions)"
    exit 1
fi

SAMPLES=$1
FQDIR=$2
OUTDIR=$3
STARIDX=$4
SLMIDX=$5
THREADS=${6:-4}

# ensure all downstream directories are present
mkdir -p "$OUTDIR"/{data/{bamFiles,bigwig,counts},tmp/trim,logs}
export PATH="/home/amonaco/miniconda3/bin:$PATH"

# Check dependencies
for cmd in trim_galore samtools STAR salmon parallel bamCoverage; do
    if ! command -v "$cmd" &>/dev/null; then
        echo "[ERROR] Required tool '$cmd' not found in PATH."
        exit 1
    fi
done

# Read sample IDs into an array
mapfile -t samples < <(grep -vE '^#|^[[:space:]]*$' "$SAMPLES")
echo "[INFO] Found ${#samples[@]} samples to process."


## 1. Trim reads
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

    # Salmon quantification
    salmon quant -i "$SLMIDX" -l A -1 "$R1" -2 "$R2" -p "$THREADS" \
    	-o "$OUTDIR/data/counts/${sample}_quant" > "$OUTDIR/logs/${sample}_salmon.log" 2>&1

    echo "[INFO] Finished ${sample}."
}

export -f process_sample
export OUTDIR FQDIR STARIDX SLMIDX THREADS


## 3. Run all samples in parallel
echo "[INFO] Running STAR + bamCoverage + Salmon in parallel..."
printf "%s\n" "${samples[@]}" | parallel -j "$THREADS" process_sample {}

echo "[DONE] Pipeline completed successfully. All outputs in $OUTDIR"
