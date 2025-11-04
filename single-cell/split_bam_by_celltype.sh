#!/usr/bin/env bash
# Split a 10X cellranger output BAM file into cell-type–specific BAMs using samtools
# Requires: samtools >= 1.13, GNU parallel, awk, sort

set -euo pipefail

if [ "$#" -lt 3 ]; then
    echo -e "Split a 10X cellranger output BAM file into cell-type–specific BAMs \
    \nUsage: $0 <barcode_to_celltype.txt> <possorted.bam> <output_dir> [num_threads]"
    exit 1
fi

ANNOTATION_FILE=$1
BAM_FILE=$2
OUTDIR=$3
THREADS=${4:-4}

mkdir -p "$OUTDIR"

# Check dependencies
for cmd in samtools parallel awk; do
    if ! command -v "$cmd" &>/dev/null; then
        echo "[ERROR] Required tool '$cmd' not found in PATH."
        exit 1
    fi
done

if [ ! -f "${BAM_FILE}.bai" ]; then
    echo "[INFO] Index not found; indexing BAM..."
    samtools index "$BAM_FILE"
fi

echo "[INFO] Collecting unique cell types..."
cell_types=($(awk 'NR>1 {print $2}' "$ANNOTATION_FILE" | sort | uniq))
echo "[INFO] Found ${#cell_types[@]} unique cell types."

# Prepare barcode list files
for ct in "${cell_types[@]}"; do
    safe_ct=$(echo "$ct" | tr ' /' '__')
    awk -v t="$ct" 'NR>1 && $2==t {print $1}' "$ANNOTATION_FILE" > "$OUTDIR/${safe_ct}_barcodes.txt"
done

# Extract BAMs using samtools
extract_bam() {
    local ct="$1"
    local safe_ct
    safe_ct=$(echo "$ct" | tr ' /' '__')

    local bc_file="$OUTDIR/${safe_ct}_barcodes.txt"
    local out_bam="$OUTDIR/${safe_ct}.bam"

    echo "[INFO] [$$] Extracting $ct ..."
    samtools view -@ 8 -b -D CB:"$bc_file" "$BAM_FILE" -o "$out_bam"
    samtools index "$out_bam"
    echo "[INFO] [$$] Finished $ct"
}

export -f extract_bam
export BAM_FILE OUTDIR

# Run in parallel
printf "%s\n" "${cell_types[@]}" | parallel -j "$THREADS" extract_bam {}

echo "[DONE] All cell-type–specific BAM files written to $OUTDIR"
