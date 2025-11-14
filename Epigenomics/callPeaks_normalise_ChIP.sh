#!/usr/bin/env bash
# ChIP peak calling and input normalisation
# Requires: macs2 and bedGraphToBigWig

set -euo pipefail

if [ "$#" -lt 5 ]; then
    echo -e "ChIP peak calling and input normalisation \n\
    Usage: $0 <sample_list.tsv> <output_dir> <genome_size> <chrom_size> [num_threads]\n\
    Note: <sample_list.tsv> should list one sample ID per line in the first column, and the input control for it in the second column (no extensions)\n\
    <genome_size> is the effective genome size [float], with shortcuts 'hs' (human), 'mm' (mouse), 'ce' (C. elegans), 'dm' (fruit fly). Common species include ateAlb (2.8e+8) and carPer (2.0e+9).\n\
    <chrom_size> is a file containing reference chromosome names and their sizes"
    exit 1
fi

SAMPLES=$1
OUTDIR=$2
GSIZE=$3
CHRSIZE=$4
THREADS=${5:-4}

# ensure all downstream directories are present
mkdir -p "$OUTDIR"/{data/{peaks,bigwig,bamFiles},logs}

# Export for GNU parallel subprocesses
export OUTDIR GSIZE CHRSIZE THREADS 

# Check dependencies
REQUIRED_CMDS=(macs2 bedGraphToBigWig)

for cmd in "${REQUIRED_CMDS[@]}"; do
    if ! command -v "$cmd" &>/dev/null; then
        echo "[ERROR] Required tool '$cmd' not found in PATH."
        exit 1
    fi
done

# Read sample IDs into an array
mapfile -t samples < <(grep -vE '^#|^[[:space:]]*$' "$SAMPLES")
echo "[INFO] Found ${#samples[@]} sample–input pairs to process."

## Define processing function
process_sample() {
    local sample=$1
    local input=$2

    # MACS peak calling
    macs2 callpeak -t "$OUTDIR/data/bamFiles/${sample}.filtered.sorted.rmDup.bam" \
    -c "$OUTDIR/data/bamFiles/${input}.filtered.sorted.rmDup.bam" -n "${sample}.filtered.sorted.rmDup.bam" -B \
    --outdir "$OUTDIR/data/peaks" -f BAM -g "$GSIZE" --call-summits

    # Subtract input
    macs2 bdgcmp -t "$OUTDIR/data/peaks/${sample}.filtered.sorted.rmDup.bam_treat_pileup.bdg" \
    -c "$OUTDIR/data/peaks/${sample}.filtered.sorted.rmDup.bamq()_control_lambda.bdg" \
    -m subtract --outdir "$OUTDIR/data/peaks" --o-prefix "${sample}"

    # Make bigwig from bedgraph
    sort -k1,1 -k2,2n "$OUTDIR/data/peaks/${sample}_subtract.bdg" > "$OUTDIR/data/peaks/${sample}_subtract.srtd.bdg"
    bedGraphToBigWig "$OUTDIR/data/peaks/${sample}_subtract.srtd.bdg" "$CHRSIZE" "$OUTDIR/data/bigwig/${sample}_noIP.bw"

    echo "[INFO] Finished ${sample}."
}

export -f process_sample

# Read both columns and pass them to GNU parallel
echo "[INFO] Running MACS2 peak caller and pileup normalisation..."
parallel --colsep '\t' -j "$THREADS" process_sample {1} {2} :::: "$SAMPLES"

echo "[DONE] Pipeline completed successfully. All outputs in $OUTDIR"
