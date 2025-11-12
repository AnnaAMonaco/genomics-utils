# genomics-utils
A collection of helper scripts for fast and modular genomics analysis.

## RNA-seq
**Mapping paired-end reads:** `trim_map_quantify.sh` trims reads, aligns with `STAR`, generates BigWigs, and optionally quantifies with `salmon`.

## Epigenomics
**Mapping ChIP-seq or ATAC-seq:** `trim_map_wSalmon.sh` trims reads, aligns with `bwa mem`, generates BigWigs, and optionally quantifies with `salmon`. Quantification with `salmon` is useful in cases of allelic assignment for downstream allele-specific analysis.

## Single-cell
**Splitting BAM file in pseduobulk by cell type:** `split_bam_by_celltype.sh` takes the 10X `cellranger` output bam and a list of cell barcodes to cell types, and generates one bam file per cell type using `samtools`.
