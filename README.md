# CpG-to-gene association benchmark across multi-omic datasets

**Work in Progress** — This repository is under active development as part of a
postgraduate Master's thesis project. The pipeline, results, and documentation are not yet
final. Interfaces, file formats, and parameters may change without notice until
the first release.

## Overview

This project benchmarks five methods for associating differentially methylated
CpG sites (DM) to differentially expressed genes (DE) across multi-omic
datasets combining bisulfite sequencing (BS-Seq) and RNA-Seq. The benchmark
evaluates how much different CpG-to-gene association strategies affect
downstream biological conclusions — a question of practical relevance when
integrating epigenomic and transcriptomic layers.

## Datasets

Four datasets are used, each chosen to probe a different regime of biological
signal strength:

| ID   | Source              | Contrast                  | Signal size    | Role in benchmark              |
|------|---------------------|---------------------------|----------------|--------------------------------|
| D0   | Synthetic           | Simulated case vs control | Controlled     | Positive control, ground truth |
| D1   | Pacis 2015          | TB-infected vs control    | Strong         | Well-defined biological effect |
| D2   | Mendizabal          | NeuN vs Olig2 (positive)  | Very strong    | Cell-type contrast (validation)|
| D3   | Fluzone vaccine     | Day 28 vs Day 0           | Sparse         | Applicability-domain limit     |
|      | (GSE263782)         |                           |                |                                |

**D0** provides ground truth for validation because both DM and DE signals are
generated with known CpG-to-gene assignments.

## CpG-to-gene association methods evaluated

Five methods are compared, spanning different assumptions about how CpG sites
should be associated to genes:

- **M1** — Nearest-TSS, no hierarchy
- **M2** — Nearest-TSS, with hierarchy
- **M3** — Broad promoter subset of M2 (-2000, +200)bp
- **M4** — Window range: ±10kb
- **M5** — Window range with basal-plus-extension domain (±5kb / +1kb basal, 1Mb extension)

## Analysis pipeline

Each dataset is processed by 6 scripts (7 for D0, which includes a data
generator):

1. **RNA-Seq differential expression** (DESeq2 with paired covariate)
2. **BS-Seq differential methylation** (methylKit for D1/D2; limma-on-M-values
   for D3)
3. **CpG-to-gene annotation** (methods M1–M5)
4. **Integration** — three complementary pipelines:
   - PIP0: gene-list intersection (Upset plot)
   - PIP1: per-CpG Spearman correlation (Delta and eQTM-style residuals)
   - PIP2: quadrant enrichment (concordance of ΔM sign and log2FC sign)
5. **Additional statistics** — hypergeometric tests for direction consistency,
   magnitude correlation, multi-DM enrichment, sites-vs-log2FC correlation
6. **Genomic region analysis** — re-annotation with ChIPseeker to compare
   method outputs at the region level (promoter, exon, intron, UTR, intergenic)

## Requirements

The pipeline is R-based (≥ 4.2). Key packages used across scripts:

- `DESeq2`, `tximport`, `matrixStats` — RNA-Seq
- `methylKit`, `limma`, `EpiDISH` — BS-Seq
- `genomation`, `ChIPseeker`, `rGREAT`, `GenomicFeatures`, `GenomicRanges` — annotation
- `TxDb.Hsapiens.UCSC.hg38.knownGene`, `org.Hs.eg.db` — gene models
- `rtracklayer`, `karyoploteR` — coordinate handling and visualization
- `dplyr`, `tidyr`, `UpSetR`, `cowplot`, `ggplot2`, `pheatmap` — data manipulation and plotting
- `openxlsx`, `writexl`, `data.table` — I/O

Reference assembly is **GRCh38** for all datasets. D1 and D2 raw data are
lifted from hg19 to hg38 via `rtracklayer::liftOver`; D3 is already in hg38.

## Methodological notes and known limitations

Because this is a work in progress, some methodological choices are
still evolving and will be corrected and changed before the write-up of the final thesis:
