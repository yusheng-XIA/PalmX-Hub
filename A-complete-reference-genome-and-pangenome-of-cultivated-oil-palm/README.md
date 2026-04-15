# A complete reference genome and pangenome of cultivated oil palm

## Overview

Analysis code for the oil palm pangenome paper.

Xia, Y.S., Li, Y., Zeng, Q.G., Li, X.Y. et al.

## Data Availability

| Database | Accession |
|----------|-----------|
| NCBI BioProject | [PRJNA1438016](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1438016) |
| NGDC BioProject | [PRJCA060109](https://ngdc.cncb.ac.cn/bioproject/browse/PRJCA060109) |
| PalmX-Hub | [https://circulargenome.com/PalmX](https://circulargenome.com/PalmX) |

## Code Structure

Scripts are organized by analysis module:

- `01_genome_assembly/` — Genome assembly, scaffolding and quality evaluation
- `02_genome_annotation/` — Gene prediction, functional annotation and allelic gene identification
- `03_repeat_annotation/` — TE, tandem repeat and segmental duplication annotation
- `04_population_genomics/` — SNP calling, population structure and diversity analysis
- `05_pangenome/` — Pan-gene clustering, graph pangenome, RGA and gene family evolution
- `06_structural_variants/` — SV detection, merging, annotation and TE-mediated rearrangement analysis
- `07_multi_omics/` — RNA-seq, DIA proteomics, metabolomics, WGCNA and allele-specific expression
- `08_snRNA_seq/` — Single-nucleus RNA-seq processing and clustering

All scripts use generic variable names (e.g., `${sample}`, `${reference}`, `${threads}`) and should be adapted to your specific file paths and computing environment. Software versions used in this study are described in the Methods section of the paper.