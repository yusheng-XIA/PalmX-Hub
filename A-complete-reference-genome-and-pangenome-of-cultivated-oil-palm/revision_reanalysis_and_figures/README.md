# Revision re-analyses and figure code

This directory contains the analysis and plotting code used for the revised version of
*Haplotype-resolved genomes and a population pangenome of oil palm*. It covers the analyses
that were re-run or added during revision and the scripts that generate or edit the final
Figures 1–5, Extended Data Figs. 1–10 and Supplementary Figs. 1–10. The code for genome assembly,
annotation, the gene-family pangenome, the graph pangenome and the original population-genomic
pipelines is in the sibling directories of this repository.

No data are included. Input tables, VCFs and matrices are described in the Methods and provided
as Supplementary Data and Source Data with the paper, or through the accessions listed under
Data availability.

## Layout

```
figures/                 plotting and figure-editing scripts, one directory per final display item
  common_style/          shared palette (palA), fonts, sizing and PDF/PNG check helpers
  Fig1 … Fig5/           main figures
  ED1 … ED10/            Extended Data figures (final numbering)
  SF1 … SF10/            Supplementary figures (final numbering)
  _post_edits/           scripts that apply in-figure text/label edits to the rendered PDFs
analyses/                analysis code, grouped by topic
  snp_filtering_308_and_population_structure/
  gwas_snp_sv/
  deleterious_burden_and_donor_path/
  homozygous_derived_load/
  ole16a/
  go_enrichment/
  ase_robustness/
  multiomics_axes/
  snrna_validation_and_population_counts/
```

Sub-directory names record where each group of scripts originated (for example `redraw__ED4`,
`snp_repair_rerun__ed6`). Several figures went through a base redraw followed by later edits;
all stages are kept so that each change can be traced.

## Figures: scripts → panels

Final numbering is used for the directories. Some scripts carry the internal numbering that was
in use when they were written; the correspondence is:

| Final | Internal name in scripts | Main scripts | Content |
|---|---|---|---|
| Fig. 1 | Fig1 | `plot_scripts__Fig1/panel_a…e`, `Figure1*.py`, `fig1a*`, `fig1b*`, `fig1c*` | map and crop shares, assembly metrics, FL-Hap2–EG11 alignment, inter-haplotype SVs, local ancestry |
| Fig. 2 | Fig2 | `plot_scripts__Fig2`, `Figure2_fix_script.py`, `fig2_ole16/` (`plot_panels.py`, `compose_fig2*.py`) | karyotype evolution, gene families, multi-omic axes, FL/TN fold changes, snRNA UMAP, *OLE16a* panels f–i |
| Fig. 3 | Fig3 | `plot_scripts__Fig3/panel_a…l`, `Figure3*.py`, `fig3k/`, `fig3l_key.py` | allelic architecture, ASE classes, cis/trans proxy classes, breeding targets |
| Fig. 4 | Fig4 | `plot_scripts__Fig4`, `Figure4*.py`, `fig4c_*`, `snp_repair_rerun__fig4/`, `draw_fig4ab_snp_repair.py` | ADMIXTURE (K = 3, 4, 8), PCA, π/F_ST networks and K = 8 heatmap, pangenome panels |
| Fig. 5 | Fig5 | `plot_scripts__Fig5`, `Figure5*.py`, `fig5c*`, `fig5hi*`, `edit_figure5hi_mask.py`, `h_labels/` | SV atlas, dSV–dSNP burden, donor mosaic (h, i) |
| ED1 | ED1 | `redraw__ED1/plot_ED1.py`, `ed1a_redraw/`, `ed1d_redraw/`, `redraw2__ED__src/` | study design, SubPhaser, read depth, junction reads, Pore-C map |
| ED2 | ED2 | `fabi__ED2/plot_ED2.py` | microsynteny of five PAV candidates |
| ED3 | ED3 | `redraw__ED3/plot_ED3.py`, `ed3a_upset_windows/` | ASE UpSet and robustness |
| ED4 | SF13 | `make_sf_new_3fgj.py` (analysis in `analyses/ase_robustness/enh_3fgj`) | Ka/Ks equivalence, expression-mode robustness |
| ED5 | SF12 | `make_sf12.py`, `sf12_data.py`, `s07_fig.py` | genome-wide allelic bias and DNA control |
| ED6 | ED4 | `redraw__ED4/plot_ED4.py`, `snp_repair_rerun__ed6/`, `r1__ed4/` | LD decay, π/F_ST, CV error, phylogeny, ancestry |
| ED7 | ED5 | `redraw__ED5/plot_ED5.py` | dSV distribution |
| ED8 | SF11 | `make_sf11.py`, `sf11ef_data.py`, `enhB_05_plot.py` | dSV definitions and confounding controls |
| ED9 | ED6 | `restructure__ed6_no_fg/plot_ED6_no_fg.py`, `snp_repair_rerun__ed9/`, `ed9c_snp_threshold/` | SHELL-region SNP/SV association |
| ED10 | SF10 | `fav233__SF10/make_sf10.py`, `plot_enh_C.py` | constrained donor designs and weight scan |
| SF1 | SF1 | `plot_SF1.py`, `redraw2__SF__SF1/plot_SF1.py` | Hi-C contact maps |
| SF2 | SF2 | `redraw__SF2/plot_SF2_p1.py`, `plot_SF2_p2.py`, `redraw__SF10/plot_SF10.py` | lipid gene dosage, gene trees, dosage null |
| SF3 | SF3 | `plot_SF3.py` | WGCNA |
| SF4 | SF4 | `plot_SF4.py`, `orphan_SF4/plot_SF4.py` | metabolome QC and metabolite heatmap |
| SF5 | SF5 | `plot_SF5.py`, `idea_A__integrate/plot_SF5_ole.py` | proteome and *OLE16* loci |
| SF6 | SF6 | `plot_SF6.py` | proteome and metabolome QC |
| SF7 | SF7 | `redraw__SF7/plot_SF7.py` | FAD2 qPCR |
| SF8 | SF8 | `plot_SF8.py`, `enh_sn_pop__SF8/` | snRNA-seq QC, markers, bulk validation |
| SF9 | SF14 | `snp_repair_rerun__sf9/…/plot_SF9_snp_repair.py`, `plot_SF9_original_G6.py` | *SHELL* alleles and homozygous derived coding load |
| SF10 | SF9 | `redraw__SF9/plot_SF9.py` | SV-GWAS |

## Analyses

| Directory | What it does | Paper |
|---|---|---|
| `snp_filtering_308_and_population_structure` | chunked hard filtering and genotype masking of the jointly genotyped VCFs restricted to the 308 accessions, PLINK sets, PCA, ADMIXTURE (K = 2–8, CV), LD decay (PopLDdecay), window π and F_ST for K = 3/4/8 groups | Fig. 4a–c; ED6; Results/Methods population genomics |
| `gwas_snp_sv` | EMMAX SNP and SV scans (two models), λGC, reporting intervals, SHELL fine-mapping and conditional tests, SV intervals far from SNP intervals, ED9 input tables | ED9; SF10; SD23 |
| `deleterious_burden_and_donor_path` | dSV/dSNP identification and re-identification, confounding controls, coverage-masked donor dynamic programming, W/P scans, constrained designs, favourable-locus derivation (233 loci) | Fig. 5e–i; ED7, ED8, ED10; SD24 |
| `homozygous_derived_load` | *SHELL* coding-allele genotyping and derived coding load in commercial versus other source groups | SF9 |
| `ole16a` | oleosin phylogeny, peptide-level quantification, RNA/ASE/snRNA summaries, public RNA-seq quantification, upstream-region comparison across 39 haplotypes | Fig. 2f–i; SF5d–j |
| `go_enrichment` | single-nucleus marker calling and clusterProfiler enrichment for clusters and WGCNA modules | SD10, SD13 |
| `ase_robustness` | Ka/Ks equivalence tests, thinning/replicate robustness, genome-wide allelic bias with DNA control | ED3–ED5 |
| `multiomics_axes` | multi-omic axis statistics, divergence-time/CAFE reruns | Fig. 2c–d; SF2 |
| `snrna_validation_and_population_counts` | bulk validation of snRNA clusters; allele-frequency contrasts between source groups | SF8i–k; SD16 |

## Paths and environment

Absolute paths were replaced by variables. Set them before running:

| Variable | Meaning |
|---|---|
| `WORK_DIR` | local working directory holding inputs/outputs of the figure scripts |
| `LOCAL_INPUT_DIR` | local directory with author-provided source tables |
| `CLUSTER_WORK`, `CLUSTER_HOME` | cluster working directories for the re-analyses |
| `DATA_DIR`, `DATA_DIR2`, `DATA_DIR3`, `DATA_ROOT` | root directories of the project data (assemblies, annotations, omics tables) |
| `ANALYSIS_DIR`, `GWAS_DIR` | project analysis and GWAS directories |
| `JOINT_SNP_DIR` | directory with the jointly genotyped per-chromosome VCFs and their filtered versions |
| `SCRATCH` | fast local scratch (tmpfs) |
| `LOGIN_HOST`, `COMPUTE_HOST` | cluster host names used in job-launch helpers |

Some comments are in Chinese; they describe the same steps as the code.

### Software

Python ≥ 3.9 with the packages in `requirements.txt`. R ≥ 4.3 with clusterProfiler 4.20.0,
limma 3.66.0, WGCNA, DESeq2, ggplot2, ggtree, treeio, ape, patchwork, cowplot, dplyr, tidyr,
readr, scales, maps, Matrix and GO.db.

Command-line tools called by the scripts (versions as in Methods): bcftools/htslib (tabix, bgzip),
VCFtools, PLINK 1.9, ADMIXTURE 1.3.0, EMMAX (emmax-kin, emmax), PopLDdecay, GATK 4.2,
minimap2 (paftools), samtools, bedtools, SyRI, SWave, IQ-TREE 2, MAFFT v7.525, PAML 4.10.10 (codeml),
wgd v2, CD-HIT, DIAMOND, BLAST+, HMMER, HISAT2, Salmon, CAFE5 and Singularity (for containerized tools).

## License

MIT; see `LICENSE` at the repository root. `figures/ED1/redraw2__ED__src/ED1e/juicerbox/HapHiC_plot.py` is
third-party code from HapHiC (Copyright (c) 2023, Xiaofei Zeng), redistributed under its BSD 3-Clause License
(`LICENSE_HapHiC` in the same directory).
