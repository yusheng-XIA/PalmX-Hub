# Haplotype-resolved genomes and a population pangenome of oil palm

Analysis code for the oil palm genome, pangenome and multi-omics study
(Xia, Y.S., Zeng, Q.G., Li, Y., Li, X.Y. *et al.*).

This folder keeps its original name, which is the link cited in the manuscript's Code availability statement.

## Data availability

| Database | Accession |
|----------|-----------|
| NCBI BioProject | [PRJNA1438016](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1438016) |
| NGDC BioProject | [PRJCA060109](https://ngdc.cncb.ac.cn/bioproject/browse/PRJCA060109) |
| PalmX-Hub | [https://circulargenome.com/palmxhub/](https://circulargenome.com/palmxhub/) (available upon publication) |

## Materials and naming

| Label | Material | Assemblies |
|-------|----------|------------|
| FL | Seedless interspecific hybrid Reyou-2 (*E. guineensis* × *E. oleifera*) | FL-Hap1 (predominantly *E. oleifera*), FL-Hap2 (predominantly *E. guineensis*; reference for all population analyses) |
| TN | Tenera hybrid Boke (dura × pisifera) | TN-Hap1, TN-Hap2 |
| TK | Deli dura | TK-Hap1, TK-Hap2 |
| NS | AVROS pisifera | NS-Hap1, NS-Hap2 |
| Nigerian | Nigerian *E. guineensis* accession | two haplotypes |
| *E. oleifera* | Independent *E. oleifera* accession | two haplotypes |
| — | 27 additional HiFi-only accessions | one chromosome-scale assembly each |

Some scripts keep the internal file and sample names used during the analysis:
`Africa_hap2` = FL-Hap2, `American_hap1` = FL-Hap1, `Dura`/`EG_dura`/`dura_hap*` = TK,
`Pisifera`/`EG_pisifera`/`pisifera_hap*` = NS, `nrly`/`EG_niriliya` = Nigerian accession,
`MZ4`/`meizhou4` = independent *E. oleifera* accession, `BK`/`Boke` = TN, `Phoenix` = date palm outgroup.
The 39 genome or haplotype assemblies represent 33 biological materials; "All38" = the 38 non-reference
assembly paths, "African35" = the 35 *E. guineensis* paths.

## Code structure

| Folder | Content | Figures / tables |
|--------|---------|------------------|
| `01_genome_assembly/` | Verkko (FL), hifiasm, HapHiC/Juicebox, RagTag, SubPhaser, CPhasing, assembly QC (BUSCO odb12, tidk telomeres, LAI, Merqury) | Fig. 1; Suppl. Data 1–3 |
| `02_genome_annotation/` | BRAKER3, EviAnn, GeneMark, Helixer, miniprot, HISAT2/StringTie/TransDecoder, EVM, eggNOG-mapper, Pfam support; allele pairing (GeneTribe, BLASTN, GMAP) | Fig. 3b, 4g; Suppl. Data 3 |
| `03_repeat_annotation/` | EDTA/RepeatMasker, LTR insertion time, TRF, BISER, TE density meta-gene profiles | Fig. 1, 4 |
| `04_population_genomics/` | Parabricks/GATK SNP calling and filter chain, ADMIXTURE, PCA, π, F<sub>ST</sub>, LD decay; `03_coding_zygosity/`: derived and homozygous coding variants per accession | Fig. 4; ED6; Suppl. Fig. 9 |
| `05_pangenome/` | OrthoFinder gene-family pangenome (39 assemblies), accumulation curves, Ks/Ka/Ks by class, minigraph-cactus graph, RGAs, WGDI karyotype, plantiSMASH BGCs | Fig. 4 |
| `06_structural_variants/` | SyRI + SVIM-asm + cuteSV per path, SURVIVOR, cross-path clustering and hybrid/dual-evidence catalogues, SV annotation, haplotype-pair SVs, SWave complex SVs, PanGenie genotyping | Fig. 5a–d; Suppl. Data 22 |
| `07_multi_omics/` | RNA-seq, WGCNA, DIA-NN/directLFQ proteomics, msconvert/xcms metabolomics with QC-RLSC, composite molecular scores, FAD2 qRT-PCR statistics, GO enrichment | Fig. 2; Suppl. Data 9–12, 25 |
| `08_snRNA_seq/` | Cell Ranger, Seurat/Harmony, Scrublet doublets, cluster stability, bulk scoring of snRNA markers | Suppl. Fig. 8; Suppl. Data 13 |
| `09_allele_specific_expression/` | Allele graphs (Minigraph, VG mpmap/surject), fragment-level allele counts, replicate-aware beta-binomial ASE test, temporal classes, FL DNA control, allelic Ka/Ks and equivalence tests | Fig. 3; ED4, ED5; Suppl. Data 14 |
| `10_phylogeny_and_gene_family/` | OrthoFinder/trimAl/IQ-TREE species tree, MCMCTree (control files and calibrations), CAFE5, copy-number (dosage) null model | Fig. 2a,b; Suppl. Fig. 2 |
| `11_gwas/` | Phenotype filtering, EMMAX SNP- and SV-GWAS, PC-adjusted and sensitivity scans, λ<sub>GC</sub>, SV intervals vs SNP signals, reciprocal conditional analysis and Wakefield fine-mapping, *SHELL* mutation test | Fig. 5; Suppl. Fig. 10; Suppl. Data 23 |
| `12_deleterious_variants_and_donor_design/` | Date-palm-polarised candidate dSVs and dSNPs, panel scope, dSV–dSNP co-localisation, *E. oleifera* outgroup sensitivity (nine definitions), neutral-class confounding controls, load matrices and exact dynamic-programming donor paths (coverage mask, W sweep, constrained designs, favourable loci) | Fig. 5h,i; ED8, ED10 |
| `revision_reanalysis_and_figures/` | Analyses re-run or added during revision (308-accession SNP filtering and population structure, GWAS, dSV/dSNP burden and donor paths, homozygous derived load, *OLE16a*, GO enrichment, ASE robustness, multi-omic axes, snRNA validation) and the plotting/figure-editing scripts for the final display items; see its own README | Fig. 1–5; ED1–10; Suppl. Fig. 1–10 |

## Two kinds of scripts

* **Analysis scripts** (Python/R; e.g. modules 04/03, 06, 09–12, 08/03–06, 07/05–08) are the programs used to
  produce the reported results, with local paths replaced by command-line arguments, environment variables or
  relative placeholders. Input formats are described in each script's header.
* **Workflow scripts** (`*.sh`) list the commands and parameters of the standard tools as reported in the Methods
  (`${sample}`, `${threads}` and file names are placeholders). They document the workflow and need to be adapted
  to local file names, job schedulers and software installations.

Software versions and all parameters are given in the Methods of the paper. Random seeds used in the analyses
are set inside the scripts (e.g. 42 for pangenome curves, 10086 for single-nucleus clustering,
20260630/20260924 for permutation and bootstrap analyses).

## Main software

Verkko 2.2.1, hifiasm 0.21.0, HapHiC, CPhasing 0.2.5, RagTag, BUSCO 6.0.0, Merqury 1.3, tidk 0.2.65,
LTR_retriever 3.0.2, EDTA 2.1.0, RepeatMasker 4.1.5, BRAKER3, EviAnn, GeneMark.hmm 4.68, Helixer, miniprot 0.12,
HISAT2 2.2.1, StringTie 2.2.1, TransDecoder 5.5.0, EvidenceModeler 2.1.0, eggNOG-mapper 2.1.13, HMMER 3.4,
Parabricks 4.3.0-1, GATK 4.2.0.0, bcftools, PLINK 1.9, ADMIXTURE 1.3.0, VCFtools 0.1.17, PopLDdecay 3.42,
OrthoFinder 3.1.1, trimAl 1.4, IQ-TREE 1.6.12, PAML 4.10.9 (MCMCTree), CAFE5 5.1.0, WGDI 0.6.5, plantiSMASH 1.0,
RGAugury, minigraph-cactus, vg, Minigraph 0.21, VG 1.67.0, minimap2 2.30, SyRI 1.7.1, SVIM-asm 1.0.3,
cuteSV 2.1.3, SURVIVOR 1.0.7, SWave 1.5.2, PanGenie 4.2.1, EMMAX (2012-02-10), KaKs_Calculator (GMYN),
GeneTribe 1.2.1, JCVI 1.5.11, BLAST+ 2.16.0, GMAP 2025-07-31, DIA-NN 2.2.0, directLFQ 0.3.3,
ProteoWizard 3.0.26121, xcms 4.8.0, DESeq2, WGCNA 1.72, clusterProfiler, Cell Ranger, Seurat 5.3.0,
Harmony 1.2.3, Scrublet 0.2.3, python-igraph 0.11.9.

Python scripts use Python ≥ 3.9 with numpy, pandas, scipy, statsmodels, pysam and Biopython.
