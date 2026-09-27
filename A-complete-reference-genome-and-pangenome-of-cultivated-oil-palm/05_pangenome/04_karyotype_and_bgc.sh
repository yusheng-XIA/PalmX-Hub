#!/bin/bash
# Ancestral palm karyotype (WGDI v0.6.5) and biosynthetic gene clusters (plantiSMASH v1.0)
set -euo pipefail
threads=32

# ---- 1. Ancestral karyotype ----------------------------------------------------------------------------------
# APK (n = 5): Phoenix dactylifera chromosomes 14, 2, 3, 8 and 6 (7,294 genes), identified by self-comparison.
# BLASTP of each palm proteome against the APK proteins (E <= 1e-5, up to 20 targets)
blastp -query ${sp}.pep.fa -db APK.pep.fa -evalue 1e-5 -max_target_seqs 20 -outfmt 6 -num_threads ${threads} \
    > ${sp}_vs_APK.blast
# WGDI collinearity (>= 5 gene pairs, P <= 0.2); [collinearity] section of ${sp}.conf: mg = 5, pvalue = 0.2
wgdi -icl ${sp}.conf
# Assign blocks to APK / ACPK (n = 10) chromosomes, merge adjacent blocks of the same ancestral chromosome
wgdi -km ${sp}.conf        # karyotype_mapping
wgdi -k  ${sp}.conf        # karyotype
# fissions = B - 10, fusions/joining = B - n (B continuous ancestral blocks, n haploid chromosome number)

# ---- 2. Biosynthetic gene clusters --------------------------------------------------------------------------
# 34 assemblies: both haplotypes of FL and TN; one assembly each of TK, NS, Nigerian and the 27 HiFi-only materials.
# BGCs of the independent E. oleifera accession transferred from FL-Hap1 with Liftoff.
for g in $(cat bgc_assemblies_34.txt); do
    plantismash --taxon plants --cpus ${threads} --outputfolder bgc_${g} ${g}.gbk   # GenBank with gene models
done
liftoff -g FL_Hap1.bgc_genes.gff3 -o Eoleifera.bgc_genes.gff3 Eoleifera.fa FL_Hap1.fa
# BGCs grouped into 52 pan-BGC families by predicted type and gene composition; two-haplotype materials are
# consolidated per material with copy number = max over haplotypes (814 material-level loci; core 11, shell 31,
# material-specific 10).
