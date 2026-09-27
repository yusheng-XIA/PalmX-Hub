#!/bin/bash
# EMMAX binary (same image and options as the published run) with intercept + SNP PC1-5, lauric acid, chr07B.
W=${CLUSTER_WORK}/enh_gwas; cd $W
IMAGE=${ANALYSIS_DIR}/05_GWAS/GWAS_tur/software/software/Reseq_genek.sif
T=${ANALYSIS_DIR}/05_GWAS/00_analysis/06_GWAS/yield/02_results/Nut_weight_g/03_GWAS_emmax/chr_split/GP_maf_chr07B
PH=${ANALYSIS_DIR}/22_answer_reviews/00_ms/03_V3/05_figure/07_GWAS_zero_extreme_excluded_20260723/phenotypes/quality/C12_0_Lauric_acid.txt
K=${ANALYSIS_DIR}/05_GWAS/00_analysis/06_GWAS/shared_genotype/GP_maf_allChr.hBN.kinf
awk -F'\t' '{printf "%s %s", $1, $1; for(i=2;i<=NF;i++) printf " %s", $i; print ""}' in/X_snp5pc.tsv > in/snp5pc.cov
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 SINGULARITYENV_OMP_NUM_THREADS=1 SINGULARITYENV_MKL_NUM_THREADS=1
mkdir -p tmp/bin
singularity exec -B ${DATA_ROOT} $IMAGE emmax-intel64 -v -d 10 -t $T -p $PH -k $K -c in/snp5pc.cov -o tmp/bin/lauric_pc_chr07B > tmp/bin/lauric_pc_chr07B.out 2>&1
echo DONE >> tmp/bin/lauric_pc_chr07B.out
