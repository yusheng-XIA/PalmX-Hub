#!/bin/bash
# merge 16 chromosome structure sets; LD prune (50 10 0.2); PCA 10 (PLINK 1.9, as 07_PCA/run_pca.sh); write all.bed for ADMIXTURE
set -euo pipefail
W=${CLUSTER_WORK}/snp_repair_r308
IMG=${ANALYSIS_DIR}/05_GWAS/GWAS_tur/software/software/Reseq_genek.sif
P="singularity exec -B ${DATA_ROOT} $IMG plink"
mkdir -p $W/struct; cd $W/struct
: > merge.txt; for i in $(seq -w 1 16); do echo $W/s308/chr${i}B.struct >> merge.txt; done
$P --merge-list merge.txt --allow-extra-chr --keep-allele-order --make-bed --out all.missing_maf
$P --bfile all.missing_maf --indep-pairwise 50 10 0.2 --allow-extra-chr --out tmp.ld
$P --bfile all.missing_maf --extract tmp.ld.prune.in --allow-extra-chr --keep-allele-order --make-bed --out all.LDfilter
$P --bfile all.LDfilter --pca 10 --allow-extra-chr --out PCA_10
# ADMIXTURE needs integer chromosome codes
awk "BEGIN{OFS=\"\t\"}{c=\$1; sub(/^chr/,\"\",c); sub(/B$/,\"\",c); \$1=c+0; print}" all.LDfilter.bim > all.bim
ln -sf all.LDfilter.bed all.bed; ln -sf all.LDfilter.fam all.fam
echo DONE > s3.done
