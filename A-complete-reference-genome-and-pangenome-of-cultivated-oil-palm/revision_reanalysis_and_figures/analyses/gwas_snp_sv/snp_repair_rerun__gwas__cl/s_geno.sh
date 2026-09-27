#!/bin/bash
# merge 16 chromosome GWAS sets; tped for emmax-kin; hBN kinship (-v -h -d 10, as published)
set -euo pipefail
W=${CLUSTER_WORK}/snp_repair_rerun
IMG=${ANALYSIS_DIR}/05_GWAS/GWAS_tur/software/software/Reseq_genek.sif
P="singularity exec -B ${DATA_ROOT} $IMG plink"
mkdir -p $W/gwas/geno; cd $W/gwas/geno
: > merge.txt; for i in $(seq -w 1 16); do echo $W/s308/chr${i}B.gwas >> merge.txt; done
$P --merge-list merge.txt --allow-extra-chr --keep-allele-order --make-bed --out GP_new
$P --bfile GP_new --allow-extra-chr --keep-allele-order --recode 12 transpose --output-missing-genotype 0 --out GP_new
singularity exec -B ${DATA_ROOT} $IMG emmax-kin-intel64 -v -h -d 10 GP_new > kin.log 2>&1
rm -f GP_new.tped
echo DONE > geno.done
