#!/bin/bash
# usage: s3_admix.sh K THREADS  (parameters as 08_structure/admixture.sh: --cv=10; default seed 43)
W=${CLUSTER_WORK}/snp_repair_r308
IMG=${ANALYSIS_DIR}/05_GWAS/GWAS_tur/software/software/Reseq_genek.sif
mkdir -p $W/struct/admix; cd $W/struct/admix
ln -sf ../all.bed all.bed; ln -sf ../all.bim all.bim; ln -sf ../all.fam all.fam
singularity exec -B ${DATA_ROOT} $IMG admixture --cv=10 -j$2 all.bed $1 > admix.$1.log 2>&1
echo "K$1 done $(date -Is)" >> admix_done.txt
