#!/bin/bash
# ADMIXTURE (no CV, default seed 43) on the 308-based LD-pruned structure set.  usage: admixK.sh K THREADS
W=${CLUSTER_WORK}/snp_repair_r308
IMG=${ANALYSIS_DIR}/05_GWAS/GWAS_tur/software/software/Reseq_genek.sif
mkdir -p $W/struct/admixK$1 && cd $W/struct/admixK$1
ln -sf ../all.bed all.bed; ln -sf ../all.bim all.bim; ln -sf ../all.fam all.fam
singularity exec -B ${DATA_ROOT} $IMG admixture -j$2 all.bed $1 > admix.$1.log 2>&1 && touch k$1.done || touch k$1.fail
