#!/bin/bash
# Genome-wide steps on the 308-based sets (run on ${COMPUTE_HOST} with nohup). Waits for the 16 stage-A chromosomes.
set -euo pipefail
W=${CLUSTER_WORK}/snp_repair_r308; R=${CLUSTER_WORK}/snp_repair_rerun/r308
IMG=${ANALYSIS_DIR}/05_GWAS/GWAS_tur/software/software/Reseq_genek.sif
P="singularity exec -B ${DATA_ROOT} $IMG plink"; PY=python
export TMPDIR=$W/tmp OPENBLAS_NUM_THREADS=4 OMP_NUM_THREADS=4
trap 'echo "FAIL line $LINENO $(date -Is)" >> $W/logs/stageB.fail' ERR
st() { echo "$1 $(date -Is)" >> $W/logs/stageB.steps; }
rm -f $W/ld_k4new/tables_smooth/* $W/fig4c_k4new/pi_fst_genomewide.tsv
# start as soon as the PLINK sets, Fig. 4c SNP sets and G6 genotypes exist (SD16 counting may still run)
until [ $(ls $W/s308/*.s2.done 2>/dev/null | wc -l) = 16 ] && [ $(ls $W/g6/geno/*.tsv.gz 2>/dev/null | wc -l) = 16 ]; do [ -s $W/logs/stageA_fail.txt ] && { echo stageA failed >> $W/logs/stageB.fail; exit 1; }; sleep 60; done
st plink_sets_complete
( cd $W/g6 && bash cl/run_g6.sh; st g6 ) > $W/logs/g6.log 2>&1 &
( until [ $(ls $W/logs/stageA_*.done 2>/dev/null | wc -l) = 16 ]; do sleep 60; done; cd $W/sd16 && python3 cl/sum.py snp > $W/sd16/sd16_counts.tsv; st sd16 ) > $W/logs/sd16.log 2>&1 &
# structure set, LD pruning, PCA
bash $W/cl/s3_struct.sh > $W/logs/s3_struct.log 2>&1; bash $W/cl/s3_trace.sh > $W/logs/s3_trace.log 2>&1 || true; st struct
# ADMIXTURE K = 4 (no CV) for the group assignments
( mkdir -p $W/struct/admix && cd $W/struct/admix && ln -sf ../all.bed all.bed && ln -sf ../all.bim all.bim && ln -sf ../all.fam all.fam \
  && singularity exec -B ${DATA_ROOT} $IMG admixture -j40 all.bed 4 > admix.4.log 2>&1 && touch k4.done && st admixture_k4 ) &
# GWAS genotypes and kinship
mkdir -p $W/gwas/geno $W/gwas/logs; cd $W/gwas/geno
: > merge.txt; for i in $(seq -w 1 16); do echo $W/s308/chr${i}B.gwas >> merge.txt; done
$P --merge-list merge.txt --allow-extra-chr --keep-allele-order --make-bed --out GP_new > /dev/null
awk '{print $1, $2}' ${ANALYSIS_DIR}/05_GWAS/00_analysis/06_GWAS/shared_genotype/GP_maf_allChr.fam > order.txt
$P --bfile GP_new --indiv-sort f order.txt --allow-extra-chr --keep-allele-order --make-bed --out GP_new_o > /dev/null
$P --bfile GP_new_o --allow-extra-chr --keep-allele-order --recode 12 transpose --output-missing-genotype 0 --out GP_new_o > /dev/null
singularity exec -B ${DATA_ROOT} $IMG /opt/emmax-20120210/emmax-kin-intel64 -v -x -d 10 -o GP_new_o.aBN.kinf GP_new_o > kin.log 2>&1
rm -f GP_new_o.tped
cd $W/gwas/cl && $PY s00_prep.py > ../logs/s00.log 2>&1
$PY -c "import numpy as np;K=np.loadtxt('$W/gwas/geno/GP_new_o.aBN.kinf');assert np.isfinite(K).all() and K.shape==(308,308);print('kin ok')" >> $W/gwas/logs/s00.log
touch $W/gwas/prep.done; st gwas_prep
# wait for the 64 scan jobs (scan_worker.sh on ${COMPUTE_HOST} + hnodes), then the summaries
until [ $(ls $W/gwas/jobs/*.done 2>/dev/null | wc -l) = 64 ]; do [ -s $W/gwas/jobs/fail.txt ] && { echo scan failed >> $W/logs/stageB.fail; exit 1; }; sleep 60; done
st scans
bash $R/post.sh > $W/logs/post.log 2>&1; st post
wait; st all
