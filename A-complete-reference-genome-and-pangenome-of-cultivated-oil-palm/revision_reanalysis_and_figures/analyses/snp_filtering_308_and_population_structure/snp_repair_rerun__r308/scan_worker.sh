#!/bin/bash
# EMMAX scan jobs (chromosome x model part), M0 and M1; waits for gwas/prep.done; exits when all jobs are claimed
W=${CLUSTER_WORK}/snp_repair_r308; PY=python
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 TMPDIR=$W/tmp
mkdir -p $W/gwas/jobs
until [ -e $W/gwas/prep.done ]; do [ -s $W/logs/stageB.fail ] && exit 1; sleep 60; done
cd $W/gwas/cl
for i in $(seq -w 1 16); do for k in 0 1 2 3; do
  j=chr${i}B_$k; mkdir $W/gwas/jobs/$j.claim 2>/dev/null || continue
  if ( ulimit -v 14000000; $PY s02_scan.py chr${i}B M0,M1 $k/4 > $W/gwas/logs/s02_$j.log 2>&1 ); then touch $W/gwas/jobs/$j.done
  else echo "$j $(hostname)" >> $W/gwas/jobs/fail.txt; fi
done; done
