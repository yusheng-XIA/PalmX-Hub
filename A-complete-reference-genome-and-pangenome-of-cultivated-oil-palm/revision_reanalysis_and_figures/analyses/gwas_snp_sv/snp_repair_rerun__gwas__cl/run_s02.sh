#!/bin/bash
# usage: run_s02.sh "chr01B chr02B ..."  -- one worker per chromosome, one BLAS thread each, 12 GB cap per worker.
cd ${CLUSTER_WORK}/snp_repair_rerun/gwas/cl
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
PY=python
for c in $1; do
  ( ulimit -v 12000000; $PY s02_scan.py $c $2 > ../logs/s02${2}_$c.log 2>&1 ) &
done
wait; echo ALLDONE $(hostname) >> ../logs/s02${2}_done.txt
