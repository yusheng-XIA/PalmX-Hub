#!/bin/bash
cd ${CLUSTER_WORK}/enh_3fgj
ulimit -v 64000000
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4 TMPDIR=${CLUSTER_WORK}/enh_3fgj/tmp
PY=python
nohup $PY s01_kaks_equiv.py > s01.out 2>&1 &
nohup $PY s02_modes_robust.py ${1:-10} > s02.out 2>&1 &
wait
