#!/bin/bash
cd ${CLUSTER_WORK}/enh_3fgj
ulimit -v 64000000
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4 TMPDIR=${CLUSTER_WORK}/enh_3fgj/tmp
python s02_modes_robust.py 10 > s02.out 2>&1
