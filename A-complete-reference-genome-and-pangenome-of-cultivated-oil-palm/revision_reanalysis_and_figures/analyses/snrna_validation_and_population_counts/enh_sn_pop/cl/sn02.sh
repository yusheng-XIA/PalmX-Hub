#!/bin/bash
W=${CLUSTER_WORK}/enh_sn_pop
ulimit -v 120000000
export TMPDIR=$W/tmp OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8 NUMBA_NUM_THREADS=8 NUMBA_CACHE_DIR=$W/tmp
nice -n 10 python $W/cl/sn02_qc.py > $W/logs/sn02.log 2>&1
echo "exit $?" >> $W/logs/sn02.log
