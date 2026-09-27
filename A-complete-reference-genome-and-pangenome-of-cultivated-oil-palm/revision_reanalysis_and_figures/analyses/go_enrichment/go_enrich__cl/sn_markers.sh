#!/bin/bash
W=${CLUSTER_WORK}/go_enrich
export TMPDIR=$W/tmp OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4
ulimit -v 60000000
cd $W
nohup nice -n 10 python $W/cl/sn_markers.py > $W/logs/sn_markers.log 2>&1 &
echo started
