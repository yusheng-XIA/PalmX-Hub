#!/bin/bash
W=${CLUSTER_WORK}/enh_sn_pop
ulimit -v 60000000
export TMPDIR=$W/tmp OMP_NUM_THREADS=4
nice -n 10 python $W/cl/sn03_stab.py > $W/logs/sn03.log 2>&1
echo "exit $?" >> $W/logs/sn03.log
