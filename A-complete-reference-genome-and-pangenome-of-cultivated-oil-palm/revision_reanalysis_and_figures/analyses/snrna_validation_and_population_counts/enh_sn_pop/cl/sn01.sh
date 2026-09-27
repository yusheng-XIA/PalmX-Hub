#!/bin/bash
W=${CLUSTER_WORK}/enh_sn_pop
ulimit -v 120000000
export TMPDIR=$W/tmp OMP_NUM_THREADS=4
cd $W/sn
nice -n 10 Rscript $W/cl/sn01_export.R > $W/logs/sn01.log 2>&1
echo "exit $?" >> $W/logs/sn01.log
