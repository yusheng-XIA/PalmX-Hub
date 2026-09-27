#!/bin/bash
# usage: 22_run_dsnp.sh LABEL SNP_TXT
W=${CLUSTER_WORK}/fig5hi_mask
ulimit -v 40000000
mkdir -p $W/dsnp
python3 $W/scripts/21_dsnp_windows.py $1 $2 $W/dsnp/$1 > $W/logs/dsnp_$1.log 2>&1
