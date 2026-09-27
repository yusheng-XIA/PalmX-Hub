#!/bin/bash
# K = 3 and K = 8: groups from the 308-based ADMIXTURE, then pi/FST (waits for the ADMIXTURE runs and the Fig. 4c SNP set)
W=${CLUSTER_WORK}/snp_repair_r308; R=${CLUSTER_WORK}/snp_repair_rerun/r308
PY=python
K=$1; P=$2
until [ -e $W/struct/admixK$K/k$K.done ]; do [ -e $W/struct/admixK$K/k$K.fail ] && exit 1; sleep 60; done
$PY $R/kgroups.py $K > $W/logs/kgroups_$K.log 2>&1 && bash $R/pifst_K.sh $K $P > $W/logs/pifst_K$K.log 2>&1 || echo "FAIL K$K $(date -Is)" >> $W/logs/stageD.fail
