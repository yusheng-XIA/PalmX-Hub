#!/bin/bash
R=${CLUSTER_WORK}/snp_repair_rerun/r308; W=${CLUSTER_WORK}/snp_repair_r308
K=$1; P=$2
python $R/kgroups_prev.py $K > $W/logs/kgroups_prev_$K.log 2>&1 && bash $R/pifst_K.sh $K $P _prevQ > $W/logs/pifst_K${K}_prevQ.log 2>&1 || echo "FAIL prevQ K$K" >> $W/logs/stageD.fail
