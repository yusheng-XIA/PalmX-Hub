#!/bin/bash
# after ADMIXTURE K = 4: new groups, then pi/FST and LD
W=${CLUSTER_WORK}/snp_repair_r308; R=${CLUSTER_WORK}/snp_repair_rerun/r308
until [ -e $W/struct/admix/k4.done ]; do [ -s $W/logs/stageB.fail ] && exit 1; sleep 60; done
python $R/k4groups.py > $W/logs/k4groups.log 2>&1 && bash $R/pifst.sh > $W/logs/pifst.log 2>&1 || echo "FAIL stageC $(date -Is)" >> $W/logs/stageC.fail
