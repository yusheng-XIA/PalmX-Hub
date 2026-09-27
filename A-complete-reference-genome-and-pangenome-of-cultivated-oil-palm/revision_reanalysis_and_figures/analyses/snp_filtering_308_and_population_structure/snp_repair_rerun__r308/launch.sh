#!/bin/bash
# usage: launch.sh NWORKERS   (run on a compute node)
W=${CLUSTER_WORK}/snp_repair_rerun/r308
for i in $(seq 1 $1); do nohup bash $W/worker.sh > /dev/null 2>&1 < /dev/null & done
