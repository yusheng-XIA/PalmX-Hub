#!/bin/bash
# run a list of chromosomes on this node, P parallel (each pipeline ~3-4 threads)
W=${CLUSTER_WORK}/fig4c_final
P=$1; shift
for c in "$@"; do echo $c; done | nice -n 10 xargs -P $P -I{} bash -c "bash $W/cl/x2_extract.sh {} > $W/logs/x2_{}.log 2>&1 || echo FAIL {} >> $W/logs/x2_fail.txt"
echo "node $(hostname) done $(date '+%F %T')" >> $W/logs/x2_nodes_done.txt
