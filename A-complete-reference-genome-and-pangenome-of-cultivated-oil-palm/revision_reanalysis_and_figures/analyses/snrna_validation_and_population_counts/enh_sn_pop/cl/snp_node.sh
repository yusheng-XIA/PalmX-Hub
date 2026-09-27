#!/bin/bash
W=${CLUSTER_WORK}/enh_sn_pop
export TMPDIR=$W/tmp
P=$1; shift
for c in "$@"; do echo $c; done | nice -n 10 xargs -P $P -I{} bash -c "bash $W/cl/snp1.sh {} > $W/logs/snp_{}.log 2>&1 || echo FAIL {} >> $W/logs/snp_fail.txt"
echo "node $(hostname) done $(date '+%F %T')" >> $W/logs/snp_nodes_done.txt
