#!/bin/bash
# usage: 11_queue.sh "LABEL1:QUERY1 LABEL2:QUERY2 ..."   (sequential on this node)
ulimit -v 199000000
W=${CLUSTER_WORK}/fig5hi_mask
for item in $1; do
  lab=${item%%:*}; q=${item#*:}
  if [[ -s $W/calls/$lab/${lab}_vs_Africa_hap2.snp.txt ]]; then continue; fi
  bash $W/scripts/10_call_dsnp_asm.sh $lab $q > $W/logs/$lab.out 2> $W/logs/$lab.err || echo "FAILED $lab" >> $W/logs/$lab.err
done
