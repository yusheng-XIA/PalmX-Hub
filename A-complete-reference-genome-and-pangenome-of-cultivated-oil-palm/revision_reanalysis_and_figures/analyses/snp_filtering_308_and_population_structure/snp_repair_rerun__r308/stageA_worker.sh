#!/bin/bash
# claims chromosomes whose 308-filter chunks are all done; exits when every chromosome is claimed
R=${CLUSTER_WORK}/snp_repair_rerun/r308; W=${CLUSTER_WORK}/snp_repair_r308
mkdir -p $W/logs
CH=$(cut -f2 $R/chunks.tsv | cut -d: -f1 | uniq)
while :; do
  left=0
  for c in $CH; do
    [ -e $W/logs/stageA_$c.claim ] && continue
    left=1
    n=$(awk -F'\t' -v c="$c" 'index($2, c":")==1' $R/chunks.tsv | wc -l); d=$(ls $R/filter/$c/*.done 2>/dev/null | wc -l)
    [ "$n" -gt 0 ] && [ "$n" = "$d" ] || continue
    mkdir $W/logs/stageA_$c.claim 2>/dev/null || continue
    bash $R/stageA_chr.sh $c > $W/logs/stageA_$c.log 2>&1 || echo "FAIL $c" >> $W/logs/stageA_fail.txt
  done
  [ $left = 0 ] && exit 0
  [ -s $R/logs/fail.tsv ] && exit 1
  sleep 60
done
