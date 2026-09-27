#!/bin/bash
# one worker: claims chunks from the shared list (atomic mkdir) until none left
W=${CLUSTER_WORK}/snp_repair_rerun/r308
export TMPDIR=$W/tmp
while IFS=$'\t' read -r raw reg out; do
  [ -s $out.done ] && continue
  mkdir -p $(dirname $out); mkdir $out.claim 2>/dev/null || continue
  s=$(date +%s)
  if bash $W/filter308.sh $raw $reg $out > $out.log 2>&1; then echo -e "$reg\t$(hostname)\t$(( $(date +%s)-s ))" >> $W/logs/times.tsv
  else echo -e "FAIL\t$reg\t$(hostname)" >> $W/logs/fail.tsv; fi
done < $W/chunks.tsv
