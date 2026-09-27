#!/bin/bash
# waits for each re-run call, then computes window coverage (20_coverage.py) and window dSNP counts (21_dsnp_windows.py)
W=${CLUSTER_WORK}/fig5hi_mask
PY=python3
ulimit -v 40000000
todo="nrly_final_hap2 nrly_final_hap1 pisifera_hap1 dura_hap1 dura_hap2 pisifera_hap2 nrly_old_hap2 nrly_old_hap1"
while [[ -n "$todo" ]]; do
  left=""
  for s in $todo; do
    P=$W/calls/$s/${s}_vs_Africa_hap2
    if grep -q '^done' $W/calls/$s/run.log 2>/dev/null; then
      $PY $W/scripts/20_coverage.py $s $P.sorted.paf $W/coverage/$s.tsv >> $W/logs/post.log 2>&1
      case $s in nrly_*) $PY $W/scripts/21_dsnp_windows.py $s $P.snp.txt $W/dsnp/$s > $W/logs/dsnp_$s.log 2>&1 ;; esac
      echo "post done $s $(date)" >> $W/logs/post.log
    else left="$left $s"; fi
  done
  todo=$(echo $left)
  [[ -n "$todo" ]] && sleep 300
done
$PY $W/scripts/24_panel_sites.py > $W/logs/panel_sites.log 2>&1
echo ALLPOST >> $W/logs/post.log
