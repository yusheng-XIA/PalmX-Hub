#!/bin/bash
# original parameters: --window-pi 100000 --window-pi-step 50000 ; --weir-fst-pop x2 --fst-window-size 100000 --fst-window-step 50000
W=${CLUSTER_WORK}/snp_repair_rerun/fig4c
VT=vcftools
P=${1:-12}
G=$W/groups; mkdir -p $W/out/strict $W/out/present $W/logs/vt
J=$W/logs/vt_jobs_$(hostname).txt; : > $J
for set in strict present; do
 for f in $W/vcf/chr*.${set}.vcf.gz; do
  c=$(basename $f .${set}.vcf.gz)
  for p in 1 2 3 4; do
   o=$W/out/$set/${c}__K4_Pop${p}_100kb.pi
   [ -s $o.windowed.pi ] || echo "$VT --gzvcf $f --keep $G/K4_Pop$p.txt --window-pi 100000 --window-pi-step 50000 --out $o > $W/logs/vt/${set}_${c}_P$p.log 2>&1" >> $J
  done
  for pr in "1 2" "1 3" "1 4" "2 3" "2 4" "3 4"; do
   set -- $pr
   o=$W/out/$set/${c}__K4_Pop$1_K4_Pop$2_100kb_fst
   [ -s $o.windowed.weir.fst ] || echo "$VT --gzvcf $f --weir-fst-pop $G/K4_Pop$1.txt --weir-fst-pop $G/K4_Pop$2.txt --fst-window-size 100000 --fst-window-step 50000 --out $o > $W/logs/vt/${set}_${c}_F$1$2.log 2>&1" >> $J
  done
 done
done
[ -n "${EXCL:-}" ] && { grep -vxF -f "$EXCL" $J > $J.tmp; mv $J.tmp $J; }
wc -l $J
nice -n 10 xargs -P $P -I{} bash -c '{}' < $J
echo "vcftools done $(hostname) $(date '+%F %T')" >> $W/logs/x3_done.txt
