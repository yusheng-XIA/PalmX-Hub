#!/bin/bash
N=${CLUSTER_WORK}/snp_repair_rerun/fig4c
VT=vcftools; G=$N/groups
J=$N/logs/jobs_5chr.txt; : > $J; mkdir -p $N/logs/vt
for c in chr04B chr05B chr07B chr10B chr16B; do for set in present strict; do f=$N/vcf/$c.$set.vcf.gz
  for p in 1 2 3 4; do echo "$VT --gzvcf $f --keep $G/K4_Pop$p.txt --window-pi 100000 --window-pi-step 50000 --out $N/out/$set/${c}__K4_Pop${p}_100kb.pi > $N/logs/vt/${set}_${c}_P$p.log 2>&1" >> $J; done
  for pr in "1 2" "1 3" "1 4" "2 3" "2 4" "3 4"; do set -- $pr; echo "$VT --gzvcf $f --weir-fst-pop $G/K4_Pop$1.txt --weir-fst-pop $G/K4_Pop$2.txt --fst-window-size 100000 --fst-window-step 50000 --out $N/out/$set/${c}__K4_Pop$1_K4_Pop$2_100kb_fst > $N/logs/vt/${set}_${c}_F$1$2.log 2>&1" >> $J; done
  echo "python3 $N/cl/pi_ma.py $f $G $N/out/pima_$set/$c > $N/logs/pima_${set}_$c.log 2>&1" >> $J
done; done
nice -n 10 xargs -P 24 -I{} bash -c "{}" < $J
echo "5chr done $(date -Is)" >> $N/logs/x34_done.txt
