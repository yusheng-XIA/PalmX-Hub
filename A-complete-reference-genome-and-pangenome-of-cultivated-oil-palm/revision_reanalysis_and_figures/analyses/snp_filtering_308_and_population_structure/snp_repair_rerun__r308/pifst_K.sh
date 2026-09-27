#!/bin/bash
# pi (missing-aware per-site estimator, 100-kb windows / 50-kb step) and pairwise W&C mean FST (vcftools windows) for the
# K = 3 or K = 8 groups on the 308-based Fig. 4c SNP set (same files, windows and estimators as K = 4).  usage: pifst_K.sh K NPAR
set -euo pipefail
K=$1; P=$2; T=${3:-}   # T: optional suffix (e.g. _prevQ) for groups_K$K$T / fig4c_K$K$T
W=${CLUSTER_WORK}/snp_repair_r308; R=${CLUSTER_WORK}/snp_repair_rerun/r308
VT=vcftools; G=$W/groups_K$K$T; N=$W/fig4c_K$K$T
mkdir -p $N/out/fst $N/out/pima $N/logs; J=$N/jobs.txt; : > $J
for f in $W/fig4c/vcf/chr*.present.vcf.gz; do c=$(basename $f .present.vcf.gz)
  echo "python3 $R/pi_ma_gen.py $f $G $N/out/pima/$c > $N/logs/pima_$c.log 2>&1" >> $J
  for a in $(seq 1 $K); do for b in $(seq $((a+1)) $K); do
    echo "$VT --gzvcf $f --weir-fst-pop $G/K${K}_Pop$a.txt --weir-fst-pop $G/K${K}_Pop$b.txt --fst-window-size 100000 --fst-window-step 50000 --out $N/out/fst/${c}__K${K}_Pop${a}_K${K}_Pop${b}_100kb_fst > $N/logs/f_${c}_$a$b.log 2>&1" >> $J
  done; done
done
nice -n 10 xargs -P $P -I{} bash -c "{}" < $J
python $R/summ_K.py $K $T > $N/logs/summ.log 2>&1
echo DONE $(date -Is) > $N/done
