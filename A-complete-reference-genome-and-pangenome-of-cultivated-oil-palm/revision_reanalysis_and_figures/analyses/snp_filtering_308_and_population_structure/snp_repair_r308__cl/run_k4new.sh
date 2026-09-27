#!/bin/bash
# Recompute π/FST (Fig.4c/ED6b) and LD decay (ED6a) with the repaired-data K=4 groups.
W=${CLUSTER_WORK}/snp_repair_r308
G=$W/groups_k4new
N=$W/fig4c_k4new; mkdir -p $N/{vcf,logs/vt,out/present,out/pima_present,meta}
ln -sfn $G $N/groups; cp $W/fig4c/cl -r $N/; cp $W/fig4c/meta/* $N/meta/; cp $W/fig4c/logs/split_*.txt $N/logs/
for f in $W/fig4c/vcf/chr*.present.vcf.gz*; do ln -sf $(readlink -f $f) $N/vcf/; done
ln -sfn present $N/out/strict; ln -sfn pima_present $N/out/pima_strict
VT=vcftools
J=$N/logs/jobs.txt; : > $J
for f in $N/vcf/chr*.present.vcf.gz; do c=$(basename $f .present.vcf.gz)
  for p in 1 2 3 4; do echo "$VT --gzvcf $f --keep $G/K4_Pop$p.txt --window-pi 100000 --window-pi-step 50000 --out $N/out/present/${c}__K4_Pop${p}_100kb.pi > $N/logs/vt/${c}_P$p.log 2>&1" >> $J; done
  for pr in "1 2" "1 3" "1 4" "2 3" "2 4" "3 4"; do set -- $pr; echo "$VT --gzvcf $f --weir-fst-pop $G/K4_Pop$1.txt --weir-fst-pop $G/K4_Pop$2.txt --fst-window-size 100000 --fst-window-step 50000 --out $N/out/present/${c}__K4_Pop$1_K4_Pop$2_100kb_fst > $N/logs/vt/${c}_F$1$2.log 2>&1" >> $J; done
  echo "python3 $N/cl/pi_ma.py $f $G $N/out/pima_present/$c > $N/logs/pima_$c.log 2>&1" >> $J
done
# LD with new lists
IMG=${ANALYSIS_DIR}/05_GWAS/GWAS_tur/software/software/Reseq_genek.sif
L=$W/ld_k4new; mkdir -p $L/stats $L/logs $L/lists
for p in 1 2 3 4; do cp $G/K4_Pop$p.txt $L/lists/K4_Pop$p.list; echo "singularity exec -B ${DATA_ROOT} $IMG /opt/PopLDdecay/PopLDdecay -InVCF $W/ld/thin2kb.vcf -SubPop $L/lists/K4_Pop$p.list -MaxDist 500 -MAF 0.005 -Het 0.9 -Miss 0.25 -OutStat $L/stats/K4_Pop$p.stat > $L/logs/K4_Pop$p.log 2>&1" >> $J; done
nice -n 10 xargs -P 24 -I{} bash -c "{}" < $J
sed 's#snp_repair_r308/ld")#snp_repair_r308/ld_k4new")#' $W/ld/plot_k4_ld_decay_thin2kb_smooth.py > $L/plot.py
sed -i 's/"sample_n"/"sample_n"/' $L/plot.py
python $L/plot.py > $L/logs/plot.log 2>&1
cd $N && python3 cl/s2_final.py > logs/s2_final.log 2>&1
echo DONE > $W/k4new.done
