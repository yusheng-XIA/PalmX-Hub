#!/bin/bash
# LD decay: 2-kb thinned biallelic SNPs of the Fig. 4c SNP set (308, QUAL>=30, MAF>=0.05, missing<=0.2), PopLDdecay as published
set -euo pipefail
W=${CLUSTER_WORK}/snp_repair_r308; E=${CLUSTER_HOME}/miniconda3/envs/genomics_a2/bin
IMG=${ANALYSIS_DIR}/05_GWAS/GWAS_tur/software/software/Reseq_genek.sif
cd $W/ld
ls $W/fig4c/vcf/chr*.present.vcf.gz | sort -V > in.list
$E/bcftools concat -f in.list -Ov | awk -v step=2000 "BEGIN{FS=OFS=\"\t\"} /^#/{print;next} {t++; if(length(\$4)!=1||length(\$5)!=1||\$5~/,/)next; if(!(\$1 in last)||\$2-last[\$1]>=step){print; last[\$1]=\$2; s++}} END{print \"total=\"t\" selected=\"s > \"/dev/stderr\"}" > thin2kb.vcf 2> thin.log
for g in K4_Pop1 K4_Pop2 K4_Pop3 K4_Pop4; do
 singularity exec -B ${DATA_ROOT} $IMG /opt/PopLDdecay/PopLDdecay -InVCF thin2kb.vcf -SubPop lists/$g.list -MaxDist 500 -MAF 0.005 -Het 0.9 -Miss 0.25 -OutStat stats/$g.stat > logs/$g.log 2>&1 &
done; wait
python3 plot_k4_ld_decay_thin2kb_smooth.py > logs/plot.log 2>&1 || python plot_k4_ld_decay_thin2kb_smooth.py > logs/plot.log 2>&1
echo DONE > s5.done
