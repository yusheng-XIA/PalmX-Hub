#!/bin/bash
# pi/FST (Fig. 4c, ED6b) and LD decay (ED6a) on the 308-based Fig. 4c SNP set with the new K = 4 groups
set -euo pipefail
W=${CLUSTER_WORK}/snp_repair_r308; E=${CLUSTER_HOME}/miniconda3/envs/genomics_a2/bin
export TMPDIR=$W/tmp
mkdir -p $W/ld/logs; cd $W/ld
ls $W/fig4c/vcf/chr*.present.vcf.gz | sort -V > in.list
$E/bcftools concat -f in.list -Ov | awk -v step=2000 "BEGIN{FS=OFS=\"\t\"} /^#/{print;next} {t++; if(length(\$4)!=1||length(\$5)!=1||\$5~/,/)next; if(!(\$1 in last)||\$2-last[\$1]>=step){print; last[\$1]=\$2; s++}} END{print \"total=\"t\" selected=\"s > \"/dev/stderr\"}" > thin2kb.vcf 2> thin.log
rm -rf $W/fig4c_k4new/out $W/fig4c_k4new/vcf $W/ld_k4new/stats
bash $W/cl/run_k4new.sh > $W/logs/run_k4new.log 2>&1
echo DONE $(date -Is) > $W/logs/pifst.done
