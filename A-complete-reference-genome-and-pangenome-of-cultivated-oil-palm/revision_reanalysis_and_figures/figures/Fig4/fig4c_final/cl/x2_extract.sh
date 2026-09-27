#!/bin/bash
# usage: x2_extract.sh CHROM
set -euo pipefail
c=$1
W=${CLUSTER_WORK}/fig4c_final
E=${CLUSTER_HOME}/miniconda3/envs/wes_env/bin
FJOINT=${JOINT_SNP_DIR}/06_filter
IN=$FJOINT/$c.population_snps.corrected.vcf.gz
O=$W/vcf; mkdir -p $O $W/logs
[ -s $O/$c.strict.vcf.gz.tbi ] && [ -s $O/$c.present.vcf.gz.tbi ] && exit 0
t0=$(date +%s)
$E/bcftools view --threads 2 -S $W/meta/s308.txt -i 'QUAL>=30' -Ou $IN \
 | $E/bcftools annotate -x INFO,^FORMAT/GT -Ou \
 | $E/bcftools +fill-tags -Ou -- -t AN,AC,MAF \
 | $E/bcftools view -i 'INFO/MAF>=0.05' -Ov \
 | python3 $W/cl/split308.py $c $W/meta/absent_308.tsv $O/$c.strict.vcf.gz.part $O/$c.present.vcf.gz.part $W/logs/split_$c.txt
mv $O/$c.strict.vcf.gz.part $O/$c.strict.vcf.gz; mv $O/$c.present.vcf.gz.part $O/$c.present.vcf.gz
$E/tabix -f -p vcf $O/$c.strict.vcf.gz; $E/tabix -f -p vcf $O/$c.present.vcf.gz
echo -e "$c\t$(( $(date +%s)-t0 ))s\t$(date '+%F %T')\t$(hostname)" >> $W/logs/extract2_times.tsv
