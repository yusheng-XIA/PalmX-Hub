#!/bin/bash
# Extract the 308 Fig.4c samples from the final corrected callset (joint-called parent set),
# refilter within 308: QUAL>=30, F_MISSING<=0.2 (== vcftools --max-missing 0.8), MAF>=0.05; GT only; '.'->'./.'
# usage: x1_extract.sh CHROM SET(population|diversity) [REGION]
set -euo pipefail
c=$1; set=$2; reg=${3:-}
W=${CLUSTER_WORK}/fig4c_final
E=${CLUSTER_HOME}/miniconda3/envs/wes_env/bin
FJOINT=${JOINT_SNP_DIR}/06_filter
IN=$FJOINT/$c.${set}_snps.corrected.vcf.gz
tag=$c; [ -n "$reg" ] && tag=${c}_test_${set}
O=$W/vcf; mkdir -p $O $W/logs
R=""; [ -n "$reg" ] && R="-r $reg"
t0=$(date +%s)
$E/bcftools view --threads 2 $R -S $W/meta/s308.txt -i 'QUAL>=30' -Ou $IN \
 | $E/bcftools annotate -x INFO,^FORMAT/GT -Ou \
 | $E/bcftools +fill-tags -Ou -- -t AN,AC,MAF,F_MISSING \
 | $E/bcftools view -i 'INFO/F_MISSING<=0.2 && INFO/MAF>=0.05' -Ov \
 | python3 $W/cl/fixgt.py $W/logs/fixgt_$tag.txt \
 | $E/bgzip -@ 2 -c > $O/$tag.308.vcf.gz.part
mv $O/$tag.308.vcf.gz.part $O/$tag.308.vcf.gz
$E/tabix -f -p vcf $O/$tag.308.vcf.gz
echo -e "$tag\t$set\t$(( $(date +%s)-t0 ))s\t$(date '+%F %T')" >> $W/logs/extract_times.tsv
