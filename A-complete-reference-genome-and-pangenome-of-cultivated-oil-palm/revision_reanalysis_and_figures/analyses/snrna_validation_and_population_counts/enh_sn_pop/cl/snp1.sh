#!/bin/bash
# usage: snp1.sh CHROM   (corrected joint-called population SNP set; 308 samples; QUAL>=30; biallelic SNPs)
set -euo pipefail
c=$1
W=${CLUSTER_WORK}/enh_sn_pop
E=${CLUSTER_HOME}/miniconda3/envs/cs213/bin
FJOINT=${JOINT_SNP_DIR}/06_filter
IN=$FJOINT/oil_palm_joint.population_snps.corrected.vcf.gz
G="ALL COMM NONC AFR HHG IDB SAEG SEAA SEAB K4P1 K4P2 K4P3 K4P4"
Q='%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n'
t0=$(date +%s)
$E/bcftools view --threads 2 -S ${CLUSTER_WORK}/fig4c_final/meta/s308.txt -i 'QUAL>=30' -m2 -M2 -v snps -r $c -Ou $IN \
 | $E/bcftools annotate -x INFO,^FORMAT/GT -Ou \
 | $E/bcftools query -H -f "$Q" \
 | python3 $W/cl/agg.py snp $c $W/cl/groups308.txt ${CLUSTER_WORK}/fig4c_final/meta/absent_308.tsv 0.8 $W/snp/$c
echo -e "$c\t$(( $(date +%s)-t0 ))s\t$(hostname)\t$(date '+%F %T')" >> $W/logs/snp_times.tsv
