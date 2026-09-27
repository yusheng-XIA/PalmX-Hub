#!/bin/bash
c=$1
O=${CLUSTER_WORK}/headline_pop/g6
SE=${JOINT_SNP_DIR}/09_variant_annotation
E=${CLUSTER_HOME}/miniconda3/envs/cs213/bin
$E/bcftools view -G -T $O/bed/$c.bed -v snps -Ou $SE/$c.diversity_snps.snpeff.vcf.gz | $E/bcftools query -f '%POS\t%REF\t%ALT\t%ANN\n' | awk -F'\t' 'BEGIN{OFS="\t"}{split($4,a,","); split(a[1],b,"|"); print $1,$2,$3,b[2],b[4]}' > $O/snpeff/$c.tsv
echo "$c ok" >> $O/snpeff/done.txt
