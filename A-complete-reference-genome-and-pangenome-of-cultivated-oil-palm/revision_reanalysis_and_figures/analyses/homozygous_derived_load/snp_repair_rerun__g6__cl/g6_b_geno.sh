#!/bin/bash
# extract 308-sample biallelic SNP genotypes (QUAL>=30) in CDS, per chromosome
set -euo pipefail
c=$1
O=${CLUSTER_WORK}/snp_repair_rerun/g6
E=${CLUSTER_HOME}/miniconda3/envs/cs213/bin
IN=${CLUSTER_WORK}/snp_repair_rerun/vcf/$c.population_snps.corrected.vcf.gz
$E/bcftools view --threads 2 -S ${CLUSTER_WORK}/fig4c_final/meta/s308.txt -R $O/bed/$c.bed -i 'QUAL>=30' -m2 -M2 -v snps -Ou $IN \
 | $E/bcftools query -H -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' | gzip > $O/geno/$c.tsv.gz
echo "$c done $(date)" >> $O/geno/done.txt
