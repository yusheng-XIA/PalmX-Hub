#!/bin/bash
# First SnpEff annotation of the population call set at CDS sites of one chromosome (sensitivity analysis;
# chr01B, chr08B, chr10B, chr14B and chr16B). Output columns: pos ref alt effect impact gene
# usage: 03_snpeff_first_annotation.sh chr01B
set -euo pipefail
c=$1
O=${WORK:-coding_zygosity}
SE=${SNPEFF_DIR:-snpeff}   # SnpEff-annotated population VCFs (<chrom>.diversity_snps.snpeff.vcf.gz)
mkdir -p $O/snpeff
bcftools view -G -T $O/bed/$c.bed -v snps -Ou $SE/$c.diversity_snps.snpeff.vcf.gz \
    | bcftools query -f '%POS\t%REF\t%ALT\t%ANN\n' \
    | awk -F'\t' 'BEGIN{OFS="\t"}{split($4,a,","); split(a[1],b,"|"); print $1,$2,$3,b[2],b[3],b[4]}' \
    | gzip > $O/snpeff/$c.ann.tsv.gz
