#!/bin/bash
# extract 308-sample biallelic SNP genotypes (QUAL>=30) in CDS, per chromosome
set -euo pipefail
c=$1
O=${WORK:-coding_zygosity}
IN=${SNP_VCF:-snp_778.population_snps.vcf.gz}   # 778-accession call set after filtering (04_population_genomics/01)
bcftools view --threads 2 -S samples_308.txt -R $O/bed/$c.bed -i 'QUAL>=30' -m2 -M2 -v snps -Ou $IN \
 | bcftools query -H -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' | gzip > $O/geno/$c.tsv.gz
echo "$c done $(date)" >> $O/geno/done.txt
