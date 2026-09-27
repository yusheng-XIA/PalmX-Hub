#!/bin/bash
# Genotyping the 308 resequenced accessions against the graph-pangenome variant set (PanGenie v4.2.1)
set -euo pipefail
threads=24

PanGenie-index -v oilpalm39.pangenie_ready.vcf -r FL_Hap2.fa -t ${threads} -o pgindex     # default k = 31
while read sample; do
    zcat ${sample}_R1.fq.gz ${sample}_R2.fq.gz > ${sample}.merged.fq
    PanGenie -f pgindex -i ${sample}.merged.fq -s ${sample} -j ${threads} -t ${threads} -o genotypes/${sample}
    rm ${sample}.merged.fq
done < samples_308.txt

bcftools merge -m none genotypes/*_genotyping.vcf.gz -Oz -o pangenie_308.vcf.gz
# SV records: REF or ALT allele >= 50 bp, MAF >= 0.05, missingness <= 0.1 (370,136 SVs for population
# analyses and SV-GWAS)
bcftools view -i '(strlen(REF)>=50 || strlen(ALT)>=50)' pangenie_308.vcf.gz -Ou \
    | bcftools view -q 0.05:minor -i 'F_MISSING<=0.1' -Oz -o pangenie_sv_308.maf05_miss10.vcf.gz
