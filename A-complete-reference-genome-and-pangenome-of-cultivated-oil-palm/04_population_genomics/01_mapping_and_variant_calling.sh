#!/bin/bash
# Population genomics: read mapping and variant calling
# 308 accessions, Illumina PE150, ~60x coverage
threads=64

# ============================================================
# 1. Read mapping (BWA-MEM2)
# ============================================================
bwa-mem2 index ${reference}.fasta

for sample in $(cat accession_list.txt); do
    bwa-mem2 mem -t ${threads} ${reference}.fasta \
        ${sample}_R1.fq.gz ${sample}_R2.fq.gz \
        | samtools sort -@ ${threads} -o ${sample}.sorted.bam

    samtools markdup ${sample}.sorted.bam ${sample}.markdup.bam
    samtools index ${sample}.markdup.bam
done

# ============================================================
# 2. Variant calling (GATK HaplotypeCaller)
# ============================================================
for sample in $(cat accession_list.txt); do
    gatk HaplotypeCaller \
        -R ${reference}.fasta \
        -I ${sample}.markdup.bam \
        -O ${sample}.g.vcf.gz \
        -ERC GVCF \
        --native-pair-hmm-threads ${threads}
done

# ============================================================
# 3. Joint genotyping
# ============================================================
# Create sample map
for sample in $(cat accession_list.txt); do
    echo -e "${sample}\t${sample}.g.vcf.gz"
done > sample_map.txt

gatk GenomicsDBImport \
    --sample-name-map sample_map.txt \
    --genomicsdb-workspace-path genomics_db \
    -L intervals.list

gatk GenotypeGVCFs \
    -R ${reference}.fasta \
    -V gendb://genomics_db \
    -O raw_variants.vcf.gz

# ============================================================
# 4. Hard filtering
# ============================================================
gatk VariantFiltration \
    -R ${reference}.fasta \
    -V raw_variants.vcf.gz \
    --filter-expression "QD < 2.0 || MQ < 40.0 || FS > 60.0 || SOR > 3.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" \
    --filter-name "HARD_FILTER" \
    -O filtered_variants.vcf.gz

# Retain biallelic SNPs, MAF > 0.01, missing < 20%
gatk SelectVariants \
    -R ${reference}.fasta \
    -V filtered_variants.vcf.gz \
    --select-type-to-include SNP \
    --restrict-alleles-to BIALLELIC \
    -O biallelic_snps.vcf.gz

vcftools --gzvcf biallelic_snps.vcf.gz \
    --maf 0.01 --max-missing 0.8 \
    --recode --recode-INFO-all \
    --out final_snps
