#!/bin/bash
# Multi-omics profiling: RNA-seq quantification and ASE analysis
# 4 varieties × 19 stages = 76 samples
threads=64

# ============================================================
# 1. RNA-seq alignment and quantification
# ============================================================
hisat2-build -p ${threads} ${reference}.fasta ${reference}_index

for sample in $(cat sample_list_76.txt); do
    # Quality control
    fastp -i ${sample}_R1.fq.gz -I ${sample}_R2.fq.gz \
        -o ${sample}_clean_R1.fq.gz -O ${sample}_clean_R2.fq.gz \
        -w ${threads} --json ${sample}.fastp.json

    # Alignment
    hisat2 -p ${threads} --dta --max-intronlen 20000 \
        -x ${reference}_index \
        -1 ${sample}_clean_R1.fq.gz \
        -2 ${sample}_clean_R2.fq.gz \
        | samtools sort -@ ${threads} -o ${sample}.bam

    # Quantification (TPM)
    stringtie -e -B -p ${threads} \
        -G ${reference}.gtf \
        -o ${sample}.gtf \
        ${sample}.bam
done

# ============================================================
# 2. Allele-specific expression (ASE) with STAR
# ============================================================
STAR --runMode genomeGenerate \
    --runThreadN ${threads} \
    --genomeDir star_allele_index \
    --genomeFastaFiles ${reference}.fasta \
    --sjdbGTFfile ${reference}.allele_aware.gtf

for sample in $(cat sample_list_76.txt); do
    STAR --runThreadN ${threads} \
        --genomeDir star_allele_index \
        --readFilesIn ${sample}_clean_R1.fq.gz ${sample}_clean_R2.fq.gz \
        --readFilesCommand zcat \
        --alignIntronMax 20000 \
        --alignMatesGapMax 25000 \
        --outFilterMultimapNmax 1 \
        --outSAMtype BAM SortedByCoordinate \
        --outFileNamePrefix ${sample}_ase_

    stringtie -e -B -p ${threads} \
        -G ${reference}.allele_aware.gtf \
        -o ${sample}_ase.gtf \
        ${sample}_ase_Aligned.sortedByCoord.out.bam
done
