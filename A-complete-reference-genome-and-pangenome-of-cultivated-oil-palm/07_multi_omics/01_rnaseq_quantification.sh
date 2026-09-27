#!/bin/bash
# Bulk RNA-seq: 114 libraries of the replicated FL/TN series (19 stages x 3 biological replicates), plus the
# single-library TK and NS series used only for exploratory comparisons. fastp + STAR; DESeq2 downstream.
set -euo pipefail
threads=16
ref=FL_Hap2

STAR --runMode genomeGenerate --runThreadN ${threads} --genomeDir star_${ref} \
    --genomeFastaFiles ${ref}.fa --sjdbGTFfile ${ref}.gtf --sjdbOverhang 149
while read lib; do
    fastp -i ${lib}_R1.fq.gz -I ${lib}_R2.fq.gz -o clean/${lib}_R1.fq.gz -O clean/${lib}_R2.fq.gz \
        -w ${threads} -j qc/${lib}.fastp.json -h qc/${lib}.fastp.html
    STAR --runThreadN ${threads} --genomeDir star_${ref} --readFilesCommand zcat \
        --readFilesIn clean/${lib}_R1.fq.gz clean/${lib}_R2.fq.gz \
        --outSAMtype BAM SortedByCoordinate --quantMode GeneCounts --outFileNamePrefix star/${lib}.
done < rnaseq_libraries.txt
# Gene-level counts -> DESeq2 normalisation, variance-stabilising transformation or differential expression
