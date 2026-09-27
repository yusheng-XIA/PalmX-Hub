#!/bin/bash
# Single-nucleus RNA-seq (10x Genomics Single Cell 3' v4): FL and TN mesocarp at 95, 125 and 185 d
# Reads are aligned to FL-Hap2 with its gene annotation; filtered feature-barcode matrices are used downstream.
set -euo pipefail
threads=32

# Reference (GTF converted from the FL-Hap2 GFF3)
gffread FL_Hap2.gff3 -T -o FL_Hap2.gtf
cellranger mkref --genome=FL_Hap2 --fasta=FL_Hap2.fa --genes=FL_Hap2.gtf --nthreads=${threads}

for sample in FL_95d FL_125d FL_185d TN_95d TN_125d TN_185d; do
    cellranger count \
        --id=${sample} \
        --transcriptome=FL_Hap2 \
        --fastqs=fastq/${sample} \
        --sample=${sample} \
        --localcores=${threads} \
        --localmem=128
done
# Called nuclei (all retained downstream): FL_95d 15,579; FL_125d 14,192; FL_185d 6,836;
# TN_95d 17,249; TN_125d 11,836; TN_185d 12,501 (78,193 in total).
