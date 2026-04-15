#!/bin/bash
# Tandem Repeats (TRF) and Segmental Duplications (BISER)

# ============================================================
# Tandem Repeats Finder
# ============================================================
trf ${sample}.fasta 2 6 6 80 10 50 2000 -h
python TRF2GFF.py -d ${sample}.fasta.2.6.6.80.10.50.2000.dat -o ${sample}.tr.gff
grep -v '^#' ${sample}.tr.gff | sed '/^$/d' \
    | awk '{OFS="\t"}{print $1,$4,$5,$9}' > ${sample}.tr.bed

# ============================================================
# Segmental Duplications (BISER)
# ============================================================
biser -o ${sample}_sd --threads ${threads} ${sample}.fasta
awk '{OFS="\t"}{print $1,$4,$5,$9}' ${sample}_sd.gff3 > ${sample}.sd.bed
