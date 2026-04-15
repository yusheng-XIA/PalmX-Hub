#!/bin/bash
# Ab initio gene prediction using multiple tools
threads=64

# ============================================================
# 1. Repeat masking
# ============================================================
RepeatMasker -nolow -species rice -pa ${threads} ${sample}.nucleus.fasta
# Output: ${sample}.nucleus.masked.fasta

# ============================================================
# 2. BRAKER
# ============================================================
braker.pl --genome=${sample}.nucleus.masked.fasta \
    --bam=${sample}.rnaseq.bam \
    --species=${sample} \
    --cores=${threads} \
    --gff3

# ============================================================
# 3. FGENESH
# ============================================================
fgenesh ${sample}.nucleus.masked.fasta -org monocot \
    > ${sample}.fgenesh.out
fgenesh_to_GFF3.pl ${sample}.fgenesh.out > ${sample}.fgenesh.gff3

# ============================================================
# 4. Helixer
# ============================================================
Helixer.py --model-filepath land_plant.h5 \
    --fasta-path ${sample}.nucleus.masked.fasta \
    --gff-output-path ${sample}.helixer.gff3

# ============================================================
# 5. EviAnn
# ============================================================
eviann -g ${sample}.nucleus.masked.fasta \
    -b ${sample}.rnaseq.bam \
    -t ${threads} \
    -o ${sample}.eviann.gff3
