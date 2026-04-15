#!/bin/bash
# TE annotation using EDTA + RepeatModeler + RepeatMasker
threads=64

# ============================================================
# 1. EDTA pipeline
# ============================================================
EDTA.pl --genome ${sample}.fasta \
    --species others \
    --sensitive 0 \
    --anno 1 \
    --threads ${threads}

# Convert GFF to BED
grep -v "Parent" ${sample}.fasta.mod.EDTA.TEanno.gff3 \
    | awk '{OFS="\t"}{print $1,$4,$5,$3}' > ${sample}.te.bed

# ============================================================
# 2. LTR insertion time estimation
# ============================================================
# Mutation rate: 6.5e-9 substitutions/site/year
# T = D / (2 * mu)
LTR_retriever -genome ${sample}.fasta \
    -inharvest ${sample}_rawLTR.scn \
    -threads ${threads} \
    -maxlenltr 15000
