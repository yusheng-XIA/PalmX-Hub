#!/bin/bash
# Homology-based gene prediction
threads=64

# ============================================================
#  miniprot (fast protein-to-genome alignment)
# ============================================================
miniprot -t ${threads} \
    --gff \
    ${sample}.nucleus.masked.fasta \
    ${palm_proteomes} \
    > ${sample}.miniprot.gff3
