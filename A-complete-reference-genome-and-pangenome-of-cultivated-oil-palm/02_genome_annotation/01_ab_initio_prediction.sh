#!/bin/bash
# Initial gene-model prediction: BRAKER3, EviAnn, GeneMark.hmm (v4.68) and Helixer
# The genome is soft-masked with the EDTA TE library of the assembly (03_repeat_annotation).
set -euo pipefail
threads=64

RepeatMasker -pa ${threads} -s -xsmall -gff -lib ${sample}.EDTA.TElib.fa ${sample}.fa
mv ${sample}.fa.masked ${sample}.softmasked.fa

# BRAKER3 with RNA-seq alignments (03_transcriptome_prediction.sh) and protein hints
braker.pl --genome=${sample}.softmasked.fa --bam=${sample}.rnaseq.merged.bam \
    --prot_seq=plant_proteins.fa --species=${sample} --threads=${threads} --gff3
# braker/Augustus.hints.gff3 -> AUGUSTUS evidence for EVM

# EviAnn (evidence-based annotation from RNA-seq and proteins)
eviann.sh -t ${threads} -g ${sample}.softmasked.fa -r rnaseq_reads.txt -p plant_proteins.fa

# GeneMark.hmm v4.68 (GeneMark-ES self-training)
gmes_petap.pl --ES --sequence ${sample}.softmasked.fa --cores ${threads}

# Helixer (land_plant model)
Helixer.py --lineage land_plant --fasta-path ${sample}.softmasked.fa --gff-output-path ${sample}.helixer.gff3
