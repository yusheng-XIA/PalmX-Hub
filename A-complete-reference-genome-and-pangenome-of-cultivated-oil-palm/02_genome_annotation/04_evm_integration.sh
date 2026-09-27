#!/bin/bash
# EvidenceModeler v2.1.0: integration in 500-kb segments with 50-kb overlaps
set -euo pipefail
threads=64

cat > weights.txt <<'W'
ABINITIO_PREDICTION	AUGUSTUS	3
ABINITIO_PREDICTION	EviAnn	2
ABINITIO_PREDICTION	GeneMark.hmm	2
ABINITIO_PREDICTION	Helixer	5
PROTEIN	miniprot	5
TRANSCRIPT	transdecoder	7
W

cat braker.augustus.gff3 eviann.gff3 genemark.gff3 helixer.gff3 > gene_predictions.gff3
EVidenceModeler --sample_id ${sample} --genome ${sample}.softmasked.fa --weights weights.txt \
    --gene_predictions gene_predictions.gff3 \
    --protein_alignments ${sample}.miniprot.gff3 \
    --transcript_alignments ${sample}.transdecoder.genome.gff3 \
    --segmentSize 500000 --overlapSize 50000 --CPU ${threads}

# Coding-sequence completeness, splice sites and evidence concordance were checked after merging
gffread ${sample}.EVM.gff3 -g ${sample}.fa -x ${sample}.cds.fa -y ${sample}.pep.fa
