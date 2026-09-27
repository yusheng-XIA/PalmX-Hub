#!/bin/bash
# Homology evidence: Swiss-Prot plant subset and published palm proteomes aligned with miniprot v0.12
set -euo pipefail
threads=64
cat uniprot_sprot_plants.fa palm_proteomes/*.fa > homology_proteins.fa
miniprot -t ${threads} --gff -I ${sample}.softmasked.fa homology_proteins.fa > ${sample}.miniprot.gff3
