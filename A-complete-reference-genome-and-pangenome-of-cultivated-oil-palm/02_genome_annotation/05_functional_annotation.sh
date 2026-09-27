#!/bin/bash
# Functional annotation and annotation support
set -euo pipefail
threads=64

# GO / KEGG with eggNOG-mapper v2.1.13
emapper.py -i ${sample}.pep.fa -o ${sample} --cpu ${threads} -m diamond

# Annotation completeness: BUSCO v6.0.0, embryophyta_odb12
busco -m protein -i ${sample}.pep.fa -l embryophyta_odb12 -o ${sample}_protein_busco -c ${threads} --offline

# Pfam support: HMMER v3.4 hmmsearch with model-specific gathering thresholds
hmmsearch --cut_ga --noali --cpu ${threads} --domtblout ${sample}.pfam.domtbl Pfam-A.hmm ${sample}.pep.fa > /dev/null
grep -v '^#' ${sample}.pfam.domtbl | awk '{print $1}' | sort -u > ${sample}.pfam_supported.txt

# FL-Hap1 vs FL-Hap2: a high-confidence gene has >= 2 of 4 evidence classes
#   Pfam support | robust RNA-seq support | protein-homology support | strong reciprocal-best-hit 1:1 orthology
# Gene-model characteristics (protein length, exon number, CDS overlap with TE annotation) are summarised
# with the same gene lists (Supplementary Data 3).
