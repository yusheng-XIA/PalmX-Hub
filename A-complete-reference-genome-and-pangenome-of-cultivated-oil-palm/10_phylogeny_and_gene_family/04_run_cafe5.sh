#!/bin/bash
# Gene-family expansion and contraction (CAFE5 v5.1.0) on the 12 monocot genomes
# gene_families.tsv: OrthoFinder Orthogroups.GeneCount.tsv reformatted to the CAFE5 input
#                    (Desc, Family ID, one column per species of the pruned tree)
set -euo pipefail
python 03_make_cafe_tree.py FigTree.tre > cafe_tree.nwk

# single birth-death rate with gamma-distributed rate variation among families (3 categories)
cafe5 --infile gene_families.tsv --tree cafe_tree.nwk -k 3 --cores 8 --output_prefix cafe_gamma
# Families with conditional P < 0.05 (Gamma_family_results.txt) are significantly expanded or contracted;
# node-level gains and losses are summarised from Gamma_change.tab / Gamma_asr.tre.
