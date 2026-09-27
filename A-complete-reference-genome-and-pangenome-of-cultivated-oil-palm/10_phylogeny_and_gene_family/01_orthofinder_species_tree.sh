#!/bin/bash
# Orthogroups and species tree for the 30 genomes of the phylogenomic analysis (Fig. 2)
# OrthoFinder v3.1.1 (MSA mode), trimAl v1.4, IQ-TREE v1.6.12
set -euo pipefail
threads=48

# 1. Orthogroups (one longest protein per gene in proteomes/)
orthofinder -f proteomes/ -M msa -t ${threads} -a ${threads} -o orthofinder_out

# 2. OrthoFinder concatenated species-tree alignment (155 orthogroups single-copy in >= 8 of 30 genomes;
#    51,285 columns); genomes without a single-copy gene in an orthogroup are coded as missing data
aln=orthofinder_out/Results_*/MultipleSequenceAlignments/SpeciesTreeAlignment.fa
trimal -in ${aln} -out SpeciesTreeAlignment.trim.fa -gt 0.6 -cons 60      # 30,771 amino-acid sites

# 3. Maximum-likelihood species tree (JTT+F+R3 selected by ModelFinder, BIC)
iqtree -s SpeciesTreeAlignment.trim.fa -m MFP -bb 1000 -nt AUTO -pre species_tree
#   final model: iqtree -s SpeciesTreeAlignment.trim.fa -m JTT+F+R3 -bb 1000 -pre species_tree

# 4. MCMCTree input: PHYLIP alignment (supergene.phy) and the calibrated topology (calibrated_tree.txt).
#    The eudicot topology of the IQ-TREE tree was corrected to established relationships (Brassica napus
#    within Brassicaceae) before adding the seven calibrations listed in calibration_points.tsv.
