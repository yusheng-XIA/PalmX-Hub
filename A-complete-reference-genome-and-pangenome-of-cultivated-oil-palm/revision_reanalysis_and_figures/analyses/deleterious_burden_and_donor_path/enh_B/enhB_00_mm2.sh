#!/bin/bash
# usage: enhB_00_mm2.sh H1|H2   ; minimap2 asm20 alignment of E. oleifera haplotype to Africa_hap2 (same preset as Phoenix)
set -euo pipefail
W=${CLUSTER_WORK}/enh_B
MM=${DATA_DIR2}/tools/TGS-GapCloser/minimap2
REF=${ANALYSIS_DIR}/21_MS/06_result/dSVs/input/Africa_hap2.fa
Q=${DATA_DIR2}/projects/1-oil_palm/08-Nature_review/7-Fig1/06_haplotype_species_pav/inputs/standardized/MZ4$1.chromosomes.fa
mkdir -p $W/aln
$MM -x asm20 -t 16 -c --cs --secondary=no $REF $Q > $W/aln/MZ4$1_vs_Africa_hap2.asm20.cs.primary.paf 2> $W/aln/MZ4$1.mm2.log
echo DONE > $W/aln/MZ4$1.done
