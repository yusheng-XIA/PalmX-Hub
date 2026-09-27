#!/bin/bash
# usage: enhB_00_mm2.sh H1|H2   ; minimap2 asm20 alignment of E. oleifera haplotype to Africa_hap2 (same preset as Phoenix)
set -euo pipefail
W=work/outgroup_sensitivity
MM=minimap2
REF=dsv_analysis/input/FL_Hap2.fa
Q=assemblies/Eoleifera_$1.chromosomes.fa   # independent E. oleifera accession, H1 or H2
mkdir -p $W/aln
$MM -x asm20 -t 16 -c --cs --secondary=no $REF $Q > $W/aln/MZ4$1_vs_Africa_hap2.asm20.cs.primary.paf 2> $W/aln/MZ4$1.mm2.log
echo DONE > $W/aln/MZ4$1.done
