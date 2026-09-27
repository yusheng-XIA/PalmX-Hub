#!/bin/bash
D=${ANALYSIS_DIR}/08_hifi_chromosome/11_final_corrected_39genomes_annotations_20260812/attempt_20260812_01
echo "== FL-Hap2 chr11 112.10-112.28"; awk -F'\t' '$1=="chr11" && $3=="gene" && $4>112100000 && $5<112280000' $D/02_annotations/African_hap2.gene.gff3 | cut -f4,5,7,9 | cut -c1-80
echo "== FL-Hap1 chr11 111.12-111.30"; awk -F'\t' '$1=="chr11" && $3=="gene" && $4>111120000 && $5<111300000' $D/02_annotations/American_hap1.gene.gff3 | cut -f4,5,7,9 | cut -c1-80
echo "== FL-Hap2 OLE16a mRNA/CDS"; grep "chr11B.1497" $D/02_annotations/African_hap2.gene.gff3 | cut -f3,4,5,7 
echo "== FL-Hap1 OLE16a"; grep "chr11A.1059" $D/02_annotations/American_hap1.gene.gff3 | cut -f3,4,5,7
