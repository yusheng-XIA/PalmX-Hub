#!/bin/bash
D=${ANALYSIS_DIR}/08_hifi_chromosome/11_final_corrected_39genomes_annotations_20260812/attempt_20260812_01
W=${CLUSTER_WORK}/ole16_cis; mkdir -p $W/tmp $W/out; cd $W
for g in African_hap2 American_hap1 BK_hap1 BK_hap2; do echo "== $g"; grep -P "\tgene\t" $D/02_annotations/$g.gene.gff3 | grep -E "chr11B\.1497|chr11A\.1059|chr9\.1211|chr9\.414|chr04B\.697|chr04A" | head -5; done
head -3 $D/02_annotations/African_hap2.gene.gff3
cut -f1 $D/01_genomes/African_hap2.fa.fai | head -20 | tr "\n" " "; echo
cut -f1 $D/01_genomes/MZ4_hap1.fa.fai | head -20 | tr "\n" " "; echo
cut -f1 $D/01_genomes/EG_008.fa.fai | head -20 | tr "\n" " "; echo
ls $D/00_manifest; head -5 $D/03_coordinate_maps/* 2>/dev/null | head -20
