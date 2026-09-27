#!/bin/bash
D=${ANALYSIS_DIR}/08_hifi_chromosome/11_final_corrected_39genomes_annotations_20260812/attempt_20260812_01/01_genomes
W=${CLUSTER_WORK}/ole16_cis; cd $W; mkdir -p q loc tmp
ST=samtools
$ST faidx $D/African_hap2.fa chr11:112170000-112215000 | sed 's/^>.*/>q45/' > q/query45k.fa
$ST faidx $D/African_hap2.fa chr11:112186961-112187344 | sed 's/^>.*/>OLE16a_FLHap2_minus/' > q/ole16a_gene_fwdstrand.fa
$ST faidx $D/African_hap2.fa chr11:112199288-112210499 | sed 's/^>.*/>G1498/' > q/g1498.fa
$ST faidx $D/African_hap2.fa chr11:112175522-112184043 | sed 's/^>.*/>G1496/' > q/g1496.fa
ls $D/*.fa | xargs -n1 basename | sed 's/\.fa$//' > genomes.txt; wc -l genomes.txt
cat genomes.txt | xargs -P 8 -I{} bash c03_extract.sh {}
cat loc/*.status > loc_status.tsv; echo done
