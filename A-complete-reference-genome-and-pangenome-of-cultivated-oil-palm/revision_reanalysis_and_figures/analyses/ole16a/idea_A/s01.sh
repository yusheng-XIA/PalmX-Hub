#!/bin/bash
W=${CLUSTER_WORK}/idea_A; mkdir -p $W/out; cd $W
B=${ANALYSIS_DIR}
P=$B/22_answer_reviews/00_ms/03_V3/03_figure3/07_proteomics_multiomics_20260724
R=$P/01_reference_database
head -3 $R/FL_TN_unified_exact_sequence_nr_map.tsv
grep -P "UFTN00833[0-4]\t|UFTN008329\t|UFTN077986\t" $R/FL_TN_unified_exact_sequence_nr_map.tsv
head -2 $R/normalized_protein_id_map.tsv
ANN=$B/20_results/Figure2/07_new_figure/05_omic/1.final_counts/GO_annotation
ls $ANN
grep -i -P "oleosin|caleosin|REF/SRPP|rubber elongation|small rubber|SRPP|lipid droplet" $ANN/Africa_hap2/Africa_hap2.emapper.annotations | cut -f1,2,5,8,9,21 > out/ann_hits_africahap2.tsv
wc -l out/ann_hits_africahap2.tsv
which hmmscan hmmsearch mafft iqtree2 FastTree blastp diamond 2>&1 | head
ls ${CLUSTER_HOME}/miniconda3/envs/ 
