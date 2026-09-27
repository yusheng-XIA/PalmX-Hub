G=${ANALYSIS_DIR}/08_hifi_chromosome/11_final_corrected_39genomes_annotations_20260812/attempt_20260812_01/02_annotations
O=${CLUSTER_WORK}/headline_pop/out
ls $G | head -50
awk -F'\t' '($1=="chr01"||$1=="chr01B") && $4<3400000 && $5>2900000' $G/African_hap2.gene.gff3 > $O/g3_region.gff3
awk -F'\t' '$3=="gene"||$3=="mRNA"' $O/g3_region.gff3 | cut -f1,3,4,5,7,9 | cut -c1-250
