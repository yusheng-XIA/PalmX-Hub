B=${ANALYSIS_DIR}/21_MS/06_result/dSVs
ls -la $B/input | head; head -3 $B/input/Africa_hap2.fa.fai
awk -F'\t' '$3=="CDS"' $B/input/Africa_hap2.EVM.gff3 | head -2
awk -F'\t' '$3=="mRNA"' $B/input/Africa_hap2.EVM.gff3 | wc -l
awk -F'\t' '$3=="gene"' $B/input/Africa_hap2.EVM.gff3 | wc -l
ls -la $B/results-8.9/06_dSNP_phoenix_polarity/alignments/
head -3 ${CLUSTER_WORK}/fig4c_final/meta/s308.txt; wc -l ${CLUSTER_WORK}/fig4c_final/meta/s308.txt
