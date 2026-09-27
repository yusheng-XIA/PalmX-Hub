#!/bin/bash
# ED1d: confirm BAM reference fasta (01_seedless_results/*.fa) == final standardized FL assembly in the windows
set -uo pipefail
W=${CLUSTER_WORK}/ed1d_redraw
ST=samtools
F=${ANALYSIS_DIR}/22_answer_reviews/00_ms/03_V3/01_figure1/00_final_input_resources/05_syri_plotsr_chain12/results/00_standardized_genomes
R=${ANALYSIS_DIR}/01_seedless_results
cd $W/seq
[ -e American_hap1.fa.fai ] || $ST faidx --fai-idx $W/seq/American_hap1.fa.fai $R/American_hap1.fa
s(){ grep -v ">" | tr -d '\n' | tr a-z A-Z | md5sum | cut -d" " -f1; }
{
echo -e "chr12A:1-1100000\tbamref\t$($ST faidx --fai-idx $W/seq/American_hap1.fa.fai $R/American_hap1.fa chr12A:1-1100000 | s)"
echo -e "chr12A:1-1100000\tfinal\t$(cat chr12A_head.fa | s)"
echo -e "chr07B:128400000-130504653\tbamref\t$($ST faidx $R/Africa_hap2.fa chr07B:128400000-130504653 | s)"
echo -e "chr07B:128400000-130504653\tfinal\t$(cat chr07B_tail.fa | s)"
} > refcheck.tsv
cat refcheck.tsv
