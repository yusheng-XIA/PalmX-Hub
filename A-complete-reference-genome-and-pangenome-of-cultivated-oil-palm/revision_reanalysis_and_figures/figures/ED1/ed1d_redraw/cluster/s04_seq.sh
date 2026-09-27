#!/bin/bash
# ED1d: window sequences from the final FL assembly + identity check against the BAM reference fasta
set -uo pipefail
W=${CLUSTER_WORK}/ed1d_redraw
ST=samtools
F=${ANALYSIS_DIR}/22_answer_reviews/00_ms/03_V3/01_figure1/00_final_input_resources/05_syri_plotsr_chain12/results/00_standardized_genomes
R=${ANALYSIS_DIR}/01_seedless_results
mkdir -p $W/seq; cd $W/seq
$ST faidx $F/FL_American_hap1.chr16.fa chr12:1-1100000 | sed 's/^>.*/>chr12A_1_1100000/' > chr12A_head.fa
$ST faidx $F/FL_Africa_hap2.chr16.fa chr07:128400000-130504653 | sed 's/^>.*/>chr07B_128400000_130504653/' > chr07B_tail.fa
ls -la $R/American_hap1.fa* $R/Africa_hap2.fa* > ref_files.txt 2>&1
if [ -e $R/American_hap1.fa.fai ]; then
  $ST faidx $R/American_hap1.fa chr12A:1-1100000 | grep -v ">" | md5sum > md5_bamref.txt
  grep -v ">" chr12A_head.fa | md5sum >> md5_bamref.txt
fi
if [ -e $R/Africa_hap2.fa.fai ]; then
  $ST faidx $R/Africa_hap2.fa chr07B:128400000-130504653 | grep -v ">" | md5sum >> md5_bamref.txt
  grep -v ">" chr07B_tail.fa | md5sum >> md5_bamref.txt
fi
echo SEQDONE
