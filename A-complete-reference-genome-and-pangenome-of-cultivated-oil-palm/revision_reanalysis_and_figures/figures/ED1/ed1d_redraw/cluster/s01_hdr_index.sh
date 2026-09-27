#!/bin/bash
# ED1d redraw step 1: check BAM headers vs final FL assembly; build ONT .bai (user dir read-only -> bai in our dir)
set -uo pipefail
W=${CLUSTER_WORK}/ed1d_redraw
ST=samtools
P=${ANALYSIS_DIR}/01_seedless_results/12_evaluate
F=${ANALYSIS_DIR}/22_answer_reviews/00_ms/03_V3/01_figure1/00_final_input_resources/05_syri_plotsr_chain12/results/00_standardized_genomes
mkdir -p $W/bai $W/logs; cd $W
for s in 01_africa 02_american; do for b in hifi_remap ont_remap; do
  echo "== $s $b"; $ST view -H $P/$s/02_remapping/$b.sorted.bam | grep -E "^@SQ" | grep -E "chr(07|12)" ; $ST view -H $P/$s/02_remapping/$b.sorted.bam | grep -E "^@PG" | head -2 | cut -c1-300
done; done > hdr.txt 2>&1
grep -E "^chr(07|12)" $F/FL_Africa_hap2.chr16.fa.fai $F/FL_American_hap1.chr16.fa.fai >> hdr.txt
cp $F/FL_Africa_hap2.chromosome_map.tsv $F/FL_American_hap1.chromosome_map.tsv $F/FL_Africa_hap2.gaps.bed $F/FL_American_hap1.gaps.bed $W/ 2>/dev/null
for s in 01_africa 02_american; do
  ( $ST index -@ 4 $P/$s/02_remapping/ont_remap.sorted.bam $W/bai/${s}.ont_remap.bai && echo "done $s ont $(date)" ) >> $W/bai/index.log 2>&1 &
done
wait; echo ALLDONE >> $W/bai/index.log
