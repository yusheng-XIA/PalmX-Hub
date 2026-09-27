#!/bin/bash
# ED1d: re-derive the draft/added junction on the final assembly by aligning the terminal 150 kb of the draft
# (Verkko) chromosome ends to the final chromosome ends (minimap2 asm5 + asm20; all hits kept).
set -uo pipefail
W=${CLUSTER_WORK}/ed1d_redraw
ST=samtools
MM=minimap2
D=${ANALYSIS_DIR}/01_seedless_results1/17_verkko_correct/chrom
mkdir -p $W/junction; cd $W/junction
# draft ends (faidx index written to our dir)
for f in hap1_10_all hap2_29; do [ -e $f.fai ] || $ST faidx --fai-idx $W/junction/$f.fai $D/$f.fa; done
cut -f1,2 hap1_10_all.fai hap2_29.fai > draft_lengths.tsv
n1=$(head -1 hap1_10_all.fai | cut -f1); n2=$(head -1 hap2_29.fai | cut -f1); L2=$(head -1 hap2_29.fai | cut -f2)
$ST faidx --fai-idx hap1_10_all.fai $D/hap1_10_all.fa $n1:1-150000 | sed "1s/.*/>draft_hap1_10_1_150000/" > draft_hap1_10_head.fa
$ST faidx --fai-idx hap2_29.fai $D/hap2_29.fa $n2:$((L2-149999))-$L2 | sed "1s/.*/>draft_hap2_29_$((L2-149999))_$L2/" > draft_hap2_29_tail.fa
for x in asm5 asm20; do
  $MM -t 4 -c -x $x --secondary=yes -N 20 $W/seq/chr12A_head.fa draft_hap1_10_head.fa > chr12A.$x.paf 2>/dev/null
  $MM -t 4 -c -x $x --secondary=yes -N 20 $W/seq/chr07B_tail.fa draft_hap2_29_tail.fa > chr07B.$x.paf 2>/dev/null
done
echo JDONE
