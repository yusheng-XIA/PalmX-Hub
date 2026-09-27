#!/bin/bash
# Own FL/TN mesocarp libraries: raw FASTQ not readable, so primary read-1 records (incl. unmapped) are
# subsampled from the hisat2 BAM to ~2M reads; then the same exact 31-mer prefilter and minimap2 sr step as public runs.
W=${CLUSTER_WORK}/ole16_cis; cd $W; mkdir -p own_hits own_rq; export TMPDIR=$W/tmp
MM=minimap2
ST=samtools
B=${ANALYSIS_DIR}/22_answer_reviews/00_ms/03_V3/03_figure3/00_minipan/03_rnaseq_mapping/linear_hisat2_Africa_hap2
G=${B%/*}/pangraphrna_hisat2_graph
for i in 31 32 33 34 35 36 37 38 39 70 71 72 73 74 75 76 77 78; do
  r=Unknown_DJ272-001T00$i; BB=$B; [ -d $B/$r ] || BB=$G; bam=$(ls $BB/$r/*.sorted.bam $BB/$r/*.bam 2>/dev/null | head -1)
  [ -s own_rq/$r.sam ] && [ -s own_hits/$r.count ] && [ "$(cat own_hits/$r.count)" -gt 0 ] && continue
  echo "$r uses $BB"
  n1=$(grep " read1" $BB/$r/*flagstat* | awk '{print $1}')
  frac=$(awk -v n=$n1 'BEGIN{f=2000000/n; if(f>1)f=0.999999; printf "%.6f", f}')
  $ST view -@4 -f 64 -F 2304 -s 7${frac#0} $bam | awk -v OFS='\t' '{print $1,$10; c++} END{print c > "/dev/stderr"}' 2> own_hits/$r.count | LC_ALL=C grep -F -f pats31.txt > own_hits/$r.hits.tsv
  awk -F'\t' '{print ">"$1"\n"$2}' own_hits/$r.hits.tsv > own_rq/$r.fa
  $MM -ax sr -N 5 --secondary=yes -t 8 targets.fa own_rq/$r.fa 2>/dev/null | grep -v "^@" | cut -f1-6,12- > own_rq/$r.sam
  rm own_rq/$r.fa; echo "$r $n1 $frac $(cat own_hits/$r.count)"
done
echo ok > own_done
