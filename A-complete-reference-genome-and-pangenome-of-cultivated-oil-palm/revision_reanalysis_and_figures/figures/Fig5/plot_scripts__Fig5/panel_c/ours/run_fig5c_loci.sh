#!/bin/bash
set -euo pipefail
W=${CLUSTER_WORK}/fix_fig5c
mkdir -p $W/tmp; cd $W
BT=bedtools
PY=python
G=${ANALYSIS_DIR}/08_hifi_chromosome/11_final_corrected_39genomes_annotations_20260812/attempt_20260812_01/02_annotations/African_hap2.gene.gff3
CAT=${ANALYSIS_DIR}/22_answer_reviews/00_ms/03_V3/05_figure/08_SV_hap39_Figure5abcd_20260808/results/final_catalog/oilpalm_hap38.final_highconfidence_sv.tsv
export LC_ALL=C
# gene / CDS beds with gene IDs (chrNN -> chrNNB)
awk -F'\t' 'BEGIN{OFS="\t"} !/^#/ && NF>=9 {if($1~/^chr[0-9][0-9]$/)$1=$1"B";
  id="";par=""; n=split($9,a,";"); for(i=1;i<=n;i++){if(a[i]~/^ID=/)id=substr(a[i],4); if(a[i]~/^Parent=/)par=substr(a[i],8)}
  if($3=="gene") print $1,$4-1,$5,id,".",$7 > "gene_id.bed";
  else if($3=="mRNA") print id,par > "mrna2gene.txt";
  else if($3=="CDS") print $1,$4-1,$5,par > "cds_raw.bed"}' $G
awk 'BEGIN{OFS="\t"} NR==FNR{m[$1]=$2;next}{split($4,p,","); g=(p[1] in m)?m[p[1]]:p[1]; print $1,$2,$3,g}' FS=' ' mrna2gene.txt FS='\t' cds_raw.bed | sort -k1,1 -k2,2n -u > cds_id.bed
sort -k1,1 -k2,2n gene_id.bed -o gene_id.bed
awk -F'\t' 'BEGIN{OFS="\t"}{if($6=="-"){s=$3;e=$3+2000}else{s=$2-2000;if(s<0)s=0;e=$2} print $1,s,e,$4}' gene_id.bed | sort -k1,1 -k2,2n > up_id.bed
awk -F'\t' 'BEGIN{OFS="\t"}{if($6=="-"){s=$2-2000;if(s<0)s=0;e=$2}else{s=$3;e=$3+2000} print $1,s,e,$4}' gene_id.bed | sort -k1,1 -k2,2n > down_id.bed
# SV context intervals: INS = breakpoint +/-50 bp ([Start-1-50, Start+50]); others [Start-1, End] (End<=Start-1 -> 1 bp at Start)
awk -F'\t' 'BEGIN{OFS="\t"} NR>1{sub(/\r$/,""); s=$3-1; e=$4; if($5=="INS"){s=$3-1-50; if(s<0)s=0; e=$3+50} else if(e<=s){e=$3} print $2,s,e,$1,$5,$7,$3,$4,"r"NR-1}' $CAT | sort -k1,1 -k2,2n > sv.bed
wc -l sv.bed gene_id.bed cds_id.bed
awk -F"\t" 'BEGIN{OFS="\t"}{print $1,$2,$3,$9}' sv.bed > sv4.bed
$BT intersect -wa -wb -nonamecheck -a sv4.bed -b cds_id.bed | cut -f4,8 > hit_cds.txt
$BT intersect -wa -wb -nonamecheck -a sv4.bed -b gene_id.bed | cut -f4,8 > hit_gene.txt
$BT intersect -wa -wb -nonamecheck -a sv4.bed -b up_id.bed | cut -f4,8 > hit_up.txt
$BT intersect -wa -wb -nonamecheck -a sv4.bed -b down_id.bed | cut -f4,8 > hit_down.txt
$BT closest -d -t all -nonamecheck -a sv4.bed -b gene_id.bed | awk -F'\t' 'BEGIN{OFS="\t"}{print $4,$8,$11}' > nearest.txt
$PY ${CLUSTER_WORK}/fix_fig5c/build_table.py
echo DONE
