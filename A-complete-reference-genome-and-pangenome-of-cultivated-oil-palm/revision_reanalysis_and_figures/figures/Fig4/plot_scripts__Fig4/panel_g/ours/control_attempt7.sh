#!/usr/bin/env bash
# Environment-equivalence control: rerun the attempt7 GeneTribe + classification step on the
# unchanged attempt2 dura inputs with env fig4g_g2, and compare with attempt7 (RBH 31768; 20671/11077/20).
set -euo pipefail
PAIR_ID=dura_new; THREADS="${1:-14}"
RUN=${ANALYSIS_DIR}/22_answer_reviews/00_ms/03_V3/04_figure4/Fig4_d_i_pan39_material33_redraw_20260808
SRC=${RUN}/work/g_allele/${PAIR_ID}/attempt2
GENETRIBE=${DATA_DIR}/youzong/software/genetribe/genetribe
export PATH=${CLUSTER_HOME}/miniconda3/envs/fig4g_g2/bin:$PATH
WD=${CLUSTER_WORK}/G2/control_dura_attempt7
mkdir -p ${WD}/GeneTribe ${WD}/final
exec > >(tee -a ${WD}/pipeline.log) 2>&1
echo "START=$(date --iso-8601=seconds) HOST=$(hostname)"
sha256sum ${SRC}/p.pep ${SRC}/a.pep ${SRC}/GeneTribe/pctg.bed ${SRC}/GeneTribe/actg.bed
cd ${WD}/GeneTribe
awk -v p='P__' '/^>/{sub(/^>/, ">" p)} {print}' ${SRC}/p.pep > pctg.fa
awk -v p='A__' '/^>/{sub(/^>/, ">" p)} {print}' ${SRC}/a.pep > actg.fa
awk -v OFS='\t' -v p='P__' '{$4=p $4; print}' ${SRC}/GeneTribe/pctg.bed > pctg.bed
awk -v OFS='\t' -v p='A__' '{$4=p $4; print}' ${SRC}/GeneTribe/actg.bed > actg.bed
echo N > pctg.chrlist; echo N > actg.chrlist
${GENETRIBE} core -l pctg -f actg -s : -n ${THREADS} > genetribe.log 2>&1
cd ${WD}/final
awk '{p=$1; a=$2; sub(/^P__/, "", p); sub(/^A__/, "", a); print a "\t" p}' ../GeneTribe/pctg_actg.RBH > rbh_pairs.tsv
awk 'NR==FNR {keys[$1 SUBSEP $2]=1; next} {key=$1 SUBSEP $2; if(key in keys && !seen[key]++) print}' rbh_pairs.tsv ${SRC}/gene_blast_gene/a-gene_to_p-gene_blast.out > rbh_blast.out
awk '$3==100 && $4==$11 {print $1"\t"$2}' rbh_blast.out | sort -u > homozygous.list
awk '$3<100 || $4!=$11 {print $1"\t"$2}' rbh_blast.out | sort -u > heterozygous.list
awk 'NR==FNR {hit[$1 SUBSEP $2]=1; next} {if(!(($1 SUBSEP $2) in hit)) print $1"\t"$2}' rbh_blast.out rbh_pairs.tsv > no_blast_hit.list
echo "CONTROL RBH=$(wc -l < rbh_pairs.tsv) HOMO=$(wc -l < homozygous.list) HETERO=$(wc -l < heterozygous.list) NOHIT=$(wc -l < no_blast_hit.list)"
A7=${RUN}/work/g_allele/${PAIR_ID}/attempt7_uniqueid/final/rbh_pairs.tsv
echo "PAIRS_IDENTICAL_TO_ATTEMPT7=$(cmp -s <(sort rbh_pairs.tsv) <(sort ${A7}) && echo yes || echo no) diff_lines=$(diff <(sort rbh_pairs.tsv) <(sort ${A7}) | grep -c '^[<>]')"
echo "DONE=$(date --iso-8601=seconds)"; touch ${WD}/SUCCESS
