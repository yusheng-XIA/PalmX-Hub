#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 3 ]]; then
  echo "Usage: $0 PAIR_ID SOURCE_ATTEMPT OUTPUT_ATTEMPT" >&2
  exit 2
fi

PAIR_ID="$1"
SOURCE_ATTEMPT="$2"
OUTPUT_ATTEMPT="$3"
RUN="${ANALYSIS_DIR}/22_answer_reviews/00_ms/03_V3/04_figure4/Fig4_d_i_pan39_material33_redraw_20260808"
GENETRIBE="${DATA_DIR}/youzong/software/genetribe/genetribe"
THREADS="${THREADS:-40}"
SRC="${RUN}/work/g_allele/${PAIR_ID}/${SOURCE_ATTEMPT}"
WD="${RUN}/work/g_allele/${PAIR_ID}/${OUTPUT_ATTEMPT}"
OUT="${RUN}/tables/g_new_pairs/${PAIR_ID}.allele_summary.${OUTPUT_ATTEMPT}.tsv"
LOG="${RUN}/logs/g_${PAIR_ID}_${OUTPUT_ATTEMPT}.pipeline.log"
PROV="${RUN}/provenance"
SCRATCH="$(mktemp -d "/tmp/fig4g-${PAIR_ID}-${OUTPUT_ATTEMPT}.XXXXXX")"

cleanup() {
  if [[ "${SCRATCH}" == "/tmp/fig4g-${PAIR_ID}-${OUTPUT_ATTEMPT}."* && -d "${SCRATCH}" ]]; then
    rm -rf -- "${SCRATCH}"
  fi
}
trap cleanup EXIT

if [[ -e "${WD}" || -e "${OUT}" ]]; then
  echo "FATAL attempt output already exists: ${WD} or ${OUT}" >&2
  exit 3
fi

for path in \
  "${SRC}/p.pep" "${SRC}/a.pep" \
  "${SRC}/GeneTribe/pctg.bed" "${SRC}/GeneTribe/actg.bed" \
  "${SRC}/gene_blast_gene/a-gene_to_p-gene_blast.out" \
  "${RUN}/tables/g_new_pairs/${PAIR_ID}.allele_summary.tsv"; do
  [[ -s "${path}" ]] || { echo "FATAL missing input ${path}" >&2; exit 11; }
done

mkdir -p "${WD}" "${RUN}/logs" "${RUN}/tables/g_new_pairs" "${PROV}"
exec > >(tee -a "${LOG}") 2>&1

echo "PLAN_ID=FIG4G-UNIQUEID-20260810-01"
echo "START=$(date --iso-8601=seconds) HOST=$(hostname) PAIR=${PAIR_ID} SOURCE=${SOURCE_ATTEMPT} OUTPUT=${OUTPUT_ATTEMPT} THREADS=${THREADS} SCRATCH=${SCRATCH}"
sha256sum "${SRC}/p.pep" "${SRC}/a.pep" "${SRC}/GeneTribe/pctg.bed" "${SRC}/GeneTribe/actg.bed" > "${WD}/input_small_files.sha256"
stat -c '%n\t%s\t%Y' "${SRC}/gene_blast_gene/a-gene_to_p-gene_blast.out" > "${WD}/large_input_identity.tsv"

set +u
source ${DATA_DIR}/miniconda3/etc/profile.d/conda.sh
conda activate allele_env
set -u
echo "GENETRIBE_VERSION=$(${GENETRIBE} --version 2>&1 | tr '\n' ' ' || true)"
blastp -version | head -1

mkdir -p "${SCRATCH}/GeneTribe"
awk -v p='P__' '/^>/{sub(/^>/, ">" p)} {print}' "${SRC}/p.pep" > "${SCRATCH}/GeneTribe/pctg.fa"
awk -v p='A__' '/^>/{sub(/^>/, ">" p)} {print}' "${SRC}/a.pep" > "${SCRATCH}/GeneTribe/actg.fa"
awk -v OFS='\t' -v p='P__' '{$4=p $4; print}' "${SRC}/GeneTribe/pctg.bed" > "${SCRATCH}/GeneTribe/pctg.bed"
awk -v OFS='\t' -v p='A__' '{$4=p $4; print}' "${SRC}/GeneTribe/actg.bed" > "${SCRATCH}/GeneTribe/actg.bed"
printf 'N\n' > "${SCRATCH}/GeneTribe/pctg.chrlist"
printf 'N\n' > "${SCRATCH}/GeneTribe/actg.chrlist"

overlap="$(comm -12 \
  <(awk '{print $4}' "${SCRATCH}/GeneTribe/pctg.bed" | sort -u) \
  <(awk '{print $4}' "${SCRATCH}/GeneTribe/actg.bed" | sort -u) | wc -l)"
[[ "${overlap}" -eq 0 ]] || { echo "FATAL prefixed ID overlap=${overlap}"; exit 12; }
echo "PREFIXED_ID_OVERLAP=${overlap}"

cd "${SCRATCH}/GeneTribe"
"${GENETRIBE}" core -l pctg -f actg -s : -n "${THREADS}" 2>&1 | tee genetribe.log
[[ -s pctg_actg.RBH ]] || { echo "FATAL empty corrected RBH"; exit 13; }

mkdir -p "${SCRATCH}/final" "${WD}/GeneTribe" "${WD}/final"
cp -a pctg_actg.RBH pctg_actg.SBH actg_pctg.SBH pctg_actg.one2one actg_pctg.one2one \
  pctg_actg.singleton actg_pctg.singleton genetribe.log "${WD}/GeneTribe/"

cd "${SCRATCH}/final"
awk '{p=$1; a=$2; sub(/^P__/, "", p); sub(/^A__/, "", a); print a "\t" p}' \
  "${SCRATCH}/GeneTribe/pctg_actg.RBH" > rbh_pairs.tsv
awk 'NR==FNR {keys[$1 SUBSEP $2]=1; next} {key=$1 SUBSEP $2; if(key in keys && !seen[key]++) print}' \
  rbh_pairs.tsv "${SRC}/gene_blast_gene/a-gene_to_p-gene_blast.out" > rbh_blast.out
awk '$3==100 && $4==$11 {print $1"\t"$2}' rbh_blast.out | sort -u > homozygous.list
awk '$3<100 || $4!=$11 {print $1"\t"$2}' rbh_blast.out | sort -u > heterozygous.list
awk 'NR==FNR {hit[$1 SUBSEP $2]=1; next} {if(!(($1 SUBSEP $2) in hit)) print $1"\t"$2}' \
  rbh_blast.out rbh_pairs.tsv > no_blast_hit.list

rbh_n="$(wc -l < rbh_pairs.tsv)"
homo_n="$(wc -l < homozygous.list)"
hetero_n="$(wc -l < heterozygous.list)"
nohit_n="$(wc -l < no_blast_hit.list)"
[[ $((homo_n + hetero_n + nohit_n)) -eq "${rbh_n}" ]] || { echo "FATAL RBH accounting mismatch"; exit 14; }

source_summary="${RUN}/tables/g_new_pairs/${PAIR_ID}.allele_summary.tsv"
read -r p_gene a_gene cds_n gmap_n inter_n < <(awk -F '\t' 'NR==2{print $2, $3, $5, $6, $7}' "${source_summary}")
printf 'Sample\tp_gene\ta_gene\tRBH\tCDS_unique\tGMAP_unique\tUnique_inter\tHomozygous\tHeterozygous\tNo_hit\n' > "${OUT}"
printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
  "${PAIR_ID}" "${p_gene}" "${a_gene}" "${rbh_n}" "${cds_n}" "${gmap_n}" "${inter_n}" \
  "${homo_n}" "${hetero_n}" "${nohit_n}" >> "${OUT}"

cp -a rbh_pairs.tsv rbh_blast.out homozygous.list heterozygous.list no_blast_hit.list "${WD}/final/"
nohit_rate="$(awk -v n="${nohit_n}" -v r="${rbh_n}" 'BEGIN{printf "%.6f", 100*n/r}')"
echo "RBH=${rbh_n} HOMO=${homo_n} HETERO=${hetero_n} NO_HIT=${nohit_n} NO_HIT_RATE_PCT=${nohit_rate}"
sha256sum "${OUT}" > "${PROV}/${PAIR_ID}.${OUTPUT_ATTEMPT}.sha256"
touch "${PROV}/${PAIR_ID}.${OUTPUT_ATTEMPT}.SUCCESS"
echo "DONE=$(date --iso-8601=seconds) OUT=${OUT}"
