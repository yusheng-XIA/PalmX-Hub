#!/bin/bash
set -eo pipefail

TASK_INDEX=${1:?Usage: 02_run_one_vcftools_task.sh TASK_INDEX}

OUTDIR=${ANALYSIS_DIR}/22_answer_reviews/00_ms/03_V3/04_figure4/K3_K4_pi_fst_pca_groups
TASKS=${OUTDIR}/metadata/vcftools_tasks.tsv
VCF=${DATA_DIR2}/projects/1-oil_palm/06-population/oil_filter.vcf.recode.vcf
VCFTOOLS=${VCFTOOLS:-${DATA_DIR}/miniconda3/envs/oil_palm_map/bin/vcftools}

line_no=$((TASK_INDEX + 2))
task_line=$(awk -v n="${line_no}" 'NR == n {print}' "${TASKS}")
if [[ -z "${task_line}" ]]; then
  echo "No task line for TASK_INDEX=${TASK_INDEX}" >&2
  exit 1
fi

IFS=$'\t' read -r task k_label pop1 pop2 pop1_file pop2_file out_prefix <<< "${task_line}"
mkdir -p "$(dirname "${out_prefix}")"

echo "HOST=$(hostname)"
echo "DATE_START=$(date '+%F %T')"
echo "TASK_INDEX=${TASK_INDEX}"
echo "TASK=${task}"
echo "K=${k_label}"
echo "POP1=${pop1}"
echo "POP2=${pop2}"
echo "OUT_PREFIX=${out_prefix}"
echo "VCFTOOLS=${VCFTOOLS}"

if [[ "${task}" == "PI" ]]; then
  "${VCFTOOLS}" --vcf "${VCF}" \
    --keep "${pop1_file}" \
    --window-pi 100000 \
    --window-pi-step 50000 \
    --out "${out_prefix}"
elif [[ "${task}" == "FST" ]]; then
  "${VCFTOOLS}" --vcf "${VCF}" \
    --weir-fst-pop "${pop1_file}" \
    --weir-fst-pop "${pop2_file}" \
    --fst-window-size 100000 \
    --fst-window-step 50000 \
    --out "${out_prefix}"
else
  echo "Unknown task: ${task}" >&2
  exit 1
fi

echo "DATE_END=$(date '+%F %T')"
