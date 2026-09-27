#!/usr/bin/env bash
set -euo pipefail

run_root="${ANALYSIS_DIR}/22_answer_reviews/00_ms/03_V3/03_figure3/11_diploid_32chrom_ancestry_breeding_20260805"
ancestry_root="${ANALYSIS_DIR}/08_hifi_chromosome/10_TN_FL_pedigree_ancestry/runs/RUN-TNFL-FAST31-20260804-001/results/k31"
fig3_root="${ANALYSIS_DIR}/22_answer_reviews/00_ms/03_V3/03_figure3"
stage_root="$(mktemp -d /tmp/diploid32_20260805.XXXXXX)"

[[ ! -e "${run_root}/PASS" ]] || { echo "[ERROR] completed run exists" >&2; exit 70; }
date --iso-8601=seconds > "${run_root}/provenance/run_started_at.txt"
${DATA_DIR}/miniconda3/bin/python -B "${run_root}/scripts/build_diploid_ancestry_design.py" \
  --ancestry-windows "${ancestry_root}/reports/four_haplotype_ancestry_windows.100kb.tsv" \
  --ase-overlap "${ancestry_root}/ase/TN_FL_ancestry_ASE_gene_overlap.tsv" \
  --lipid-candidates "${fig3_root}/08_Fig3_multitrait_ASE_redesign_20260725/04_Fig3d_candidate_ASE_enrichment/source_data/Fig3d_frozen_Fig2_candidates.tsv" \
  --mesocarp-candidates "${fig3_root}/08_Fig3_multitrait_ASE_redesign_20260725/05_Fig3e_mesocarp_thickening/source_data/Fig3e_thickening_candidate_summary.tsv" \
  --rancidity-candidates "${fig3_root}/10_Figure3_ASE_complementarity_flat_20260727/Fig3g_TN_NS_direct_rancidity_candidates.tsv" \
  --heterosis-gene-stage "${fig3_root}/01_ASE/00_shared/runs/RUN-TN-EXPR-HETEROSIS-V2-001/output/gene_stage_heterosis.tsv.gz" \
  --trait-gene-catalog "${fig3_root}/01_ASE/00_shared/runs/RUN-ASE-HETEROSIS-DOWNSTREAM-V2-001/output/trait_gene_catalog.tsv" \
  --run-root "${stage_root}" --bins 100 --long-block-percent 5 \
  > "${run_root}/logs/stdout.log" 2> "${run_root}/logs/stderr.log"
cp "${stage_root}/results/"* "${run_root}/results/"
cp "${stage_root}/figures/"* "${run_root}/figures/"
cp "${stage_root}/provenance/"* "${run_root}/provenance/"
cp "${stage_root}/PASS" "${run_root}/PASS"
printf '%s\n' "${stage_root}" > "${run_root}/provenance/local_stage_root.txt"
find "${run_root}/results" "${run_root}/figures" -maxdepth 1 -type f -print0 \
  | sort -z | xargs -0 sha256sum > "${run_root}/provenance/accepted_outputs.sha256"
date --iso-8601=seconds > "${run_root}/provenance/run_finished_at.txt"
