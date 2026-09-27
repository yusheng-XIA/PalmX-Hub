#!/usr/bin/env bash
set -euo pipefail

readonly RUN_ROOT="${ANALYSIS_DIR}/22_answer_reviews/00_ms/03_V3/03_figure3/13_comprehensive_primitive_ancestry_breeding_20260805/runs/RUN-COMP-BREEDING-DPM-20260805-001"
readonly RUN_RESULTS="${RUN_ROOT}/results"
readonly RUN_LOGS="${RUN_ROOT}/logs"
readonly RUN_PROVENANCE="${RUN_ROOT}/provenance"
readonly RUN_CHECKPOINTS="${RUN_ROOT}/checkpoints"
readonly RUN_SCRIPTS="${RUN_ROOT}/scripts"
readonly PYTHON="${DATA_DIR}/miniconda3/bin/python"

readonly CATALOG="${ANALYSIS_DIR}/22_answer_reviews/00_ms/03_V3/03_figure3/05_multiomics_integration/runs/RUN-MULTIOMICS-INTEGRATION-20260721-001/outputs/stage27_final_candidate_prioritization_attempt002/complete_ranked_25480_gene_candidate_catalog.tsv"
readonly MECHANISM="${ANALYSIS_DIR}/22_answer_reviews/00_ms/03_V3/03_figure3/05_multiomics_integration/runs/RUN-MULTIOMICS-HAPLOTYPE-HETEROSIS-SUPPLEMENT-20260723-001/outputs/stage8_mechanism_integration_attempt004/complete_ranked_gene_haplotype_metabolite_mechanism_catalog.tsv"
readonly PROTEIN="${ANALYSIS_DIR}/22_answer_reviews/00_ms/03_V3/03_figure3/07_proteomics_multiomics_20260724/07_genomic_breeding_update/runs/RUN-PROTEIN-GENOMIC-UPDATE-20260725-005/outputs/candidate_OAU_proteomics_evidence_updated.tsv"
readonly HETEROSIS="${ANALYSIS_DIR}/22_answer_reviews/00_ms/03_V3/03_figure3/11_diploid_32chrom_ancestry_breeding_20260805/results/TN_expression_heterosis_gene_summary.tsv"
readonly ASE="${ANALYSIS_DIR}/08_hifi_chromosome/10_TN_FL_pedigree_ancestry/runs/RUN-TNFL-FAST31-20260804-001/results/k31/ase/TN_FL_ancestry_ASE_gene_overlap.tsv"
readonly BINS="${ANALYSIS_DIR}/22_answer_reviews/00_ms/03_V3/03_figure3/11_diploid_32chrom_ancestry_breeding_20260805/results/diploid_ancestry_bins.tsv"
readonly GTF="${ANALYSIS_DIR}/22_answer_reviews/00_ms/03_V3/03_figure3/00_minipan/05_linear_hisat2/Africa_hap2.gtf"

assert_run_root() {
    [[ "$(pwd -P)" == "${RUN_ROOT}" ]] || {
        echo "[ERROR] Run from ${RUN_ROOT}; current=$(pwd -P)" >&2
        exit 64
    }
}
