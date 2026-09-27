#!/bin/bash
# Fig. 2 divergence times on the final MCMCTree run (user job 847506, 2026-06-17).
# Run from the scratchpad root (same convention as final_rebuild.sh). Not run by the candidate build.
#   bash fix/fig2_time/fix_sd_fig2_time.sh          -> Option A (time axis only; authors' CAFE5 numbers kept)
#   bash fix/fig2_time/fix_sd_fig2_time.sh B        -> Option B (also CAFE5 re-run on the final-run tree: Fig2b, SF2b)
set -euo pipefail
D=fix/fig2_time/sd
python3 deliver_build/fix_sd_add.py deliver/Source_Data/Source_Data.xlsx \
  "Fig2a_divergence=$D/Fig2a_divergence.tsv|replace=Fig2a_divergence|desc=MCMCTree posterior mean ages and 95% HPD intervals (Ma) of the five Fig. 2a nodes (PAML 4.10.9, seven calibrations)" \
  "Fig2b_node_ages=$D/Fig2b_node_ages.tsv|after=Fig2b_gene_families|fig=Figure 2|panel=b|desc=Node ages (Ma) of the time-calibrated tree in Fig. 2b (MCMCTree posterior means and 95% HPD intervals; CAFE5 node ids as in Fig2b_gene_families)"
if [ "${1:-A}" = "B" ]; then
  python3 deliver_build/fix_sd_add.py deliver/Source_Data/Source_Data.xlsx \
    "Fig2b_gene_families=$D/optionB/Fig2b_gene_families.tsv|replace=Fig2b_gene_families|desc=CAFE5 (gamma model, 3 rate categories) expanded, contracted and unchanged families per branch on the MCMCTree time-calibrated tree of 12 monocot genomes" \
    "SF2b_ancestral_changes=$D/optionB/SF2b_ancestral_changes.tsv|replace=SF2b_ancestral_changes|desc=CAFE5 gains/losses at the oil palm ancestral node (re-run on the MCMCTree time-calibrated tree)"
fi
