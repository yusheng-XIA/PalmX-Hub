#!/bin/bash
# Rebuild Figure4_candidate.pdf from the submitted Figure4.pdf (read-only) - run from fix/fig4/scripts
# 00_contig_nx_ragtag.py runs on ${COMPUTE_HOST} (output copied to ../src/nx_ragtag_full.tsv)
set -euo pipefail
cd "$(dirname "$0")"
cp ../src/Figure4_submitted_copy.pdf ../steps/step0.pdf
cmp ../steps/step0.pdf ../../../deliver/Main_Figures_revised/Figure4.pdf
python3 01_fig4d_contig_nx.py
python3 02_fig4i_te_density.py
python3 03_fig4b_pca_labels.py
python3 04_fig4c_mean_fst.py
python3 05_fig4k_class_type.py
prev=../steps/step0.pdf
for s in step1_4d step2_4i step3_4b step4_4c step5_4k; do python3 step_check.py $prev ../steps/$s.pdf; prev=../steps/$s.pdf; done
python3 -c "import fitz; fitz.open('../steps/step5_4k.pdf').save('../Figure4_candidate.pdf', garbage=3, deflate=True)"
python3 06_verify_render.py > /dev/null
echo "done: ../Figure4_candidate.pdf ; checks in ../compare/verify_report.json"
