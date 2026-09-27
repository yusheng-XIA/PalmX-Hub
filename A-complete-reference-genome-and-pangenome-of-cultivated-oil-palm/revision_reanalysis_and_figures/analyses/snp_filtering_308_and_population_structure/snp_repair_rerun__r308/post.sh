#!/bin/bash
# GWAS summaries, SHELL fine-mapping/conditional tests, SV far intervals, SD23 tables, ED9 data (after all scans)
set -euo pipefail
W=${CLUSTER_WORK}/snp_repair_r308/gwas; PY=python
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4
cd $W/cl
$PY s03_summarize.py norepro > ../logs/s03.log 2>&1
$PY s04_sv_complement.py > ../logs/s04.log 2>&1
$PY s04b_sv_far_M0.py > ../logs/s04b.log 2>&1
$PY s05_focal_r2.py > ../logs/s05.log 2>&1
$PY s06_regional_sv_r2.py > ../logs/s06.log 2>&1
$PY s20_shell_finemap.py > ../logs/s20.log 2>&1
LEAD=$($PY -c "import pandas as pd; d=pd.read_csv('../shell/out/finemap_shell_new.tsv',sep='\t'); print(int(d[d.trait=='Nut_weight_g'].snp_lead_pos.iloc[0]))")
echo LEAD=$LEAD >> ../logs/s20.log
LEAD=$LEAD $PY s21b_shell_alleles.py > ../logs/s21b.log 2>&1
$PY s23_sd23.py > ../logs/s23.log 2>&1
$PY s22_far_finemap.py > ../logs/s22.log 2>&1
$PY s30_ed9_data.py > ../logs/s30.log 2>&1
$PY s31_ed9a_beta.py > ../logs/s31.log 2>&1
echo DONE $(date -Is) >> ../logs/post.done
