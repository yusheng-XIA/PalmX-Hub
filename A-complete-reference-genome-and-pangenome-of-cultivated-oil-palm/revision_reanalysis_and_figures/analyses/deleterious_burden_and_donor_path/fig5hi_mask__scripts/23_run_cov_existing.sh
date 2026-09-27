#!/bin/bash
# coverage for the 29 donors whose original PAFs are readable (EG_* and bk_hap*)
W=${CLUSTER_WORK}/fig5hi_mask
M=${ANALYSIS_DIR}/14_pan_genome/04_SNP_calling/minimap2_results
PY=python3
mkdir -p $W/coverage
for s in EG_008 EG_015 EG_017 EG_025 EG_033 EG_035 EG_037 EG_041 EG_057 EG_058 EG_062 EG_065 EG_067 EG_071 EG_072 EG_075 EG_083 EG_086 EG_090 EG_095 EG_102 EG_107 EG_113 EG_146 EG_176 EG_183 EG_houke bk_hap1 bk_hap2; do
  [[ -s $W/coverage/$s.tsv ]] || $PY $W/scripts/20_coverage.py $s $M/$s/${s}_vs_Africa_hap2.sorted.paf $W/coverage/$s.tsv >> $W/logs/cov_existing.log 2>&1
done
echo ALLDONE >> $W/logs/cov_existing.log
