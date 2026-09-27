#!/bin/bash
set -euo pipefail

TIDK=${DATA_DIR}/miniconda3/envs/haphic/bin/tidk
BASE=${ANALYSIS_DIR}/05_GWAS/00_analysis/00_data/public_data
OUT=$BASE/telomere_tidk
LOG=$OUT/run_tidk.log

{
echo "### run_tidk.sh started: $(date -Is)"
echo "### host: $(hostname)"
echo "### tidk version: $($TIDK --version)"
echo "### input checksums:"
md5sum "$BASE/EG11_chromosomes.fa" "$BASE/EO12_chromosomes.fa"
echo "### command: $TIDK search -s TTTAGGG -w 10000 -o EG11 -d $OUT/EG11 --log $BASE/EG11_chromosomes.fa"
echo
} | tee "$LOG"

$TIDK search -s TTTAGGG -w 10000 -o EG11 -d "$OUT/EG11" --log "$BASE/EG11_chromosomes.fa" 2>&1 | tee -a "$LOG"
$TIDK search -s TTTAGGG -w 10000 -o EO12 -d "$OUT/EO12" --log "$BASE/EO12_chromosomes.fa" 2>&1 | tee -a "$LOG"

{
echo "### outputs:"
md5sum "$OUT/EG11/"* "$OUT/EO12/"* 2>/dev/null || true
echo "### run_tidk.sh finished: $(date -Is)"
} | tee -a "$LOG"
