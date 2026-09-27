#!/bin/bash
# Build .bai for FL DNA BAMs (read-only in user dir) into our work dir.
set -uo pipefail
W=${CLUSTER_WORK}/enh_A
ST=samtools
P=${ANALYSIS_DIR}/01_seedless_results/12_evaluate
mkdir -p $W/bai; cd $W/bai
for s in 01_africa 02_american; do
  for b in seedless3-10 hifi_remap; do
    ( $ST index -@ 4 $P/$s/02_remapping/$b.sorted.bam $W/bai/${s}.${b}.bai && echo "done $s $b $(date)" ) >> $W/bai/index.log 2>&1 &
  done
done
wait
echo ALLDONE >> $W/bai/index.log
