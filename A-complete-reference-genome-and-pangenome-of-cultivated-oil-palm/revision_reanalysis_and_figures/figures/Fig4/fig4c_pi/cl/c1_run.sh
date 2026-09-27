#!/bin/bash
# Step 1: split the source VCF (read-only, user) into GT-only per-chromosome files, 8 parallel byte ranges
V=${DATA_DIR2}/projects/1-oil_palm/05-GWAS/oil_filter.vcf.recode.vcf
W=${CLUSTER_WORK}/fix_fig4c
mkdir -p $W/tmp/parts $W/logs
NW=8
echo "start $(date '+%F %T')" > $W/logs/c1_status.txt
for i in $(seq 0 $((NW-1))); do
  nice -n 10 python3 $W/c1_gtonly.py $V $i $NW $W/tmp/parts chr16B,chr14B > $W/logs/c1_w$i.log 2>&1 &
done
wait
echo "end $(date '+%F %T')" >> $W/logs/c1_status.txt
