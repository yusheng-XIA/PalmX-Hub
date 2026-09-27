#!/bin/bash
# PanGenie 308-sample biallelic SV (>=50 bp) genotypes, before MAF filtering; called fraction >= 0.9 per group
set -euo pipefail
W=${CLUSTER_WORK}/enh_sn_pop
E=${CLUSTER_HOME}/miniconda3/envs/cs213/bin
IN=${ANALYSIS_DIR}/14_pan_genome/06_Minigraph/Pangenie/02_sv_combined/步骤四_过滤分类统计/sv.biallelic.sv50.tags.vcf.gz
G="ALL COMM NONC AFR HHG IDB SAEG SEAA SEAB K4P1 K4P2 K4P3 K4P4"
Q='%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n'
$E/bcftools annotate --threads 2 -x INFO,^FORMAT/GT -Ou $IN \
 | $E/bcftools query -H -f "$Q" \
 | python3 $W/cl/agg.py sv all $W/cl/groups308.txt - 0.9 $W/sv/sv_all
echo "sv done $(date '+%F %T')" >> $W/logs/sv_done.txt
