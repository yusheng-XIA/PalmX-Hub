#!/bin/bash
W=${CLUSTER_WORK}/fig4c_final
O=${CLUSTER_WORK}/fix_fig4c/tmp/chrom
mkdir -p $W/out/pima_oldvcf
ls $O/fixed_chr*.vcf | xargs -P 14 -I{} bash -c 'c=$(basename {} .vcf); c=${c#fixed_}; nice -n 10 python3 '$W'/cl/pi_ma.py {} '$W'/groups '$W'/out/pima_oldvcf/$c'
echo "old done $(date '+%F %T')" >> $W/logs/x5_done.txt
