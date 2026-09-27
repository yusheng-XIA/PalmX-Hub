#!/bin/bash
W=${CLUSTER_WORK}/snp_repair_rerun/fig4c
mkdir -p $W/out/pima_strict $W/out/pima_present
: > $W/logs/x4_jobs.txt
for set in present strict; do for f in $W/vcf/chr*.$set.vcf.gz; do c=$(basename $f .$set.vcf.gz)
 echo "python3 $W/cl/pi_ma.py $f $W/groups $W/out/pima_$set/$c > $W/logs/pima_${set}_$c.log 2>&1" >> $W/logs/x4_jobs.txt; done; done
nice -n 10 xargs -P ${1:-14} -I{} bash -c '{}' < $W/logs/x4_jobs.txt
echo "pima done $(date '+%F %T')" >> $W/logs/x4_done.txt
