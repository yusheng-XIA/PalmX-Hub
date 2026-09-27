#!/bin/bash
# per repaired chromosome: wait for concat, then 308 subset (s2), Fig.4c extraction, SD16 counts, G6 genotypes
c=$1
W=${CLUSTER_WORK}/snp_repair_r308
until grep -qP "^$c\tdiversity\t" $W/vcf/concat_summary.tsv 2>/dev/null; do sleep 30; done
bash $W/cl/s2_sub308.sh $c > $W/logs/s2_$c.log 2>&1 || echo FAIL s2 $c >> $W/logs/chain_fail.txt &
until grep -qP "^$c\tpopulation\t" $W/vcf/concat_summary.tsv 2>/dev/null; do sleep 30; done
bash $W/fig4c/cl/x2_extract.sh $c > $W/fig4c/logs/x2_$c.log 2>&1 || echo FAIL x2 $c >> $W/logs/chain_fail.txt &
mkdir -p $W/sd16/snp; bash $W/sd16/cl/snp1.sh $c > $W/sd16/logs/snp_$c.log 2>&1 || echo FAIL sd16 $c >> $W/logs/chain_fail.txt &
bash $W/g6/cl/g6_b_geno.sh $c > $W/g6/geno_$c.log 2>&1 || echo FAIL g6 $c >> $W/logs/chain_fail.txt &
wait
echo "$c chain done $(date -Is)" >> $W/logs/chain_done.txt
