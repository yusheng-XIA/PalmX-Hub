#!/bin/bash
# Step 2: build per-chromosome GT-only VCFs and run the original vcftools commands per chromosome
# (same parameters as 02_run_one_vcftools_task.sh: --window-pi 100000 --window-pi-step 50000;
#  --weir-fst-pop x2 --fst-window-size 100000 --fst-window-step 50000)
set -u
W=${CLUSTER_WORK}/fix_fig4c
VT=vcftools
P=$W/tmp/parts; C=$W/tmp/chrom; mkdir -p $C $W/out/fixed $W/out/raw $W/logs/vt
CHRS="chr14B chr04B chr06B chr05B chr02B chr12B chr03B chr09B chr08B chr13B chr15B chr11B chr01B chr16B chr10B chr07B"
{ echo '##fileformat=VCFv4.2'
  for c in $CHRS; do echo "##contig=<ID=$c>"; done
  echo '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">'
  cat $W/probe/header_chrom.txt | awk 'BEGIN{OFS="\t"}{$9="FORMAT"; print}'
} > $C/header.vcf
for c in $CHRS; do
  { cat $C/header.vcf; cat $(ls $P/part??_$c.body | sort); } > $C/fixed_$c.vcf
done
for c in chr14B chr16B; do
  { cat $C/header.vcf; cat $(ls $P/raw??_$c.body | sort); } > $C/raw_$c.vcf
done
echo "chrom files built $(date '+%F %T')" >> $W/logs/c2_status.txt
G=$W/groups
: > $W/tmp/jobs.txt
for set in fixed raw; do
  if [ $set = raw ]; then L="chr14B chr16B"; else L=$CHRS; fi
  for c in $L; do
    for p in 1 2 3 4; do
      echo "$VT --vcf $C/${set}_$c.vcf --keep $G/K4_Pop$p.txt --window-pi 100000 --window-pi-step 50000 --out $W/out/$set/${c}__K4_Pop${p}_100kb.pi > $W/logs/vt/${set}_${c}_P$p.log 2>&1" >> $W/tmp/jobs.txt
    done
    for pr in "1 2" "1 3" "1 4" "2 3" "2 4" "3 4"; do
      set -- $pr
      echo "$VT --vcf $C/${set}_$c.vcf --weir-fst-pop $G/K4_Pop$1.txt --weir-fst-pop $G/K4_Pop$2.txt --fst-window-size 100000 --fst-window-step 50000 --out $W/out/$set/${c}__K4_Pop$1_K4_Pop$2_100kb_fst > $W/logs/vt/${set}_${c}_F$1$2.log 2>&1" >> $W/tmp/jobs.txt
    done
  done
done
# validation jobs first
( grep ' raw_' $W/tmp/jobs.txt; grep -v ' raw_' $W/tmp/jobs.txt ) > $W/tmp/jobs_ordered.txt
nice -n 10 xargs -P 8 -I{} bash -c '{}' < $W/tmp/jobs_ordered.txt
echo "vcftools done $(date '+%F %T')" >> $W/logs/c2_status.txt
