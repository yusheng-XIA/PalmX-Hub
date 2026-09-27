#!/bin/bash
# concat chunks per chromosome (skips a set already concatenated+indexed); validate sort order and sample count
set -uo pipefail
export OPENBLAS_NUM_THREADS=1
E=${CLUSTER_HOME}/miniconda3/envs/genomics_a2/bin
W=${CLUSTER_WORK}/snp_repair_r308
c=$1; mkdir -p $W/vcf
for s in diversity population; do
  F=$W/vcf/$c.${s}_snps.corrected.vcf.gz
  if [ ! -s $F.tbi ]; then
    ls $W/filter/$c/${c}_??.$s.vcf.gz | sort > $W/vcf/$c.$s.list
    $E/bcftools concat --threads 4 -f $W/vcf/$c.$s.list -Oz -o $F || { echo "FAIL concat $c $s"; exit 1; }
    $E/bcftools index -t --threads 4 $F || exit 1
  fi
  n=$($E/bcftools index -n $F); ns=$($E/bcftools query -l $F | wc -l)
  read unsorted first last <<< $($E/bcftools query -f "%POS\n" $F | awk "NR==1{f=\$1} \$1<p{b++}{p=\$1}END{print b+0, f, p}")
  echo -e "$c\t$s\t$n\t$ns\t$first\t$last\tunsorted=$unsorted" >> $W/vcf/concat_summary.tsv
  sha256sum $F >> $W/vcf/sha256.txt
done
