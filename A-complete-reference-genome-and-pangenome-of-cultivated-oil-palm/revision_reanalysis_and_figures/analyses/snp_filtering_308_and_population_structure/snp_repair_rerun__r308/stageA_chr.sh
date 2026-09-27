#!/bin/bash
# per chromosome: concat 308-filtered chunks -> vcf/{div308,pop308}.bcf; PLINK structure/GWAS sets; Fig. 4c SNP set;
# G6 coding genotypes; SD16 group-frequency counts.   usage: stageA_chr.sh CHR
set -euo pipefail
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
c=$1
F=${CLUSTER_WORK}/snp_repair_rerun/r308/filter
W=${CLUSTER_WORK}/snp_repair_r308
E=${CLUSTER_HOME}/miniconda3/envs/genomics_a2/bin
export TMPDIR=$W/tmp; mkdir -p $W/vcf $W/tmp $W/logs $W/g6/geno $W/sd16/snp $W/sd16/logs $W/fig4c/logs
for s in div308 pop308; do
  ls $F/$c/${c}_??.$s.bcf | sort > $W/vcf/$c.$s.list
  $E/bcftools concat --threads 2 -f $W/vcf/$c.$s.list -Ob -o $W/vcf/$c.$s.bcf
  $E/bcftools index --threads 2 $W/vcf/$c.$s.bcf
  echo -e "$c\t$s\t$($E/bcftools index -n $W/vcf/$c.$s.bcf)" >> $W/vcf/counts.tsv
done
bash $W/cl/s2_sub308.sh $c > $W/logs/s2_$c.log 2>&1
bash $W/fig4c/cl/x2_extract.sh $c > $W/fig4c/logs/x2_$c.log 2>&1
bash $W/g6/cl/g6_b_geno.sh $c > $W/logs/g6geno_$c.log 2>&1
bash $W/sd16/cl/snp1.sh $c > $W/logs/sd16_$c.log 2>&1
echo "DONE $(date -Is) $(hostname)" > $W/logs/stageA_$c.done
