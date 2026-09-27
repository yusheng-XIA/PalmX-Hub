#!/bin/bash
# Exact re-implementation of the repair-run filter chain recorded in the headers of
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 GOTO_NUM_THREADS=1
# snp_repair/06_filter/*.corrected.vcf.gz (bcftools 1.21 there; 1.22 here).
set -euo pipefail
E=${CLUSTER_HOME}/miniconda3/envs/genomics_a2/bin
REF=${GWAS_DIR}/02_analysis/01_ref/genome.fasta
RAW=$1; CHR=$2; OUT=$3; REG=${4:-$CHR}; T=${5:-4}
mkdir -p $(dirname $OUT)
$E/bcftools norm --threads $T -r $REG --regions-overlap pos -f $REF -m -any -Ou $RAW \
 | $E/bcftools view --threads $T -v snps -m2 -M2 -Ou \
 | $E/bcftools annotate --threads $T --set-id %CHROM:%POS:%REF:%FIRST_ALT -Ou \
 | $E/bcftools filter --threads $T -s SNP_HARD_FILTER -e "INFO/QD<2.0 || INFO/MQ<40.0 || INFO/FS>60.0 || INFO/SOR>3.0 || INFO/MQRankSum<-12.5 || INFO/ReadPosRankSum<-8.0" -Ou \
 | $E/bcftools view --threads $T -f PASS,. -Ou \
 | $E/bcftools filter --threads $T -S . -e "FMT/GQ<20 | FMT/DP<5 | FMT/DP>100" -Ou \
 | $E/bcftools +fill-tags -Ou -- -t AC,AN,AF,MAF,F_MISSING,NS \
 | $E/bcftools view --threads $T -i "INFO/F_MISSING<=0.25 && INFO/AC>=1 && INFO/AC<INFO/AN" -Oz -o $OUT.diversity.vcf.gz
$E/bcftools index -t --threads $T $OUT.diversity.vcf.gz
$E/bcftools view --threads $T -i "INFO/MAF>=0.01" -Oz -o $OUT.population.vcf.gz $OUT.diversity.vcf.gz
$E/bcftools index -t --threads $T $OUT.population.vcf.gz
echo DONE $(date -Is) > $OUT.done
