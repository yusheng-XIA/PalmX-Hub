#!/bin/bash
# Filter chain with every site-level filter computed among the 308 accessions only:
# samples restricted to the 308 right after normalisation; hard filters; genotype masking GQ<20|DP<5|DP>100;
# FORMAT reduced to GT after masking; then F_MISSING<=0.25 and AC>=1 (AC<AN) among the 308 -> div308; MAF>=0.01 among the 308 -> pop308.
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 GOTO_NUM_THREADS=1
set -euo pipefail
E=${CLUSTER_HOME}/miniconda3/envs/genomics_a2/bin
REF=${GWAS_DIR}/02_analysis/01_ref/genome.fasta
S308=${CLUSTER_WORK}/fig4c_final/meta/s308.txt
RAW=$1; REG=$2; OUT=$3
mkdir -p $(dirname $OUT)
$E/bcftools norm -r $REG --regions-overlap pos -f $REF -m -any -Ou $RAW \
 | $E/bcftools view -S $S308 -v snps -m2 -M2 -Ou \
 | $E/bcftools annotate --set-id %CHROM:%POS:%REF:%FIRST_ALT -Ou \
 | $E/bcftools filter -s SNP_HARD_FILTER -e "INFO/QD<2.0 || INFO/MQ<40.0 || INFO/FS>60.0 || INFO/SOR>3.0 || INFO/MQRankSum<-12.5 || INFO/ReadPosRankSum<-8.0" -Ou \
 | $E/bcftools view -f PASS,. -Ou \
 | $E/bcftools filter -S . -e "FMT/GQ<20 | FMT/DP<5 | FMT/DP>100" -Ou \
 | $E/bcftools annotate -x ^FORMAT/GT -Ou \
 | $E/bcftools +fill-tags -Ou -- -t AC,AN,AF,MAF,F_MISSING,NS \
 | $E/bcftools view -i "INFO/F_MISSING<=0.25 && INFO/AC>=1 && INFO/AC<INFO/AN" -Ob -o $OUT.div308.bcf
$E/bcftools index $OUT.div308.bcf
$E/bcftools view -i "INFO/MAF>=0.01" -Ob -o $OUT.pop308.bcf $OUT.div308.bcf
$E/bcftools index $OUT.pop308.bcf
echo "DONE $(date -Is) $(hostname)" > $OUT.done
