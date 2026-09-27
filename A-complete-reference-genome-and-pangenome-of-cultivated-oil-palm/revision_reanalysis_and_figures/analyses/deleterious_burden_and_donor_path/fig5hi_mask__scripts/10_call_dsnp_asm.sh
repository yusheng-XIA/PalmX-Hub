#!/bin/bash
# Re-run of results-8.9/scripts/20_call_dsnp_new8_array.slurm (same binaries, same commands) for one query.
# usage: 10_call_dsnp_asm.sh SAMPLE_LABEL QUERY_FASTA
set -euo pipefail
export LC_ALL=C
SAMPLE="$1"; QRY_FASTA="$2"
W=${CLUSTER_WORK}/fig5hi_mask
REF_FASTA=${ANALYSIS_DIR}/14_pan_genome/04_SNP_calling/00_renamed_genomes/Africa_hap2.fasta
MINIMAP2=${DATA_DIR2}/anaconda3/envs/HiTE/bin/minimap2
K8=${DATA_DIR2}/anaconda3/envs/HiTE/bin/k8
PAFTOOLS=${DATA_DIR2}/anaconda3/envs/HiTE/bin/paftools.js
THREADS=16
OUT=$W/calls/$SAMPLE; mkdir -p $OUT $W/tmp
export TMPDIR=$W/tmp
[[ "$($MINIMAP2 --version)" == "2.28-r1209" ]]
[[ "$(sha256sum $MINIMAP2 | cut -d' ' -f1)" == 8a561f2c29f5a250eeecb6915100d23194dd38a41bb230d0b6afe69f78d5adf4 ]]
[[ "$(sha256sum $PAFTOOLS | cut -d' ' -f1)" == 0539eb3e4217920abedca07dd255a7955825a26c3f3dbdcaa2e2631ed4c99c3c ]]
[[ "$(awk 'END{print NR}' $QRY_FASTA.fai)" -eq 16 ]]
P=$OUT/${SAMPLE}_vs_Africa_hap2
echo "start $(date) $(hostname) $SAMPLE $QRY_FASTA" > $OUT/run.log
/usr/bin/time -v -o $P.minimap2.time.txt $MINIMAP2 -cx asm5 -t $THREADS --cs $REF_FASTA $QRY_FASTA > $P.paf 2> $P.minimap2.log
sort -T $W/tmp -S 8G -k6,6 -k8,8n $P.paf > $P.sorted.paf
$K8 $PAFTOOLS call -l 1000 -L 1000 $P.sorted.paf > $P.var.txt 2> $P.call.log
awk 'BEGIN{OFS="\t"} $1=="V" && $7!="" && $8!="" {ref=$7; alt=$8; if(ref=="-") ref="."; if(alt=="-") alt="."; if(ref!="." || alt!=".") print $2,$3,".",ref,alt,$6,"PASS","."}' $P.var.txt > $P.vcfbody
awk -F '\t' '$4 ~ /^[ACGTacgt]$/ && $5 ~ /^[ACGTacgt]$/' $P.vcfbody > $P.snp.txt
grep '^R' $P.var.txt > $P.R.txt || true
rm -f $P.vcfbody $P.var.txt $P.paf
echo "done $(date) snp=$(wc -l < $P.snp.txt)" >> $OUT/run.log
