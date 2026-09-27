#!/bin/bash
# Assembly-based SNV calls of one donor haplotype against FL-Hap2 (minimap2 asm5 + paftools call),
# as used for the donor-path SNVs (e.g. the two haplotypes of the final Nigerian assembly).
# usage: 02_call_donor_snvs.sh SAMPLE_LABEL QUERY_FASTA
set -euo pipefail
export LC_ALL=C
SAMPLE="$1"; QRY_FASTA="$2"
W=${WORK:-donor_snv_calls}
REF_FASTA=${REF_FASTA:-FL_Hap2.fasta}
MINIMAP2=minimap2   # v2.28-r1209
K8=k8
PAFTOOLS=$(dirname $(which minimap2))/paftools.js
THREADS=16
OUT=$W/calls/$SAMPLE; mkdir -p $OUT $W/tmp
export TMPDIR=$W/tmp
[[ "$($MINIMAP2 --version)" == "2.28-r1209" ]]
[[ "$(awk 'END{print NR}' $QRY_FASTA.fai)" -eq 16 ]]
P=$OUT/${SAMPLE}_vs_FL_Hap2
echo "start $(date) $(hostname) $SAMPLE $QRY_FASTA" > $OUT/run.log
/usr/bin/time -v -o $P.minimap2.time.txt $MINIMAP2 -cx asm5 -t $THREADS --cs $REF_FASTA $QRY_FASTA > $P.paf 2> $P.minimap2.log
sort -T $W/tmp -S 8G -k6,6 -k8,8n $P.paf > $P.sorted.paf
$K8 $PAFTOOLS call -l 1000 -L 1000 $P.sorted.paf > $P.var.txt 2> $P.call.log
awk 'BEGIN{OFS="\t"} $1=="V" && $7!="" && $8!="" {ref=$7; alt=$8; if(ref=="-") ref="."; if(alt=="-") alt="."; if(ref!="." || alt!=".") print $2,$3,".",ref,alt,$6,"PASS","."}' $P.var.txt > $P.vcfbody
awk -F '\t' '$4 ~ /^[ACGTacgt]$/ && $5 ~ /^[ACGTacgt]$/' $P.vcfbody > $P.snp.txt
grep '^R' $P.var.txt > $P.R.txt || true
rm -f $P.vcfbody $P.var.txt $P.paf
echo "done $(date) snp=$(wc -l < $P.snp.txt)" >> $OUT/run.log
