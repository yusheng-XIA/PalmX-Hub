#!/bin/bash
# DNA allele counts at FL diagnostic sites: short reads (MGI, bwa) and HiFi (minimap2), each mapped separately
# to Africa_hap2 (A) and American_hap1 (B). Region-restricted via index (samtools view -M -L), then mpileup.
set -uo pipefail
W=${CLUSTER_WORK}/enh_A
ST=samtools
PY=python
P=${ANALYSIS_DIR}/01_seedless_results/12_evaluate
mkdir -p $W/pileup
awk 'BEGIN{OFS="\t"}{print $1,$2-1,$2}' $W/sites/FL_A_positions.tsv > $W/sites/FL_A.bed
awk 'BEGIN{OFS="\t"}{print $1,$2-1,$2}' $W/sites/FL_B_positions.tsv > $W/sites/FL_B.bed
run() { # s(01_africa/02_american) b(bam prefix) side(A/B) extra-flags
  s=$1; b=$2; side=$3; ff=$4
  bam=$P/$s/02_remapping/$b.sorted.bam; bai=$W/bai/$s.$b.bai
  $ST view -@2 -u -M -L $W/sites/FL_$side.bed "$bam##idx##$bai" \
   | $ST mpileup -B -d 100000 -q 20 -Q 20 $ff --no-output-ins --no-output-ins --no-output-del --no-output-del --no-output-ends -l $W/sites/FL_${side}_positions.tsv - 2> $W/logs/mp.$side.$b.err \
   | $PY $W/parse_mp.py | gzip > $W/pileup/$side.$b.counts.tsv.gz
  echo "done $side $b $(date)" >> $W/logs/s04.done
}
run 01_africa seedless3-10 A "--ff UNMAP,SECONDARY,QCFAIL,DUP,SUPPLEMENTARY" &
run 01_africa hifi_remap A "--ff UNMAP,SECONDARY,QCFAIL,DUP,SUPPLEMENTARY" &
run 02_american seedless3-10 B "--ff UNMAP,SECONDARY,QCFAIL,DUP,SUPPLEMENTARY" &
run 02_american hifi_remap B "--ff UNMAP,SECONDARY,QCFAIL,DUP,SUPPLEMENTARY" &
wait; echo ALLDONE >> $W/logs/s04.done
