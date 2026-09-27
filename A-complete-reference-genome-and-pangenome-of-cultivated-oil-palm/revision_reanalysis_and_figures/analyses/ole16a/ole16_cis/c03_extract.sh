#!/bin/bash
# per-genome: map FL-Hap2 OLE16a neighbourhood to chr11, extract locus (+/-15 kb)
set -uo pipefail
g=$1
D=${ANALYSIS_DIR}/08_hifi_chromosome/11_final_corrected_39genomes_annotations_20260812/attempt_20260812_01/01_genomes
W=${CLUSTER_WORK}/ole16_cis; cd $W; export TMPDIR=$W/tmp
MM=minimap2
ST=samtools
mkdir -p loc
$ST faidx $D/$g.fa chr11 > tmp/$g.chr11.fa
$MM -c -x asm20 -t 2 tmp/$g.chr11.fa q/query45k.fa 2>/dev/null > loc/$g.q45.paf
# best primary chain
read ts te strand <<<$(awk '$12>=0 {print $8,$9,$5,$10}' loc/$g.q45.paf | sort -k4,4nr | head -1 | awk '{print $1,$2,$3}')
if [ -z "${ts:-}" ]; then echo "$g NOHIT" > loc/$g.status; rm -f tmp/$g.chr11.fa*; exit; fi
a=$((ts-15000)); [ $a -lt 1 ] && a=1; b=$((te+15000))
$ST faidx tmp/$g.chr11.fa chr11:$((a+1))-$b | sed "s/^>.*/>${g}|chr11:$((a+1))-$b|$strand/" > loc/$g.locus.fa
echo "$g chr11 $((a+1)) $b $strand" > loc/$g.status
rm -f tmp/$g.chr11.fa tmp/$g.chr11.fa.fai
