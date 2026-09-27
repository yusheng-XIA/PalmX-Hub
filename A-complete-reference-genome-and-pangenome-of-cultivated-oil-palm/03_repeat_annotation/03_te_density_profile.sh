#!/bin/bash
# TE density across a composite meta-gene: 2 kb upstream of the TSS (20 bins), gene body (100 bins),
# 2 kb downstream of the TES (20 bins); BEDTools v2.31.0.
# For display only, profiles were smoothed with Savitzky-Golay windows of 7, 11 and 7 bins.
set -euo pipefail

awk 'BEGIN{OFS="\t"} $3=="gene"{match($9,/ID=[^;]+/); id=substr($9,RSTART+3,RLENGTH-3); print $1,$4-1,$5,id,".",$7}' \
    ${sample}.gff3 > ${sample}.genes.bed
awk 'BEGIN{OFS="\t"} !/^#/ {print $1,$4-1,$5}' ${sample}.TEanno.gff3 | sort -k1,1 -k2,2n | bedtools merge > ${sample}.te.bed

# 140 bins per gene, oriented 5' -> 3'
awk 'BEGIN{OFS="\t"} {
    s=$2; e=$3; L=e-s; st=$6
    for (i=0; i<20; i++) {                                  # upstream 2 kb
        a = (st=="+") ? s-2000+i*100 : e+2000-(i+1)*100
        if (a >= 0) print $1, a, a+100, $4, i, st }
    for (i=0; i<100; i++) {                                 # gene body
        a = (st=="+") ? s+int(i*L/100) : e-int((i+1)*L/100)
        b = (st=="+") ? s+int((i+1)*L/100) : e-int(i*L/100)
        if (b > a) print $1, a, b, $4, 20+i, st }
    for (i=0; i<20; i++) {                                  # downstream 2 kb
        a = (st=="+") ? e+i*100 : s-(i+1)*100
        if (a >= 0) print $1, a, a+100, $4, 120+i, st }
}' ${sample}.genes.bed | sort -k1,1 -k2,2n > ${sample}.metagene_bins.bed

bedtools coverage -a ${sample}.metagene_bins.bed -b ${sample}.te.bed \
    | awk 'BEGIN{OFS="\t"} {cov[$5]+=$NF; n[$5]++} END{for (b in cov) print b, cov[b]/n[b]}' \
    | sort -k1,1n > ${sample}.te_density_profile.tsv
