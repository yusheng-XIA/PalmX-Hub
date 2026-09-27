#!/bin/bash
# Genomic context and repeat association of the linear SV catalogue (BEDTools v2.31.1)
# One class per SV in the priority order CDS > intron > 2-kb upstream > 2-kb downstream > intergenic
# (gene models carry no UTR, so CDS and exons coincide). INS: +-50 bp around the breakpoint; other types: interval.
# TE association: >= 50% overlap with EDTA TEs (INS: within +-50 bp of the breakpoint).
set -euo pipefail
G=FL_Hap2

awk 'BEGIN{OFS="\t"} $3=="CDS"{print $1,$4-1,$5}' ${G}.gff3 | sort -k1,1 -k2,2n | bedtools merge > cds.bed
awk 'BEGIN{OFS="\t"} $3=="gene"{print $1,$4-1,$5,".",".",$7}' ${G}.gff3 | sort -k1,1 -k2,2n > genes.bed
bedtools subtract -a genes.bed -b cds.bed > intron.bed
bedtools flank -i genes.bed -g ${G}.genome -l 2000 -r 0 -s > up2k.bed
bedtools flank -i genes.bed -g ${G}.genome -l 0 -r 2000 -s > down2k.bed

# sv.bed: chrom start end id type (INS as breakpoint +-50 bp)
for f in cds intron up2k down2k; do
    bedtools intersect -u -a sv.bed -b ${f}.bed | cut -f4 > hit_${f}.txt
done
awk 'FNR==1{k++} k==1{c[$1]=1; next} k==2{i[$1]=1; next} k==3{u[$1]=1; next} k==4{d[$1]=1; next}
     {cls = ($4 in c) ? "CDS" : ($4 in i) ? "intron" : ($4 in u) ? "upstream_2kb" : ($4 in d) ? "downstream_2kb" : "intergenic"
      print $4 "\t" $5 "\t" cls}' hit_cds.txt hit_intron.txt hit_up2k.txt hit_down2k.txt sv.bed > sv_context.tsv

bedtools intersect -u -f 0.5 -a sv.intervals.bed -b ${G}.EDTA.TEanno.bed > sv_te_associated.bed
bedtools intersect -u -a sv.ins_window.bed -b ${G}.EDTA.TEanno.bed >> sv_te_associated.bed
# Context within 5 kb of large SyRI breakpoints: TE, BISER segmental duplications and TRF tandem repeats
bedtools slop -i syri_breakpoints.bed -g ${G}.genome -b 5000 \
    | bedtools annotate -i - -files ${G}.EDTA.TEanno.bed ${G}.biser.bed ${G}.trf.bed > breakpoint_context.tsv
