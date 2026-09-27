#!/bin/bash
# SVs between the paired haplotypes of each haplotype-resolved material
# minimap2 v2.30-r1287 (asm5, --eqx), SyRI v1.7.1, plotsr v1.1.1
set -euo pipefail
threads=32
for m in FL TN TK NS Nigerian Eoleifera; do
    minimap2 -ax asm5 --eqx -t ${threads} ${m}_Hap1.fa ${m}_Hap2.fa > ${m}.hap1_hap2.sam
    syri -c ${m}.hap1_hap2.sam -r ${m}_Hap1.fa -q ${m}_Hap2.fa -k -F S --nc ${threads} --prefix ${m}.
    printf "%s\t%s\n%s\t%s\n" ${m}_Hap1.fa ${m}_Hap1 ${m}_Hap2.fa ${m}_Hap2 > ${m}.genomes.txt
    plotsr --sr ${m}.syri.out --genomes ${m}.genomes.txt -o ${m}.plotsr.pdf
done
# TE overlap with EDTA annotations and genic/intergenic classification with BEDTools as in 07_sv_annotation.sh
