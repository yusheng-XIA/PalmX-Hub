#!/bin/bash
# hifiasm (v0.21.0-r686) assemblies
#   TN (tenera Boke)          : HiFi + ONT ultra-long + Hi-C
#   TK (dura), NS (pisifera)  : HiFi + ONT + Hi-C
#   Nigerian accession, independent E. oleifera accession : HiFi + Hi-C
#   27 additional pangenome accessions : HiFi only, primary contigs (--primary)
set -euo pipefail
threads=64
gfa2fa() { awk '/^S/{print ">"$2; print $3}' "$1"; }

# ---- Hi-C phased assemblies (TN, TK, NS: with ONT reads) ------------------------------------
hifiasm -o ${sample}.asm -t ${threads} \
    --ul ${sample}.ont.fq.gz \
    --h1 ${sample}_HiC_R1.fq.gz --h2 ${sample}_HiC_R2.fq.gz \
    ${sample}.hifi.fq.gz
gfa2fa ${sample}.asm.hic.hap1.p_ctg.gfa > ${sample}.hap1.p_ctg.fa
gfa2fa ${sample}.asm.hic.hap2.p_ctg.gfa > ${sample}.hap2.p_ctg.fa

# ---- Hi-C phased assemblies without ONT (Nigerian, E. oleifera) -----------------------------
hifiasm -o ${sample}.asm -t ${threads} \
    --h1 ${sample}_HiC_R1.fq.gz --h2 ${sample}_HiC_R2.fq.gz \
    ${sample}.hifi.fq.gz

# ---- 27 HiFi-only accessions: primary / alternate contigs -----------------------------------
hifiasm -o ${sample}.asm -t ${threads} --primary ${sample}.hifi.fq.gz
gfa2fa ${sample}.asm.p_ctg.gfa > ${sample}.p_ctg.fa     # main assembly sequence
gfa2fa ${sample}.asm.a_ctg.gfa > ${sample}.a_ctg.fa     # alternate-allele contigs (allele pairing only)
