#!/bin/bash
# Genome assembly using hifiasm with different data combinations
# Reference: Cheng et al. (2021) Nature Methods

threads=64

# ============================================================
# Mode 1: HiFi + ONT + Hi-C (for tenera hybrid Boke)
# ============================================================
hifiasm -o ${sample}.asm \
    -t ${threads} \
    --h1 ${sample}_HiC_R1.fq.gz \
    --h2 ${sample}_HiC_R2.fq.gz \
    --ul ${sample}.ont.fq.gz \
    ${sample}.hifi.fq.gz

# Extract primary and alternate assemblies
awk '/^S/{print ">"$2;print $3}' ${sample}.asm.hic.hap1.p_ctg.gfa > ${sample}.hap1.fa
awk '/^S/{print ">"$2;print $3}' ${sample}.asm.hic.hap2.p_ctg.gfa > ${sample}.hap2.fa

# ============================================================
# Mode 2: HiFi + ONT (for dura and pisifera)
# ============================================================
hifiasm -o ${sample}.asm \
    -t ${threads} \
    --ul ${sample}.ont.fq.gz \
    ${sample}.hifi.fq.gz

awk '/^S/{print ">"$2;print $3}' ${sample}.asm.bp.p_ctg.gfa > ${sample}.p_ctg.fa
awk '/^S/{print ">"$2;print $3}' ${sample}.asm.bp.a_ctg.gfa > ${sample}.a_ctg.fa

# ============================================================
# Mode 3: HiFi only (for 28 pangenome accessions)
# ============================================================
hifiasm -o ${sample}.asm \
    -t ${threads} \
    ${sample}.hifi.fq.gz

awk '/^S/{print ">"$2;print $3}' ${sample}.asm.bp.p_ctg.gfa > ${sample}.p_ctg.fa
