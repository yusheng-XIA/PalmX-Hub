#!/bin/bash
# Genome assembly using Verkko with HiFi + ONT + Pore-C
# Used for the seedless interspecific hybrid (E. guineensis × E. oleifera)
# Reference: Rautiainen et al. (2023) Nature Biotechnology

threads=64

# ============================================================
# Verkko assembly with HiFi + ONT+Pore-C
# ============================================================
verkko -d asm \
  --hifi hifi/*.fastq.gz \
  --nano ont/*.fastq.gz \
  --porec porec/*fastq.gz
# ============================================================
# Haplotype phasing using SubPhaser
# Assign contigs to E. guineensis (Hap B) and E. oleifera (Hap A)
# ============================================================
SubPhaser -i ${sample}_verkko_asm/assembly.fasta \
    -labels HapA HapB \
    -o ${sample}_subphaser_out \
    -t ${threads}

# ============================================================
# Pore-C scaffolding
# Process Pore-C data using pairtools pipeline
# ============================================================

# Index assembly
bwa-mem2 index ${sample}.assembly.fasta

# Align Pore-C reads
pairtools parse \
    --min-mapq 40 \
    --walks-policy 5unique \
    --max-inter-align-gap 30 \
    --output ${sample}.parsed.pairsam \
    <(bwa-mem2 mem -t ${threads} -SP ${sample}.assembly.fasta ${sample}.porec.fq.gz)

pairtools sort --output ${sample}.sorted.pairsam ${sample}.parsed.pairsam
pairtools dedup --output ${sample}.dedup.pairsam ${sample}.sorted.pairsam
pairtools select '(pair_type == "UU") or (pair_type == "UR") or (pair_type == "RU")' \
    --output ${sample}.filtered.pairsam ${sample}.dedup.pairsam

# ============================================================
# Gap closing with ONT reads
# ============================================================
# Identify reads spanning gap junctions
minimap2 -ax map-ont -t ${threads} ${sample}.scaffolds.fasta ${sample}.ont.fq.gz \
    | samtools sort -@ ${threads} -o ${sample}.ont_aligned.bam
samtools index ${sample}.ont_aligned.bam

# Manual validation of gap-filling reads in IGV

# ============================================================
# Telomere extension
# ============================================================
# Extract reads containing telomeric repeats
grep -B1 "TTTAGGGTTTAGGG" ${sample}.hifi.fq | grep "^@" | sed 's/@//' > telomere_reads.list
seqtk subseq ${sample}.hifi.fq.gz telomere_reads.list > telomere_reads.fq
