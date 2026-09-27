#!/bin/bash
# FL (seedless interspecific hybrid Reyou-2): Verkko v2.2.1 with HiFi, ONT ultra-long and Pore-C
# Haplotypes assigned with SubPhaser v1.1 (subgenome-specific k-mers):
#   FL-Hap1 predominantly E. oleifera ancestry, FL-Hap2 predominantly E. guineensis ancestry
# Scaffolding with CPhasing v0.2.5.r291 (Pore-C), gap closing with TGS-GapCloser, telomere extension
set -euo pipefail
threads=64
sample=FL

# 1. Assembly
verkko -d ${sample}_verkko --hifi hifi/*.fastq.gz --nano ont/*.fastq.gz --porec porec/*.fastq.gz
cp ${sample}_verkko/assembly.fasta ${sample}.verkko.fa

# 2. Haplotype assignment by subgenome-specific k-mers (config lists contigs per homologous group)
SubPhaser.py -i ${sample}.verkko.fa -c subphaser_groups.cfg -pre ${sample} -t ${threads}

# 3. Pore-C scaffolding into 16 pseudochromosomes per haplotype
for hap in Hap1 Hap2; do
    cphasing pipeline -f ${sample}_${hap}.contigs.fa -pct porec/${sample}.porec.fq.gz -n 16 -t ${threads} \
        -o cphasing_${hap}
done

# 4. Gap closing with ONT reads
for hap in Hap1 Hap2; do
    tgsgapcloser --scaff ${sample}_${hap}.scaffolds.fa --reads ont/${sample}.ont.fa \
        --output ${sample}_${hap}.gapclosed --ne --thread ${threads}
done

# 5. Telomere extension: reads carrying (TTTAGGG)n / (CCCTAAA)n are extracted and assembled locally
#    at chromosome termini; every correction was checked on read pileups in IGV.
seqkit grep -s -r -p '(TTTAGGG){5,}|(CCCTAAA){5,}' hifi/*.fastq.gz ont/*.fastq.gz > telomeric_reads.fq
minimap2 -ax map-hifi -t ${threads} ${sample}_Hap2.gapclosed.scaff_seqs telomeric_reads.fq \
    | samtools sort -@ ${threads} -o telomeric_reads.bam
samtools index telomeric_reads.bam

# 6. Consistency checks: Pore-C contacts (CPhasing) and long-read support
for hap in Hap1 Hap2; do
    minimap2 -ax map-ont -t ${threads} ${sample}_${hap}.final.fa ont/${sample}.ont.fq.gz \
        | samtools sort -@ ${threads} -o ${sample}_${hap}.ont.bam
done
