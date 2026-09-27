#!/bin/bash
# Transcript evidence
#   Public RNA-seq: leaf SRR6835890, root SRR6835891, embryo SRR11912814, female flower SRR1583660,
#                   male flower SRR1583661, stem core SRR25120000; in-house young leaf and fruit RNA-seq
#   Full-length: PacBio Iso-Seq (7.04 Gb) and ONT (14.37 Gb) from leaves
# HISAT2 v2.2.1, StringTie v2.2.1, TransDecoder v5.5.0
set -euo pipefail
threads=64

hisat2-build -p ${threads} ${sample}.fa ${sample}.hisat2
for lib in $(cat rnaseq_libraries.txt); do
    hisat2 -p ${threads} --dta -x ${sample}.hisat2 -1 ${lib}_R1.fq.gz -2 ${lib}_R2.fq.gz \
        | samtools sort -@ ${threads} -o ${lib}.bam
    stringtie -p ${threads} -o ${lib}.gtf ${lib}.bam
done
samtools merge -@ ${threads} ${sample}.rnaseq.merged.bam $(sed 's/$/.bam/' rnaseq_libraries.txt)
stringtie --merge -p ${threads} -o ${sample}.stringtie.merged.gtf $(sed 's/$/.gtf/' rnaseq_libraries.txt)

# Long reads (Iso-Seq and ONT) aligned with minimap2 splice mode and collapsed with StringTie long-read mode
minimap2 -ax splice:hq -uf -t ${threads} ${sample}.fa isoseq.fq.gz | samtools sort -o ${sample}.isoseq.bam
minimap2 -ax splice -uf -k14 -t ${threads} ${sample}.fa ont_rna.fq.gz | samtools sort -o ${sample}.ont.bam
stringtie -L -p ${threads} -o ${sample}.isoseq.gtf ${sample}.isoseq.bam
stringtie -L -p ${threads} -o ${sample}.ont.gtf ${sample}.ont.bam
stringtie --merge -p ${threads} -o ${sample}.transcripts.gtf \
    ${sample}.stringtie.merged.gtf ${sample}.isoseq.gtf ${sample}.ont.gtf

# Candidate coding regions (TransDecoder v5.5.0)
gtf_genome_to_cdna_fasta.pl ${sample}.transcripts.gtf ${sample}.fa > transcripts.fa
gtf_to_alignment_gff3.pl ${sample}.transcripts.gtf > transcripts.gff3
TransDecoder.LongOrfs -t transcripts.fa
TransDecoder.Predict -t transcripts.fa
cdna_alignment_orf_to_genome_orf.pl transcripts.fa.transdecoder.gff3 transcripts.gff3 transcripts.fa \
    > ${sample}.transdecoder.genome.gff3
