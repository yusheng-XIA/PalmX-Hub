#!/bin/bash
# Transcriptome-based gene prediction
threads=64

# ============================================================
# 1. Short-read RNA-seq alignment (HISAT2 + StringTie)
# ============================================================
hisat2-build -p ${threads} ${sample}.nucleus.masked.fasta ${sample}_index

for tissue in leaf root embryo male_flower female_flower stem fruit; do
    hisat2 -p ${threads} --dta --max-intronlen 20000 \
        -x ${sample}_index \
        -1 ${tissue}_R1.fq.gz \
        -2 ${tissue}_R2.fq.gz \
        -S ${sample}.${tissue}.sam

    samtools sort -@ ${threads} -o ${sample}.${tissue}.bam ${sample}.${tissue}.sam
    stringtie -p ${threads} -o ${sample}.${tissue}.gtf ${sample}.${tissue}.bam
done

# Merge all BAM and GTF files
samtools merge -@ ${threads} ${sample}.merged.bam ${sample}.*.bam
samtools sort -@ ${threads} -o ${sample}.merged.sorted.bam ${sample}.merged.bam

ls ${sample}.*.gtf > mergelist.txt
stringtie --merge -p ${threads} -o ${sample}.merged.gtf mergelist.txt

# ============================================================
# 2. PacBio Iso-Seq processing (FLAIR)
# ============================================================
minimap2 -ax splice -t ${threads} \
    ${sample}.nucleus.fasta ${sample}.isoseq.fq.gz \
    | samtools sort -o ${sample}.isoseq.bam

flair align -g ${sample}.nucleus.fasta \
    -r ${sample}.isoseq.fq.gz \
    -t ${threads} \
    -o ${sample}.isoseq

flair correct -g ${sample}.nucleus.fasta \
    -q ${sample}.isoseq.bed \
    -f ${sample}.merged.gtf \
    -t ${threads} \
    -o ${sample}.isoseq.corrected

flair collapse -g ${sample}.nucleus.fasta \
    -r ${sample}.isoseq.fq.gz \
    -q ${sample}.isoseq.corrected.bed \
    -t ${threads} \
    -o ${sample}.isoseq.collapsed

# ============================================================
# 3. ONT direct RNA-seq processing
# ============================================================
minimap2 -ax splice -uf -k14 -t ${threads} \
    ${sample}.nucleus.fasta ${sample}.directRNA.fq.gz \
    | samtools sort -o ${sample}.directRNA.bam

# ============================================================
# 4. Trinity de novo and genome-guided assembly
# ============================================================
# De novo assembly
Trinity --normalize_reads --seqType fq --SS_lib_type FR \
    --left ${sample}_RNA_R1.fq.gz \
    --right ${sample}_RNA_R2.fq.gz \
    --CPU ${threads} --max_memory 500G \
    --output trinity_denovo

# Genome-guided assembly
Trinity --normalize_reads --seqType fq --SS_lib_type FR \
    --left ${sample}_RNA_R1.fq.gz \
    --right ${sample}_RNA_R2.fq.gz \
    --genome_guided_bam ${sample}.merged.sorted.bam \
    --genome_guided_max_intron 25000 \
    --CPU ${threads} --max_memory 500G \
    --output trinity_guided

cat trinity_denovo/Trinity.fasta trinity_guided/Trinity-GG.fasta \
    > all_transcripts.fasta

# ============================================================
# 5. PASA transcript assembly
# ============================================================
Launch_PASA_pipeline.pl \
    -c alignAssembly.config \
    --ALIGNERS gmap \
    -I 25000 \
    --cufflinks_gtf ${sample}.merged.gtf \
    -C -R \
    -g ${sample}.nucleus.masked.fasta \
    -t all_transcripts.fasta \
    --transcribed_is_aligned_orient \
    --stringent_alignment_overlap 30 \
    --CPU ${threads}

# ============================================================
# 6. TransDecoder for CDS prediction
# ============================================================
TransDecoder.LongOrfs -t all_transcripts.fasta
TransDecoder.Predict -t all_transcripts.fasta
