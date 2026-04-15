#!/bin/bash
# Allelic gene identification between haplotypes
threads=64

# ============================================================
# 1. Synteny-based allelic gene identification (JCVI)
# ============================================================
python -m jcvi.formats.gff bed --key=ID ${sample}.hap1.gff3 -o hap1.bed
python -m jcvi.formats.gff bed --key=ID ${sample}.hap2.gff3 -o hap2.bed

python -m jcvi.compara.catalog ortholog \
    --cscore=0.99 hap1 hap2

# ============================================================
# 2. Coordinate-based allelic assignment (BLASTN + GMAP)
# ============================================================
# BLASTN
makeblastdb -in hap1.cds.fna -dbtype nucl
blastn -db hap1.cds.fna -query hap2.cds.fna \
    -evalue 1e-5 -outfmt 6 \
    -out hap2_vs_hap1.blastn \
    -num_threads ${threads}

# GMAP
gmap_build -d hap1_db ${sample}.hap1.fasta
gmap -t ${threads} -d hap1_db ${sample}.hap2.cds.fna -f samse \
    > hap2_to_hap1.gmap.sam

# ============================================================
# 3. Classify allelic pairs
# ============================================================
# Retain pairs with >50% coordinate overlap, >80% identity, >80% alignment length
# Same-CDS alleles: identical CDS sequences
# Haplotype-specific: unpaired genes

# ============================================================
# 4. Allele-specific expression (ASE)
# ============================================================
# Align RNA-seq to allele-aware gene models
STAR --runThreadN ${threads} \
    --genomeDir ${sample}_star_index \
    --readFilesIn ${tissue}_R1.fq.gz ${tissue}_R2.fq.gz \
    --readFilesCommand zcat \
    --alignIntronMax 20000 \
    --alignMatesGapMax 25000 \
    --outFilterMultimapNmax 1 \
    --outSAMtype BAM SortedByCoordinate

# Quantify expression
stringtie -e -B -p ${threads} \
    -G ${sample}.allele_aware.gtf \
    -o ${tissue}.gtf \
    ${tissue}.sorted.bam

# Alternative: Salmon quantification
salmon index -t ${sample}.allele_transcript.fa \
    -i ${sample}_salmon_index \
    --decoys ${sample}.decoys.txt -k 31

salmon quant -i ${sample}_salmon_index -l A \
    -1 ${tissue}_R1.fq.gz -2 ${tissue}_R2.fq.gz \
    --validateMappings --gcBias \
    -o ${tissue}_quant -p ${threads}
