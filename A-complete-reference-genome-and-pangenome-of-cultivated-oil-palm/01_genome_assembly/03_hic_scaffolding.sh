#!/bin/bash
# Hi-C scaffolding using HapHiC + YaHS + Juicebox
threads=64

# ============================================================
# Step 1: Hi-C read alignment using chromap
# ============================================================
samtools faidx ${sample}.fasta
chromap -i -r ${sample}.fasta -o contigs.index

chromap \
    --preset hic \
    -r ${sample}.fasta \
    -x contigs.index \
    --remove-pcr-duplicates \
    -1 ${sample}_HiC_R1.fq.gz \
    -2 ${sample}_HiC_R2.fq.gz \
    --SAM \
    -o ${sample}.sam \
    -t ${threads}

samtools view -bh -u -F0xF0C -q 10 ${sample}.sam \
    | samtools sort -@ ${threads} -o ${sample}.hic.bam
samtools index ${sample}.hic.bam

# Convert to BED for YaHS
bedtools bamtobed -i ${sample}.hic.bam \
    | awk -v OFS='\t' '{$4=substr($4,1,length($4)-2); print}' > ${sample}.bed

# ============================================================
# Step 2: HapHiC allele-aware scaffolding
# ============================================================
HapHiC pipeline ${sample}.fasta ${sample}.hic.bam ${nchrom}

# ============================================================
# Step 3: YaHS scaffolding
# ============================================================
yahs ${sample}.fasta ${sample}.bed -o ${sample}_yahs

# ============================================================
# Step 4: Juicebox visualization and manual curation
# ============================================================
juicer pre -a -o out_JBAT \
    ${sample}_yahs.bin \
    ${sample}_yahs_scaffolds_final.agp \
    ${sample}.fasta.fai

JUICER=juicer_tools.jar
asm_size=$(awk '{s+=$2} END{print s}' ${sample}.fasta.fai)
java -Xmx36G -jar ${JUICER} \
    pre out_JBAT.txt out_JBAT.hic assembly ${asm_size}

# After manual curation in Juicebox:
juicer post -o out_JBAT \
    out_JBAT.review.assembly \
    ${sample}_yahs_scaffolds_final.agp \
    ${sample}.fasta
