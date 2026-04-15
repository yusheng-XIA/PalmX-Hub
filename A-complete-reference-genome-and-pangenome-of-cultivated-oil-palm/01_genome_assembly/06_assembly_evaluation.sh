#!/bin/bash
# Genome assembly quality assessment
threads=64

# ============================================================
# 1. Contig/scaffold statistics (N50)
# ============================================================
assembly-stats ${sample}.fasta

# ============================================================
# 2. BUSCO completeness
# ============================================================
busco -m genome \
    -i ${sample}.fasta \
    -l embryophyta_odb10 \
    -o ${sample}_busco \
    -c ${threads} \
    --offline -f

# ============================================================
# 3. LTR Assembly Index (LAI)
# ============================================================
ltr_finder_parallel -seq ${sample}.fasta -harvest_out -threads ${threads}
gt suffixerator -db ${sample}.fasta -indexname ${sample} -tis -suf -lcp -des -ssp -sds -dna
gt ltrharvest -index ${sample} -maxlenltr 7000 > ${sample}.harvest.scn

cat ${sample}.harvest.scn ${sample}.finder.combine.scn > ${sample}_rawLTR.scn
LTR_retriever -genome ${sample}.fasta \
    -inharvest ${sample}_rawLTR.scn \
    -threads ${threads} \
    -maxlenltr 15000

# ============================================================
# 4. QV and k-mer completeness (Merqury)
# ============================================================
meryl k=18 count ${sample}.hifi.fq.gz output ${sample}.meryl
merqury.sh ${sample}.meryl ${sample}.fasta ${sample}_merqury

# ============================================================
# 5. Assembly accuracy (Inspector)
# ============================================================
inspector.py -c ${sample}.fasta \
    -r ${sample}.hifi.fq.gz \
    -o ${sample}_inspector \
    --datatype hifi \
    -t ${threads}

# ============================================================
# 6. QUAST
# ============================================================
quast.py ${sample}.fasta \
    -R ${reference} \
    -g gene:${reference}.gff3 \
    -l ${sample} \
    -o ${sample}_quast \
    -t ${threads}

# ============================================================
# 7. Read mapping rate
# ============================================================
minimap2 -ax map-hifi -t ${threads} ${sample}.fasta ${sample}.hifi.fq.gz \
    | samtools sort -@ ${threads} -o ${sample}.hifi.bam
samtools flagstat ${sample}.hifi.bam > ${sample}.hifi.flagstat

minimap2 -ax map-ont -t ${threads} ${sample}.fasta ${sample}.ont.fq.gz \
    | samtools sort -@ ${threads} -o ${sample}.ont.bam
samtools flagstat ${sample}.ont.bam > ${sample}.ont.flagstat

# ============================================================
# 8. Syntenic dot plot
# ============================================================
nucmer -p ${sample}_vs_ref -t ${threads} ${reference} ${sample}.fasta
delta-filter -i 85 -l 1000 -1 ${sample}_vs_ref.delta > ${sample}_vs_ref.filtered.delta
mummerplot -p ${sample}_vs_ref ${sample}_vs_ref.filtered.delta \
    --color --filter --fat --png
