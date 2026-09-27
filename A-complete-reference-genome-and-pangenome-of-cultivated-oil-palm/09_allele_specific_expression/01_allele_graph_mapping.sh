#!/bin/bash
# Allele-specific expression (ASE): splice-aware allele graphs and RNA-seq graph mapping
# Run separately for the two hybrids:
#   FL : backbone = FL-Hap2 (E. guineensis), alternative path = FL-Hap1 (E. oleifera)
#   TN : backbone = dura (TK-like) haplotype, alternative path = pisifera (NS-like) haplotype
# Tools: Minigraph v0.21-r606, VG v1.67.0, bcftools, samtools
# Contig names of both haplotypes are prefixed "<assembly>__" before graph construction.
set -euo pipefail

threads=16
backbone=${backbone}          # e.g. FL_Hap2
alternative=${alternative}    # e.g. FL_Hap1
prefix=${prefix}              # output prefix, e.g. FL
gff_backbone=${gff_backbone}  # backbone gene models (prefixed contig names; exon features with Parent=transcript)

# ============================================================
# 1. Diploid rGFA graph (backbone first)
# ============================================================
minigraph -cxggs -t ${threads} ${backbone}.prefixed.fa ${alternative}.prefixed.fa > ${prefix}.minigraph.gfa

# ============================================================
# 2. Biallelic diagnostic SNPs relative to the backbone
# ============================================================
vg deconstruct -P "${backbone}__" -C -a -t ${threads} ${prefix}.minigraph.gfa > ${prefix}.deconstruct.raw.vcf
bcftools norm -f ${backbone}.prefixed.fa -m -any ${prefix}.deconstruct.raw.vcf -Ou \
    | bcftools view -m2 -M2 -v snps -Oz -o ${prefix}.diagnostic_snps.vcf.gz
tabix -p vcf ${prefix}.diagnostic_snps.vcf.gz

# ============================================================
# 3. Splice-aware graph index (backbone + diagnostic SNPs + exon annotation)
# ============================================================
vg autoindex --workflow mpmap \
    --ref-fasta ${backbone}.prefixed.fa \
    --vcf ${prefix}.diagnostic_snps.vcf.gz \
    --tx-gff ${gff_backbone} --gff-feature exon --gff-tx-tag Parent \
    --prefix ${prefix}.rna_index --threads ${threads} --target-mem 120G

# Backbone path names used for surjection
awk -v p="${backbone}__" '{print p $1}' ${backbone}.fa.fai > ${prefix}.backbone_paths.txt

# ============================================================
# 4. Map each RNA-seq library (114 libraries: 19 stages x 3 replicates x FL/TN)
#    and project alignments onto backbone coordinates keeping pair information (-i)
# ============================================================
while read sample r1 r2; do
    vg mpmap -x ${prefix}.rna_index.xg -g ${prefix}.rna_index.gcsa -d ${prefix}.rna_index.dist \
        -f ${r1} -f ${r2} -n RNA -t ${threads} -N ${sample} -R ${sample} > ${sample}.gam

    vg surject -x ${prefix}.rna_index.xg -b -i -m -S -F ${prefix}.backbone_paths.txt \
        -N ${sample} -R ${sample} -t ${threads} ${sample}.gam > ${sample}.unsorted.bam

    samtools sort -@ ${threads} -o ${sample}.bam ${sample}.unsorted.bam
    samtools index ${sample}.bam
    samtools flagstat ${sample}.bam > ${sample}.flagstat.txt
done < ${prefix}_rnaseq_samples.tsv

# ============================================================
# 5. Unique-gene exonic diagnostic sites and fragment-level allele counts
# ============================================================
python 02_prepare_informative_sites.py --analysis ${prefix} \
    --vcf ${prefix}.diagnostic_snps.vcf.gz --gff ${gff_backbone} \
    --output ${prefix}.informative_sites.tsv.gz --summary ${prefix}.informative_sites.json

while read sample r1 r2; do
    # read pairs are grouped by read name with samtools collate before counting
    samtools collate -@ ${threads} -O -u ${sample}.bam \
    | python 03_count_fragments.py --analysis ${prefix} --sample ${sample} \
        --sites ${prefix}.informative_sites.tsv.gz --bam - \
        --mapq 20 --baseq 20 \
        --output counts/${sample}.gene_counts.tsv --summary counts/${sample}.qc.json
done < ${prefix}_rnaseq_samples.tsv

# Concatenate per-library counts for 04_ase_test.py
awk 'FNR>1 || NR==1' counts/*.gene_counts.tsv > ${prefix}.gene_counts.all.tsv
