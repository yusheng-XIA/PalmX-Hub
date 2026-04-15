#!/bin/bash
# Karyotype evolution and biosynthetic gene cluster analysis

# ============================================================
# 1. Ancestral karyotype reconstruction (WGDI)
# ============================================================
# Syntenic blocks between oil palm and other palm species
wgdi -icl ${sample}.conf  # Identify collinear blocks
wgdi -ks ${sample}.conf   # Ks calculation for syntenic pairs
wgdi -bi ${sample}.conf   # Block information
wgdi -c ${sample}.conf    # Chromosome painting

# ============================================================
# 2. Biosynthetic gene cluster prediction (plantiSMASH)
# ============================================================
for sample in $(cat reference_assemblies.txt); do
    plantiSMASH ${sample}.fasta \
        --genefinding-gff3 ${sample}.gff3 \
        --output-dir ${sample}_bgc
done

# Mid-parent heterosis: MPH% = (hybrid - mid_parent) / mid_parent * 100
# Classify: over-dominance, dominance, under-dominance, below low-parent

# ============================================================
# 3. GO enrichment (used throughout the study)
# ============================================================
# R script using clusterProfiler
# Rscript -e '
# library(clusterProfiler)
# ego <- enrichGO(gene = gene_list, OrgDb = go_db,
#                  keyType = "GID", ont = "BP",
#                  pAdjustMethod = "BH", pvalueCutoff = 0.05)
# '
