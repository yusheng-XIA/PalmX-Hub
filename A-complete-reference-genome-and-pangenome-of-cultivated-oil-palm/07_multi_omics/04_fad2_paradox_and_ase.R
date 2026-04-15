# FAD2 paradox analysis and ASE classification

library(DESeq2)
library(ggplot2)

# ============================================================
# 1. FAD2 paradox score
# ============================================================
# paradox_score = log2(RNA_rank / protein_rank)
# RNA_rank: percentile ranking of mean expression across 19 stages
# protein_rank: percentile ranking of mean protein abundance

rna_expr <- read.csv("rna_mean_expression.csv", row.names = 1)
protein_abund <- read.csv("protein_mean_abundance.csv", row.names = 1)

for (variety in c("FL", "BK", "TK", "NS")) {
    rna_rank <- rank(rna_expr[, variety]) / nrow(rna_expr)
    prot_rank <- rank(protein_abund[, variety]) / nrow(protein_abund)

    common_genes <- intersect(rownames(rna_expr), rownames(protein_abund))
    paradox_score <- log2(rna_rank[common_genes] / prot_rank[common_genes])

    # FL-specific translationally suppressed: paradox > 1.0 only in FL
    assign(paste0("paradox_", variety), paradox_score)
}

fl_suppressed <- names(paradox_FL[paradox_FL > 1.0 &
                                    paradox_BK <= 1.0 &
                                    paradox_TK <= 1.0 &
                                    paradox_NS <= 1.0])

# ============================================================
# 2. ASE analysis (DESeq2)
# ============================================================
# Allele-specific read counts for each gene pair
ase_counts <- read.csv("allele_counts.csv", row.names = 1)

dds_ase <- DESeqDataSetFromMatrix(
    countData = ase_counts,
    colData = ase_coldata,
    design = ~ haplotype)

dds_ase <- DESeq(dds_ase)
res_ase <- results(dds_ase, alpha = 0.05,
                    lfcThreshold = log2(2))  # FC > 2

# ============================================================
# 3. ASE classification
# ============================================================
# HapDom: one allele dominant in >=1/3 of stages
# Sub: each allele dominant in >=1 stage with switching
# NoDiff: remaining ASE genes

classify_ase <- function(ase_results, stages) {
    hap1_dominant <- sum(ase_results$log2FC > 1 & ase_results$padj < 0.05)
    hap2_dominant <- sum(ase_results$log2FC < -1 & ase_results$padj < 0.05)
    total <- length(stages)

    if (hap1_dominant >= total / 3 & hap2_dominant == 0) return("HapDom_H1")
    if (hap2_dominant >= total / 3 & hap1_dominant == 0) return("HapDom_H2")
    if (hap1_dominant >= 1 & hap2_dominant >= 1) return("Sub")
    return("NoDiff")
}

# ============================================================
# 4. Four-phase temporal ASE
# ============================================================
# Early: 0-50 d, Mid: 65-110 d, Late: 125-185 d, Postharvest: 12-72 h
phases <- list(
    Early = c("0d", "10d", "20d", "30d", "40d", "50d"),
    Mid = c("65d", "80d", "95d", "110d"),
    Late = c("125d", "140d", "185d"),
    Postharvest = c("12h", "24h", "36h", "48h", "60h", "72h")
)

# ============================================================
# 5. Ka/Ks analysis
# ============================================================
# Pairwise CDS alignment: ParaAT + KaKs_Calculator (NG method)
# system("ParaAT.pl -h allele_pairs.txt -n cds.fna -a protein.faa -p ${threads} -o kaks_aln")
# system("KaKs_Calculator -i kaks_aln -o kaks_results.txt -m NG")
# Filter: Ks > 5 or Ka/Ks > 5 excluded

# ============================================================
# 6. Cross-variety ASE comparison
# ============================================================
# OrthoFinder orthologous pairs between FL and BK
# Consistent: same dominant haplotype
# Inconsistent: switched dominant haplotype
