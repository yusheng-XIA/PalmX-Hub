# WGCNA co-expression network analysis
# 76 samples (4 varieties x 19 stages), top 5000 genes

library(WGCNA)
library(DESeq2)

options(stringsAsFactors = FALSE)
allowWGCNAThreads()

# ============================================================
# 1. Data preparation
# ============================================================
# Read count matrix (genes x samples)
counts <- read.csv("gene_counts_matrix.csv", row.names = 1)
coldata <- read.csv("sample_metadata.csv", row.names = 1)

# Variance-stabilized transformation
dds <- DESeqDataSetFromMatrix(countData = counts,
                               colData = coldata,
                               design = ~ variety + stage)
vsd <- varianceStabilizingTransformation(dds)
expr_data <- assay(vsd)

# Select top 5000 genes by median absolute deviation
mad_values <- apply(expr_data, 1, mad)
top_genes <- names(sort(mad_values, decreasing = TRUE))[1:5000]
datExpr <- t(expr_data[top_genes, ])

# ============================================================
# 2. Network construction
# ============================================================
# Soft-thresholding power selection
powers <- c(1:20)
sft <- pickSoftThreshold(datExpr, powerVector = powers,
                          networkType = "signed", verbose = 5)
# Selected power = 12 (R^2 > 0.85)

# Signed network construction
net <- blockwiseModules(datExpr,
                         power = 12,
                         networkType = "signed",
                         TOMType = "signed",
                         minModuleSize = 30,
                         reassignThreshold = 0,
                         mergeCutHeight = 0.25,
                         numericLabels = TRUE,
                         pamRespectsDendro = FALSE,
                         verbose = 3)

moduleLabels <- net$colors
moduleColors <- labels2colors(moduleLabels)

# ============================================================
# 3. Module-trait correlation
# ============================================================
traits <- read.csv("trait_data.csv", row.names = 1)
# Includes: variety, stage, metabolite concentrations, protein abundances

MEs <- moduleEigengenes(datExpr, moduleColors)$eigengenes
moduleTraitCor <- cor(MEs, traits, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nrow(datExpr))

# BH correction
moduleTraitPadj <- p.adjust(moduleTraitPvalue, method = "BH")

# ============================================================
# 4. Module membership and gene significance
# ============================================================
geneModuleMembership <- cor(datExpr, MEs, use = "p")
geneTraitSignificance <- cor(datExpr, traits$oil_content, use = "p")

# Export hub genes per module
for (mod in unique(moduleColors)) {
    inModule <- moduleColors == mod
    hub_genes <- names(sort(abs(geneModuleMembership[inModule, paste0("ME", mod)]),
                            decreasing = TRUE))[1:30]
    write.csv(hub_genes, paste0("hub_genes_", mod, ".csv"))
}
