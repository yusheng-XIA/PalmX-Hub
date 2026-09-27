# Weighted gene co-expression network (WGCNA v1.72) on the replicated 114-library FL/TN series
# (the single-library TK and NS series are excluded from the confirmatory network)
library(WGCNA)
library(DESeq2)
options(stringsAsFactors = FALSE)
allowWGCNAThreads()

counts  <- as.matrix(read.delim("gene_counts_114.tsv", row.names = 1, check.names = FALSE))
coldata <- read.delim("samples_114.tsv", row.names = 1)            # material, stage, replicate
counts  <- counts[, rownames(coldata)]

# 1. Genes with counts >= 10 in >= 10% of samples; DESeq2 variance-stabilising transformation
keep <- rowSums(counts >= 10) >= ceiling(0.10 * ncol(counts))
dds  <- DESeqDataSetFromMatrix(counts[keep, ], coldata, design = ~ 1)
vst_mat <- assay(vst(dds, blind = TRUE))

# 2. Top 5,000 genes by the mean rank of variance and median absolute deviation
rk <- (rank(-apply(vst_mat, 1, var)) + rank(-apply(vst_mat, 1, mad))) / 2
datExpr <- t(vst_mat[order(rk)[1:5000], ])

# 3. Sample outliers: standardized connectivity z < -2.5 flagged
A <- adjacency(t(datExpr), type = "signed")
k <- colSums(A) - 1
z_k <- (k - mean(k)) / sd(k)
write.table(data.frame(sample = rownames(datExpr), z_connectivity = z_k, outlier = z_k < -2.5),
            "sample_connectivity.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

# 4. Soft threshold: first power reaching a scale-free fit of 0.749 -> 22
sft <- pickSoftThreshold(datExpr, powerVector = 1:30, networkType = "signed")
net <- blockwiseModules(datExpr, power = 22, networkType = "signed", TOMType = "signed",
                        minModuleSize = 30, mergeCutHeight = 0.25, reassignThreshold = 0,
                        pamRespectsDendro = FALSE, maxBlockSize = 6000, numericLabels = TRUE)
moduleColors <- labels2colors(net$colors)
write.table(data.frame(gene = colnames(datExpr), module = moduleColors), "gene_modules.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)

# 5. Module eigengenes vs material, stage, protein and metabolite abundance (Pearson, BH correction)
MEs    <- orderMEs(moduleEigengenes(datExpr, moduleColors)$eigengenes)
traits <- read.delim("traits_114.tsv", row.names = 1)[rownames(datExpr), ]   # numeric-coded material/stage,
                                                                              # protein and metabolite abundances
r <- cor(MEs, traits, use = "p")
p <- corPvalueStudent(r, nrow(datExpr))
padj <- matrix(p.adjust(p, method = "BH"), nrow(p), dimnames = dimnames(p))
write.table(cbind(as.data.frame(as.table(r)), P = as.vector(p), P_adj = as.vector(padj)),
            "module_trait_correlation.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
