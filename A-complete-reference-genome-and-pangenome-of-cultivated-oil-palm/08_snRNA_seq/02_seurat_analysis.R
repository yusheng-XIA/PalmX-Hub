# Single-nucleus RNA-seq analysis: FL and TN mesocarp at 95, 125 and 185 d
# R v4.2.3, Seurat v5.3.0, Harmony v1.2.3, clusterProfiler v4.20.0

library(Seurat)
library(harmony)
library(dplyr)

set.seed(10086)
samples <- c("FL_95d", "FL_125d", "FL_185d", "TN_95d", "TN_125d", "TN_185d")

# ------------------------------------------------------------
# 1. One Seurat object per sample: genes detected in >= 5 nuclei, nuclei with >= 300 genes.
#    No additional 200-5,000 gene-range or mitochondrial-proportion filter was applied.
# ------------------------------------------------------------
obj_list <- lapply(samples, function(s) {
    counts <- Read10X(data.dir = file.path(s, "outs/filtered_feature_bc_matrix"))
    obj <- CreateSeuratObject(counts = counts, project = s, min.cells = 5, min.features = 300)
    obj$variety <- sub("_.*", "", s)
    obj$timepoint <- sub(".*_", "", s)
    obj
})
sce <- merge(obj_list[[1]], obj_list[-1], add.cell.ids = samples)
sce <- JoinLayers(sce)

# ------------------------------------------------------------
# 2. Normalization, PCA and Harmony batch correction by sample (orig.ident)
# ------------------------------------------------------------
sce <- NormalizeData(sce, normalization.method = "LogNormalize", scale.factor = 10000)
sce <- FindVariableFeatures(sce)
sce <- ScaleData(sce)
sce <- RunPCA(sce)
sce <- RunHarmony(sce, group.by.vars = "orig.ident")

# ------------------------------------------------------------
# 3. Neighbour graph, UMAP and Louvain clustering on Harmony dimensions 1-15
# ------------------------------------------------------------
sce <- FindNeighbors(sce, reduction = "harmony", dims = 1:15)
sce <- RunUMAP(sce, reduction = "harmony", dims = 1:15, seed.use = 10086)
sce <- FindClusters(sce, algorithm = 1, resolution = 1.0, random.seed = 10086)   # 20 cell states
saveRDS(sce, "sce.all_int.rds")

# ------------------------------------------------------------
# 4. Cell-state markers and FL-versus-TN differences within each cell state (Wilcoxon)
# ------------------------------------------------------------
Idents(sce) <- "seurat_clusters"
markers <- FindAllMarkers(sce, only.pos = TRUE, test.use = "wilcox")
write.table(markers, "cluster_markers.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

de_all <- list()
for (cl in levels(sce)) {
    sub <- subset(sce, idents = cl)
    if (min(table(sub$variety)) < 3) next
    de <- FindMarkers(sub, group.by = "variety", ident.1 = "FL", ident.2 = "TN", test.use = "wilcox")
    de$gene <- rownames(de); de$cluster <- cl
    de_all[[cl]] <- de
}
de_all <- bind_rows(de_all)
de_sig <- subset(de_all, p_val_adj < 0.05 & abs(avg_log2FC) > 0.25)
write.table(de_sig, "FL_vs_TN_within_cluster.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

# ------------------------------------------------------------
# 5. Cell-state composition by material and developmental stage
# ------------------------------------------------------------
comp <- as.data.frame(table(cluster = sce$seurat_clusters, variety = sce$variety, timepoint = sce$timepoint))
comp <- comp %>% group_by(variety, timepoint) %>% mutate(fraction = Freq / sum(Freq))
write.table(comp, "cluster_composition.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

# GO enrichment of cluster markers: see 07_multi_omics/09_go_enrichment.R (background = tested genes with GO terms)
