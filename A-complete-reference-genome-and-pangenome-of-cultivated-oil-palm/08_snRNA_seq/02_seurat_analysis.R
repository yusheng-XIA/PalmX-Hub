# Single-nucleus RNA-seq analysis using Seurat v5
# FL and TN mesocarp at 95d, 125d, 185d

library(Seurat)
library(ggplot2)
library(dplyr)

# ============================================================
# 1. Load and create Seurat objects
# ============================================================
samples <- c("FL_95d", "FL_125d", "FL_185d",
             "TN_95d", "TN_125d", "TN_185d")

obj_list <- lapply(samples, function(s) {
    data <- Read10X(data.dir = file.path(s, "outs/filtered_feature_bc_matrix"))
    obj <- CreateSeuratObject(counts = data, project = s, min.cells = 3, min.features = 200)
    obj$variety <- sub("_.*", "", s)
    obj$stage <- sub(".*_", "", s)
    obj
})

merged <- merge(obj_list[[1]], obj_list[-1])

# ============================================================
# 2. Quality filtering
# ============================================================
merged[["percent.mt"]] <- PercentageFeatureSet(merged, pattern = "^Mt-|^mt-")

merged <- subset(merged,
                  subset = nFeature_RNA > 200 &
                           nFeature_RNA < 5000 &
                           percent.mt < 5)

# ============================================================
# 3. Normalization and integration
# ============================================================
merged <- SCTransform(merged, vars.to.regress = "percent.mt",
                       variable.features.n = 3000, verbose = FALSE)

# ============================================================
# 4. Dimensionality reduction and clustering
# ============================================================
merged <- RunPCA(merged, npcs = 30, verbose = FALSE)
merged <- FindNeighbors(merged, dims = 1:30)
merged <- FindClusters(merged, resolution = 0.8, algorithm = 1)  # Louvain
merged <- RunUMAP(merged, dims = 1:30)

# ============================================================
# 5. Differential expression between FL and TN
# ============================================================
Idents(merged) <- "seurat_clusters"

for (cluster in levels(merged)) {
    sub <- subset(merged, idents = cluster)
    de_results <- FindMarkers(sub,
                               group.by = "variety",
                               ident.1 = "FL",
                               ident.2 = "TN",
                               test.use = "wilcox",
                               logfc.threshold = 0.25,
                               min.pct = 0.1)
    de_sig <- de_results[de_results$p_val_adj < 0.05, ]
    write.csv(de_sig, paste0("DE_FL_vs_TN_cluster", cluster, ".csv"))
}

# ============================================================
# 6. Cell-type proportion analysis
# ============================================================
prop_table <- table(merged$seurat_clusters, merged$variety, merged$stage)
prop_df <- as.data.frame(prop.table(prop_table, margin = c(2, 3)))

# Variety-enriched clusters: >70% contribution from one variety
cluster_variety <- table(merged$seurat_clusters, merged$variety)
cluster_prop <- prop.table(cluster_variety, margin = 1)
enriched <- which(apply(cluster_prop, 1, max) > 0.7)

# ============================================================
# 7. GO enrichment of cluster markers
# ============================================================
library(clusterProfiler)

all_markers <- FindAllMarkers(merged, only.pos = TRUE,
                               min.pct = 0.25, logfc.threshold = 0.25)

for (cluster in unique(all_markers$cluster)) {
    genes <- all_markers$gene[all_markers$cluster == cluster]
    ego <- enrichGO(gene = genes,
                     OrgDb = go_annotation,
                     keyType = "GID",
                     ont = "BP",
                     pAdjustMethod = "BH",
                     pvalueCutoff = 0.05)
    write.csv(ego, paste0("GO_cluster", cluster, ".csv"))
}

# ============================================================
# 8. Visualization
# ============================================================
p1 <- DimPlot(merged, group.by = "seurat_clusters", label = TRUE)
p2 <- DimPlot(merged, group.by = "variety")
p3 <- DimPlot(merged, split.by = "stage")

ggsave("umap_clusters.pdf", p1, width = 10, height = 8)
ggsave("umap_variety.pdf", p2, width = 10, height = 8)
ggsave("umap_stages.pdf", p3, width = 24, height = 8)
