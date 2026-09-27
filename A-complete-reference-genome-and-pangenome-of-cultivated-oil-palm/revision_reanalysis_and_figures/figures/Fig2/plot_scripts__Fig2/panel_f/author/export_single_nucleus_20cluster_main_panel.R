suppressPackageStartupMessages({
  library(Matrix)
})

outdir <- "${ANALYSIS_DIR}/22_answer_reviews/00_ms/03_V3/03_figure3/03_single"
rds_file <- "${ANALYSIS_DIR}/20_results/10_database/14_scrna/raw/1-harmony/sce.all_int.rds"
fa_table <- "${ANALYSIS_DIR}/20_results/Figure3/03_omic/01_FA_identify/FA_core_enzyme_FINAL_detail.tsv"

obj <- readRDS(rds_file)
meta <- slot(obj, "meta.data")
umap <- slot(slot(obj, "reductions")[["umap"]], "cell.embeddings")
rna <- slot(obj, "assays")[["RNA"]]
layers <- slot(rna, "layers")

expr_layer <- if ("counts" %in% names(layers)) layers[["counts"]] else layers[[1]]
feature_names <- attributes(slot(rna, "features"))$dimnames[[1]]
cell_names <- attributes(slot(rna, "cells"))$dimnames[[1]]
rownames(expr_layer) <- feature_names
colnames(expr_layer) <- cell_names

common_cells <- Reduce(intersect, list(rownames(meta), rownames(umap), colnames(expr_layer)))
meta <- meta[common_cells, , drop = FALSE]
umap <- umap[common_cells, , drop = FALSE]
expr_layer <- expr_layer[, common_cells, drop = FALSE]

if (!"seurat_clusters" %in% colnames(meta)) {
  stop("Expected 20-cluster column 'seurat_clusters' was not found.")
}

fa <- read.delim(fa_table, stringsAsFactors = FALSE, check.names = FALSE)
fad2_genes <- unique(fa$GeneID[fa$Genome == "Africa_hap2" & fa$Enzyme == "FAD2"])
fad2_genes <- fad2_genes[fad2_genes %in% rownames(expr_layer)]
if (length(fad2_genes) == 0) {
  stop("No Africa_hap2 FAD2 genes matched the single-nucleus expression matrix.")
}

fad2_detected <- Matrix::colSums(expr_layer[fad2_genes, , drop = FALSE] > 0) > 0

dat <- data.frame(
  cell = common_cells,
  umap_1 = umap[, 1],
  umap_2 = umap[, 2],
  variety = as.character(meta[["variety"]]),
  timepoint = as.character(meta[["timepoint"]]),
  cluster = paste0("C", as.character(meta[["seurat_clusters"]])),
  fad2_detected = as.integer(fad2_detected),
  stringsAsFactors = FALSE
)

dat$timepoint <- sub("d$", "", dat$timepoint)
dat$timepoint_label <- paste0(dat$timepoint, " d")

cluster_summary <- as.data.frame(table(dat$variety, dat$timepoint, dat$cluster), stringsAsFactors = FALSE)
colnames(cluster_summary) <- c("variety", "timepoint", "cluster", "n_cells")
cluster_summary <- cluster_summary[cluster_summary$n_cells > 0, , drop = FALSE]
sample_totals <- aggregate(n_cells ~ variety + timepoint, cluster_summary, sum)
cluster_summary <- merge(cluster_summary, sample_totals, by = c("variety", "timepoint"), suffixes = c("", "_sample"))
cluster_summary$sample_fraction <- cluster_summary$n_cells / cluster_summary$n_cells_sample

mature <- dat[dat$timepoint == "185", , drop = FALSE]
focus <- mature[mature$cluster %in% c("C6", "C9", "C10"), , drop = FALSE]
focus_composition <- as.data.frame(table(focus$cluster, focus$variety), stringsAsFactors = FALSE)
colnames(focus_composition) <- c("cluster", "variety", "n_cells")
cluster_totals <- aggregate(n_cells ~ cluster, focus_composition, sum)
focus_composition <- merge(focus_composition, cluster_totals, by = "cluster", suffixes = c("", "_cluster"))
focus_composition$cluster_composition_fraction <- focus_composition$n_cells / focus_composition$n_cells_cluster

fad2_summary <- aggregate(fad2_detected ~ variety + timepoint, dat, function(x) mean(x > 0))
colnames(fad2_summary)[3] <- "fad2_detection_fraction"
fad2_counts <- aggregate(fad2_detected ~ variety + timepoint, dat, length)
colnames(fad2_counts)[3] <- "n_nuclei"
fad2_summary <- merge(fad2_summary, fad2_counts, by = c("variety", "timepoint"))

gene_set <- data.frame(
  gene_set = "FAD2_Africa_hap2",
  gene_id = fad2_genes,
  source = fa_table,
  stringsAsFactors = FALSE
)

write.table(dat, file.path(outdir, "single_nucleus_20cluster_main_panel_cells.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(cluster_summary, file.path(outdir, "single_nucleus_20cluster_sample_composition.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(focus_composition, file.path(outdir, "single_nucleus_20cluster_focus_cluster_composition_185d.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(fad2_summary, file.path(outdir, "single_nucleus_20cluster_fad2_detection_by_sample.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(gene_set, file.path(outdir, "single_nucleus_20cluster_fad2_gene_set.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

cat("Wrote 20-cluster single-nucleus main-panel tables to", outdir, "\n")
cat("Cells:", nrow(dat), "\n")
cat("FAD2 genes matched:", paste(fad2_genes, collapse = ","), "\n")
print(focus_composition[order(focus_composition$cluster, focus_composition$variety), ])
print(fad2_summary[fad2_summary$timepoint == "185", ])
