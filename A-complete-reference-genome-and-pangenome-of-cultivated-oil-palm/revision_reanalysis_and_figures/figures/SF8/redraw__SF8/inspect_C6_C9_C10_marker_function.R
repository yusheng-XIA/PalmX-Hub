suppressPackageStartupMessages({
  library(Matrix)
})

outdir <- "${ANALYSIS_DIR}/22_answer_reviews/00_ms/03_V3/03_figure3/03_single"
rds_file <- "${ANALYSIS_DIR}/20_results/10_database/14_scrna/raw/1-harmony/sce.all_int.rds"
fa_table <- "${ANALYSIS_DIR}/20_results/Figure3/03_omic/01_FA_identify/FA_core_enzyme_FINAL_detail.tsv"
annotation_file <- "${ANALYSIS_DIR}/20_results/10_database/04_annotations/Africa_hap2_annotations.tsv"

target_clusters <- c("C6", "C9", "C10")

obj <- readRDS(rds_file)
meta <- slot(obj, "meta.data")
rna <- slot(obj, "assays")[["RNA"]]
layers <- slot(rna, "layers")
expr <- if ("data" %in% names(layers)) layers[["data"]] else layers[[1]]
feature_names <- attributes(slot(rna, "features"))$dimnames[[1]]
cell_names <- attributes(slot(rna, "cells"))$dimnames[[1]]
rownames(expr) <- feature_names
colnames(expr) <- cell_names

common_cells <- intersect(rownames(meta), colnames(expr))
meta <- meta[common_cells, , drop = FALSE]
expr <- expr[, common_cells, drop = FALSE]
cluster <- paste0("C", as.character(meta[["seurat_clusters"]]))

anno <- read.delim(annotation_file, stringsAsFactors = FALSE, check.names = FALSE)
anno$gene_tu <- sub("^evm\\.model\\.", "evm.TU.", anno$gene_id)
keep_cols <- intersect(c("gene_tu", "gene_id", "product", "go_terms", "interpro_domains", "eggnog_ortholog", "kegg_ko", "description"), colnames(anno))
anno_small <- anno[, keep_cols, drop = FALSE]

marker_rows <- list()
for (cl in target_clusters) {
  in_idx <- cluster == cl
  out_idx <- cluster != cl
  mean_in <- Matrix::rowMeans(expr[, in_idx, drop = FALSE])
  mean_out <- Matrix::rowMeans(expr[, out_idx, drop = FALSE])
  pct_in <- Matrix::rowSums(expr[, in_idx, drop = FALSE] > 0) / sum(in_idx)
  pct_out <- Matrix::rowSums(expr[, out_idx, drop = FALSE] > 0) / sum(out_idx)
  stat <- data.frame(
    cluster = cl,
    gene_id = rownames(expr),
    mean_in = as.numeric(mean_in),
    mean_out = as.numeric(mean_out),
    log2fc_mean = log2((as.numeric(mean_in) + 1e-6) / (as.numeric(mean_out) + 1e-6)),
    pct_in = as.numeric(pct_in),
    pct_out = as.numeric(pct_out),
    pct_diff = as.numeric(pct_in - pct_out),
    stringsAsFactors = FALSE
  )
  stat$marker_rank_score <- stat$log2fc_mean * pmax(stat$pct_diff, 0)
  stat <- stat[is.finite(stat$marker_rank_score), , drop = FALSE]
  stat <- stat[order(-stat$marker_rank_score, -stat$log2fc_mean, -stat$pct_in), , drop = FALSE]
  marker_rows[[cl]] <- head(stat, 120)
}
markers <- do.call(rbind, marker_rows)
markers <- merge(markers, anno_small, by.x = "gene_id", by.y = "gene_tu", all.x = TRUE)
markers <- markers[order(match(markers$cluster, target_clusters), -markers$marker_rank_score), , drop = FALSE]
write.table(markers, file.path(outdir, "C6_C9_C10_top120_markers_with_annotation.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

score_genes <- function(expr, genes) {
  genes <- unique(genes[genes %in% rownames(expr)])
  if (length(genes) == 0) {
    return(list(score = rep(NA_real_, ncol(expr)), genes = genes))
  }
  m <- as.matrix(expr[genes, , drop = FALSE])
  gene_mean <- rowMeans(m)
  gene_sd <- apply(m, 1, sd)
  gene_sd[is.na(gene_sd) | gene_sd == 0] <- 1
  z <- sweep(sweep(m, 1, gene_mean, "-"), 1, gene_sd, "/")
  list(score = as.numeric(scale(colMeans(z))), genes = genes)
}

fa <- read.delim(fa_table, stringsAsFactors = FALSE, check.names = FALSE)
fa <- fa[fa$Genome == "Africa_hap2", , drop = FALSE]
genes_by_enzyme <- split(fa$GeneID, fa$Enzyme)
enz <- function(names) unique(unlist(genes_by_enzyme[names], use.names = FALSE))

lipid_accumulation_genes <- enz(c(
  "ACCase", "ACP", "FabD (MCAT)", "FabG (KAR)", "FabI (ENR)",
  "KAS I/II", "KASIII", "FATA/B", "SAD", "LACS", "KCS",
  "GPAT", "LPAT", "PAP", "DGAT", "PDAT", "PDCT"
))
oleic_genes <- enz(c("SAD", "FATA/B"))
tag_genes <- enz(c("GPAT", "LPAT", "PAP", "DGAT", "PDAT", "PDCT"))

anno_text <- paste(anno$go_terms, anno$interpro_domains, anno$product, sep = " ")
rancidity_idx <- grepl(
  "triacylglycerol lipase|triglyceride catabolic|neutral lipid catabolic|lipase activity|phospholipase|lipoxygenase|GDSL|Patatin|Lipase|carboxylesterase|chlorophyll catabolic|senescence|oxidoreductase",
  anno_text,
  ignore.case = TRUE
)
rancidity_genes <- sub("^evm\\.model\\.", "evm.TU.", anno$gene_id[rancidity_idx])
storage_idx <- grepl(
  "oleosin|caleosin|lipid droplet|oil body|seed oil|diacylglycerol O-acyltransferase|triglyceride biosynthetic|triacylglycerol biosynthetic|neutral lipid biosynthetic|acylglycerol biosynthetic",
  anno_text,
  ignore.case = TRUE
)
storage_genes <- unique(c(tag_genes, sub("^evm\\.model\\.", "evm.TU.", anno$gene_id[storage_idx])))

scores <- list(
  lipid_accumulation = score_genes(expr, lipid_accumulation_genes),
  oleic_acid_axis = score_genes(expr, oleic_genes),
  oil_storage_TAG = score_genes(expr, storage_genes),
  rancidity_senescence = score_genes(expr, rancidity_genes)
)

score_dat <- data.frame(
  cell = common_cells,
  cluster = cluster,
  variety = as.character(meta[["variety"]]),
  timepoint = sub("d$", "", as.character(meta[["timepoint"]])),
  stringsAsFactors = FALSE
)
for (nm in names(scores)) {
  score_dat[[nm]] <- scores[[nm]]$score
}

score_rows <- list()
for (nm in names(scores)) {
  for (cl in sort(unique(score_dat$cluster))) {
    vals <- score_dat[score_dat$cluster == cl, nm]
    score_rows[[paste(nm, cl, sep = "_")]] <- data.frame(
      score = nm,
      cluster = cl,
      mean_z = mean(vals, na.rm = TRUE),
      median_z = median(vals, na.rm = TRUE),
      n_cells = length(vals),
      stringsAsFactors = FALSE
    )
  }
}
score_summary <- do.call(rbind, score_rows)
score_summary <- score_summary[order(score_summary$score, -score_summary$mean_z), , drop = FALSE]
write.table(score_summary, file.path(outdir, "C6_C9_C10_and_all_20cluster_functional_score_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

target_score_summary <- score_summary[score_summary$cluster %in% target_clusters, , drop = FALSE]
write.table(target_score_summary, file.path(outdir, "C6_C9_C10_functional_score_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

cat("Wrote marker and functional score summaries to", outdir, "\n")
cat("Top marker products:\n")
for (cl in target_clusters) {
  cat("\n", cl, "\n", sep = "")
  x <- markers[markers$cluster == cl, , drop = FALSE]
  cols <- intersect(c("gene_id", "log2fc_mean", "pct_in", "pct_out", "product"), colnames(x))
  print(head(x[, cols, drop = FALSE], 12))
}
cat("\nTarget functional score summary:\n")
print(target_score_summary)
