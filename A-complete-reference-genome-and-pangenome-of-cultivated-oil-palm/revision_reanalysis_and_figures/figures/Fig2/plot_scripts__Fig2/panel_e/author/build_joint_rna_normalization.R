#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(DESeq2))

full_args <- commandArgs(trailingOnly = FALSE)
script_arg <- full_args[grep("^--file=", full_args)]
script_path <- sub("^--file=", "", script_arg)
here <- normalizePath(dirname(script_path), mustWork = TRUE)
analysis <- "${ANALYSIS_DIR}"
run_dir <- file.path(
  analysis,
  "22_answer_reviews/00_ms/03_V3/03_figure3/00_minipan",
  "07_two_variety_rna_timecourse/03_work/RUN-RNAFT-001"
)
counts_path <- file.path(
  analysis,
  "22_answer_reviews/00_ms/03_V3/03_figure3/00_minipan",
  "03_rnaseq_mapping/count_matrices/pangraphrna_hisat2_graph.gene_counts.tsv"
)
design_path <- file.path(run_dir, "00_manifest/design_all_114.tsv")
gene_list_path <- file.path(analysis, "21_MS/03_result/01_omic/FA_gene_list.tsv")
out_dir <- file.path(here, "ST9_latest_reanalysis_audit")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

safe_write <- function(x, path, row.names = FALSE) {
  temp <- paste0(path, ".tmp.", Sys.getpid())
  write.table(
    x, temp, sep = "\t", quote = FALSE, row.names = row.names,
    col.names = if (row.names) NA else TRUE
  )
  if (!file.rename(temp, path)) stop("Atomic rename failed: ", path)
}

design <- read.delim(design_path, check.names = FALSE, stringsAsFactors = FALSE)
if (nrow(design) != 114L || anyDuplicated(design$sample)) {
  stop("The RNA design is not 114 unique samples")
}

counts_df <- read.delim(counts_path, check.names = FALSE, stringsAsFactors = FALSE)
if (anyDuplicated(counts_df$gene_id)) stop("Raw count matrix has duplicate gene IDs")
rownames(counts_df) <- counts_df$gene_id
counts_df$gene_id <- NULL
counts_all <- as.matrix(counts_df)
storage.mode(counts_all) <- "integer"
if (!setequal(colnames(counts_all), design$sample)) stop("Count/design sample mismatch")
counts_all <- counts_all[, design$sample, drop = FALSE]

keep_for_size_factors <- rowSums(counts_all >= 10L) >= 3L
dds <- DESeqDataSetFromMatrix(
  counts_all[keep_for_size_factors, , drop = FALSE],
  data.frame(row.names = design$sample),
  ~ 1
)
dds <- estimateSizeFactors(dds)
size_factors <- sizeFactors(dds)
if (length(size_factors) != 114L || any(!is.finite(size_factors)) || any(size_factors <= 0)) {
  stop("Invalid DESeq2 size factors")
}

gene_list <- read.delim(gene_list_path, check.names = FALSE, stringsAsFactors = FALSE)
fatty_acid_genes <- unique(gene_list$GeneID)
missing_genes <- setdiff(fatty_acid_genes, rownames(counts_all))
if (length(fatty_acid_genes) != 176L || length(missing_genes) != 0L) {
  stop(
    "Unexpected fatty-acid gene catalogue: unique=", length(fatty_acid_genes),
    ", missing_from_raw_counts=", length(missing_genes)
  )
}

normalized <- sweep(
  counts_all[fatty_acid_genes, , drop = FALSE],
  2,
  size_factors,
  "/"
)
normalized_output <- data.frame(Gene_ID = rownames(normalized), normalized, check.names = FALSE)
safe_write(
  normalized_output,
  file.path(out_dir, "RNA_114_joint_normalized_FA_genes.tsv")
)

size_factor_output <- data.frame(
  sample = design$sample,
  client_id = design$client_id,
  variety = design$variety,
  replicate = design$replicate,
  delivery_phase = design$delivery_phase,
  delivery_time_raw = design$delivery_time_raw,
  DESeq2_size_factor = as.numeric(size_factors),
  stringsAsFactors = FALSE
)
safe_write(
  size_factor_output,
  file.path(out_dir, "RNA_114_joint_size_factors.tsv")
)

validation <- data.frame(
  check = c(
    "unique_samples", "genes_used_for_size_factors", "fatty_acid_genes",
    "missing_fatty_acid_genes", "finite_positive_size_factors"
  ),
  observed = c(
    nrow(design), sum(keep_for_size_factors), length(fatty_acid_genes),
    length(missing_genes), sum(is.finite(size_factors) & size_factors > 0)
  ),
  expected = c(114L, sum(keep_for_size_factors), 176L, 0L, 114L),
  status = c("PASS", "PASS", "PASS", "PASS", "PASS")
)
safe_write(validation, file.path(out_dir, "RNA_114_joint_normalization_validation.tsv"))

session <- capture.output(sessionInfo())
writeLines(session, file.path(out_dir, "RNA_114_joint_normalization_sessionInfo.txt"))
cat("Joint RNA normalization complete\n")
cat("Samples:", nrow(design), "\n")
cat("Fatty-acid genes:", length(fatty_acid_genes), "\n")
