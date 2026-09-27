#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop("Usage: run_qc_preprocessing.R MODE XCMS_OUTPUT_DIR QC_OUTPUT_DIR")
}
mode <- args[[1L]]
input_dir <- normalizePath(args[[2L]], mustWork = TRUE)
output_dir <- args[[3L]]
if (!mode %in% c("neg", "pos")) stop("MODE must be neg or pos")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "figures"), showWarnings = FALSE)

suppressPackageStartupMessages({
  library(ggplot2)
  library(cowplot)
  library(matrixStats)
})

read_feature_matrix <- function(path) {
  x <- read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
  if (!identical(names(x)[[1L]], "feature_id") || anyDuplicated(x$feature_id)) {
    stop("Invalid feature matrix: ", path)
  }
  rn <- x$feature_id
  x$feature_id <- NULL
  ans <- as.matrix(x)
  storage.mode(ans) <- "double"
  rownames(ans) <- rn
  ans
}

manifest <- read.delim(file.path(input_dir, "sample_manifest.tsv"),
                       check.names = FALSE, stringsAsFactors = FALSE)
if (!("injection_order" %in% names(manifest)) &&
    "analysis_injection_order" %in% names(manifest)) {
  manifest$injection_order <- as.numeric(manifest$analysis_injection_order)
}
if (!("injection_order" %in% names(manifest)) ||
    any(!is.finite(manifest$injection_order)) ||
    anyDuplicated(manifest$injection_order)) {
  stop("Missing, non-finite, or duplicate analysis injection order")
}
definitions <- read.delim(file.path(input_dir, "aligned_feature_definitions.tsv"),
                          check.names = FALSE, stringsAsFactors = FALSE)
detected <- read_feature_matrix(file.path(input_dir,
                                           "aligned_feature_areas_detected.tsv"))
filled <- read_feature_matrix(file.path(input_dir,
                                         "aligned_feature_areas_filled.tsv"))
peak_counts <- read.delim(file.path(input_dir, "initial_chrom_peak_counts.tsv"),
                          check.names = FALSE, stringsAsFactors = FALSE)
rt_shift <- read.delim(file.path(input_dir, "retention_time_shift_summary.tsv"),
                       check.names = FALSE, stringsAsFactors = FALSE)
if (!("sample_name" %in% names(peak_counts)) ||
    !("sample_name" %in% names(rt_shift))) {
  stop("Peak-count or RT-shift table lacks sample_name")
}
peak_counts$injection_order <- manifest$injection_order[
  match(peak_counts$sample_name, manifest$sample_name)]
rt_shift$injection_order <- manifest$injection_order[
  match(rt_shift$sample_name, manifest$sample_name)]
if (any(!is.finite(peak_counts$injection_order)) ||
    any(!is.finite(rt_shift$injection_order))) {
  stop("Peak-count or RT-shift sample identity cannot be mapped to injection order")
}

sample_names <- manifest$sample_name
is_legacy <- "pair_id" %in% names(manifest) &&
  !("group_id" %in% names(manifest))
expected_samples <- if (is_legacy) 168L else 128L
if (nrow(manifest) != expected_samples ||
    ("mode" %in% names(manifest) && any(manifest$mode != mode)) ||
    !identical(colnames(detected), sample_names) ||
    !identical(colnames(filled), sample_names) ||
    !identical(rownames(detected), rownames(filled)) ||
    !setequal(definitions$feature_id, rownames(filled))) {
  stop("Input identities or dimensions do not match the validated full xcms output")
}

analysis_qc <- manifest$analysis_qc == "yes"
conditioning_qc <- manifest$conditioning_qc == "yes"
biological <- manifest$sample_type == "biological"
expected_analysis_qc <- if (is_legacy) 15L else 11L
expected_conditioning_qc <- if (is_legacy) 1L else 3L
expected_biological <- if (is_legacy) 152L else 114L
if (sum(analysis_qc) != expected_analysis_qc ||
    sum(conditioning_qc) != expected_conditioning_qc ||
    sum(biological) != expected_biological) {
  stop("Unexpected analysis-QC, conditioning-QC, or biological sample count")
}
biological_group <- if (is_legacy) manifest$pair_id else manifest$group_id
qc_detection_minimum <- ceiling(0.8 * sum(analysis_qc))
biological_group_detection_minimum <- 2L

present <- !is.na(detected) & detected > 0
qc_detected_n <- rowSums(present[, analysis_qc, drop = FALSE])
bio_group_detected_max <- vapply(seq_len(nrow(present)), function(i) {
  max(vapply(split(which(biological), biological_group[biological]), function(j) {
    sum(present[i, j])
  }, integer(1)))
}, integer(1))
presence_pass <- qc_detected_n >= qc_detection_minimum &
  bio_group_detected_max >= biological_group_detection_minimum
if (sum(presence_pass) < 100L) stop("Too few features pass presence filtering")

work <- filled[presence_pass, , drop = FALSE]
work[!is.finite(work) | work <= 0] <- NA_real_
qc_reference <- rowMedians(work[, analysis_qc, drop = FALSE], na.rm = TRUE)
valid_reference <- is.finite(qc_reference) & qc_reference > 0
size_factor <- vapply(seq_len(ncol(work)), function(j) {
  ratios <- work[valid_reference, j] / qc_reference[valid_reference]
  ratios <- ratios[is.finite(ratios) & ratios > 0]
  if (length(ratios) < 100L) return(NA_real_)
  median(ratios)
}, numeric(1))
if (any(!is.finite(size_factor) | size_factor <= 0)) {
  stop("Pooled-QC median-fold normalization failed")
}
size_factor <- size_factor / median(size_factor[analysis_qc])
normalized <- sweep(work, 2L, size_factor, "/")

qc_order <- manifest$injection_order[analysis_qc]
all_order <- manifest$injection_order
drift_one <- function(y) {
  qy <- log2(y[analysis_qc])
  good <- is.finite(qy)
  if (sum(good) < 6L) return(y)
  fit_data <- data.frame(x = qc_order[good], y = qy[good])
  fit <- try(loess(y ~ x, data = fit_data, span = 0.75, degree = 2,
                   family = "symmetric", surface = "direct",
                   control = loess.control(iterations = 3)), silent = TRUE)
  if (inherits(fit, "try-error")) return(y)
  fitted_qc <- as.numeric(predict(fit, newdata = data.frame(x = qc_order[good])))
  if (sum(is.finite(fitted_qc)) < 4L) return(y)
  trend <- approx(qc_order[good][is.finite(fitted_qc)],
                  fitted_qc[is.finite(fitted_qc)], xout = all_order,
                  rule = 2, ties = "ordered")$y
  target <- median(qy[good], na.rm = TRUE)
  y / (2^(trend - target))
}
workers <- min(8L, parallel::detectCores(logical = FALSE))
corrected_results <- parallel::mclapply(seq_len(nrow(normalized)), function(i) {
  tryCatch(
    list(value = drift_one(normalized[i, ]), error = ""),
    error = function(e) list(value = normalized[i, ], error = conditionMessage(e))
  )
}, mc.cores = workers, mc.preschedule = TRUE)
drift_errors <- vapply(corrected_results, `[[`, character(1), "error")
if (any(nzchar(drift_errors))) {
  error_table <- as.data.frame(table(drift_errors[nzchar(drift_errors)]),
                               stringsAsFactors = FALSE)
  names(error_table) <- c("error_message", "feature_count")
} else {
  error_table <- data.frame(error_message = character(0),
                            feature_count = integer(0))
}
write.table(error_table, file.path(output_dir, "drift_correction_error_summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
error_n <- sum(nzchar(drift_errors))
if (error_n > max(10L, ceiling(0.01 * nrow(normalized)))) {
  stop("Drift correction errors exceed 1%: ", error_n, "/", nrow(normalized))
}
corrected <- do.call(rbind, lapply(corrected_results, `[[`, "value"))
rownames(corrected) <- rownames(normalized)
colnames(corrected) <- colnames(normalized)

rsd_percent <- function(x) {
  good <- is.finite(x) & x > 0
  if (sum(good) < 3L || mean(x[good]) <= 0) return(NA_real_)
  100 * sd(x[good]) / mean(x[good])
}
qc_rsd_before <- apply(normalized[, analysis_qc, drop = FALSE], 1L, rsd_percent)
qc_rsd_after <- apply(corrected[, analysis_qc, drop = FALSE], 1L, rsd_percent)
rsd_pass <- is.finite(qc_rsd_after) & qc_rsd_after <= 30
final_feature_ids <- rownames(corrected)[rsd_pass]
if (length(final_feature_ids) < 100L) stop("Too few features pass QC RSD filtering")
final_matrix <- corrected[rsd_pass, , drop = FALSE]

metrics <- data.frame(
  feature_id = rownames(filled),
  qc_detected_n = qc_detected_n,
  qc_detected_fraction = qc_detected_n / sum(analysis_qc),
  maximum_biological_group_detected_n = bio_group_detected_max,
  presence_pass = presence_pass,
  qc_rsd_before_percent = NA_real_,
  qc_rsd_after_percent = NA_real_,
  qc_rsd_pass = FALSE,
  final_pass = FALSE,
  stringsAsFactors = FALSE
)
idx <- match(rownames(normalized), metrics$feature_id)
metrics$qc_rsd_before_percent[idx] <- qc_rsd_before
metrics$qc_rsd_after_percent[idx] <- qc_rsd_after
metrics$qc_rsd_pass[idx] <- rsd_pass
metrics$final_pass <- metrics$presence_pass & metrics$qc_rsd_pass
metrics <- merge(definitions, metrics, by = "feature_id", all.y = TRUE, sort = FALSE)
metrics <- metrics[match(rownames(filled), metrics$feature_id), , drop = FALSE]

write.table(metrics, file.path(output_dir, "feature_qc_metrics.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE, na = "")
write.table(data.frame(feature_id = rownames(final_matrix), final_matrix,
                       check.names = FALSE),
            file.path(output_dir,
                      "normalized_drift_corrected_qc_filtered_matrix.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE, na = "")
safe_manifest_column <- function(name, fallback = "") {
  if (name %in% names(manifest)) manifest[[name]] else rep(fallback, nrow(manifest))
}
write.table(data.frame(sample_name = sample_names,
                       injection_order = manifest$injection_order,
                       sample_type = manifest$sample_type,
                       conditioning_qc = manifest$conditioning_qc,
                       analysis_qc = manifest$analysis_qc,
                       group_id = biological_group,
                       material = if (is_legacy) manifest$series else manifest$material,
                       time_value = if (is_legacy) manifest$series_number else manifest$time_value,
                       time_unit = safe_manifest_column("time_unit"),
                       normalization_factor = size_factor,
                       detected_missing_fraction = colMeans(!present),
                       final_missing_fraction = colMeans(!is.finite(final_matrix) |
                                                           final_matrix <= 0),
                       final_total_area = colSums(final_matrix, na.rm = TRUE)),
            file.path(output_dir, "sample_qc_metrics.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE, na = "")

flow <- data.frame(
  stage = c("aligned_features",
            sprintf("QC detected in at least %d of %d",
                    qc_detection_minimum, sum(analysis_qc)),
            sprintf("detected in at least %d samples in any biological group",
                    biological_group_detection_minimum),
            "combined presence filter", "post-correction QC RSD <= 30%"),
  feature_count = c(nrow(filled), sum(qc_detected_n >= qc_detection_minimum),
                    sum(bio_group_detected_max >= biological_group_detection_minimum),
                    sum(presence_pass),
                    nrow(final_matrix))
)
write.table(flow, file.path(output_dir, "feature_filter_flow.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE)
write.table(data.frame(
  parameter = c("analysis_QC_count", "conditioning_QC_count",
                "QC_detection_minimum", "biological_group_detection_minimum",
                "normalization", "drift_correction", "drift_loess_span",
                "post_correction_QC_RSD_maximum_percent", "PCA_transform",
                "PCA_scaling"),
  value = c(sum(analysis_qc), sum(conditioning_qc), qc_detection_minimum,
            biological_group_detection_minimum,
            "median fold to pooled analysis-QC reference",
            "robust QC-RLSC; conditioning QC excluded from fit", 0.75, 30,
            "log2 after half-minimum imputation for PCA only", "Pareto")
), file.path(output_dir, "preprocessing_parameters.tsv"), sep = "\t",
quote = FALSE, row.names = FALSE)

pca_matrix <- final_matrix
for (i in seq_len(nrow(pca_matrix))) {
  good <- is.finite(pca_matrix[i, ]) & pca_matrix[i, ] > 0
  replacement <- if (any(good)) min(pca_matrix[i, good]) / 2 else 1
  pca_matrix[i, !good] <- replacement
}
log_matrix <- log2(pca_matrix)
feature_sd <- rowSds(log_matrix)
pca_keep <- is.finite(feature_sd) & feature_sd > 0
x <- t(log_matrix[pca_keep, , drop = FALSE])
x <- scale(x, center = TRUE, scale = sqrt(feature_sd[pca_keep]))
pca <- prcomp(x, center = FALSE, scale. = FALSE)
variance <- 100 * pca$sdev^2 / sum(pca$sdev^2)
biological_class <- if (is_legacy) manifest$series else paste0("M", manifest$material)
sample_class <- ifelse(conditioning_qc, "Conditioning QC",
                       ifelse(analysis_qc, "Analysis QC", biological_class))
pca_scores <- data.frame(
  sample_name = sample_names,
  injection_order = manifest$injection_order,
  sample_class = sample_class,
  group_id = biological_group,
  PC1 = pca$x[, 1L], PC2 = pca$x[, 2L], PC3 = pca$x[, 3L]
)
write.table(pca_scores, file.path(output_dir, "pca_scores.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE, na = "")
loadings <- data.frame(feature_id = rownames(pca$rotation),
                       PC1 = pca$rotation[, 1L], PC2 = pca$rotation[, 2L],
                       PC3 = pca$rotation[, 3L])
write.table(loadings, file.path(output_dir, "pca_loadings.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE)
write.table(data.frame(component = paste0("PC", seq_along(variance)),
                       variance_percent = variance),
            file.path(output_dir, "pca_variance.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE)

theme_qc <- theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(), legend.position = "right")
qc_kind <- ifelse(conditioning_qc, "Conditioning QC",
                  ifelse(analysis_qc, "Analysis QC", "Biological"))
plot_meta <- data.frame(manifest, qc_kind = qc_kind,
                        normalization_factor = size_factor,
                        final_missing_fraction = colMeans(!is.finite(final_matrix) |
                                                            final_matrix <= 0))
peak_counts$qc_kind <- qc_kind[match(peak_counts$sample_name, sample_names)]
rt_shift$qc_kind <- qc_kind[match(rt_shift$sample_name, sample_names)]
pal <- c("Biological" = "#0072B2", "Analysis QC" = "#222222",
         "Conditioning QC" = "#E69F00")

p_peak <- ggplot(peak_counts, aes(injection_order, chrom_peak_count,
                                  color = qc_kind)) +
  geom_line(aes(group = 1), color = "grey75", linewidth = 0.35) +
  geom_point(size = 1.5) + scale_color_manual(values = pal) + theme_qc +
  labs(x = "Injection order", y = "Detected chromatographic peaks",
       color = NULL, title = paste(toupper(mode), "initial peak counts"))
p_rt <- ggplot(rt_shift, aes(injection_order, median_shift_seconds,
                             color = qc_kind)) +
  geom_hline(yintercept = 0, color = "grey70", linewidth = 0.35) +
  geom_line(aes(group = 1), color = "grey75", linewidth = 0.35) +
  geom_point(size = 1.5) + scale_color_manual(values = pal) + theme_qc +
  labs(x = "Injection order", y = "Median RT correction (s)", color = NULL,
       title = "Retention-time alignment")
p_norm <- ggplot(plot_meta, aes(injection_order, normalization_factor,
                                color = qc_kind)) +
  geom_hline(yintercept = 1, color = "grey70", linewidth = 0.35) +
  geom_line(aes(group = 1), color = "grey75", linewidth = 0.35) +
  geom_point(size = 1.5) + scale_color_manual(values = pal) + theme_qc +
  labs(x = "Injection order", y = "Median-fold normalization factor",
       color = NULL, title = "Signal normalization")
rsd_plot <- rbind(data.frame(RSD = qc_rsd_before, stage = "Before QC-RLSC"),
                  data.frame(RSD = qc_rsd_after, stage = "After QC-RLSC"))
rsd_plot <- rsd_plot[is.finite(rsd_plot$RSD) & rsd_plot$RSD <= 100, ]
p_rsd <- ggplot(rsd_plot, aes(RSD, fill = stage)) +
  geom_histogram(position = "identity", bins = 40, alpha = 0.55) +
  geom_vline(xintercept = 30, linetype = 2, color = "#D55E00") +
  scale_fill_manual(values = c("Before QC-RLSC" = "grey70",
                               "After QC-RLSC" = "#009E73")) + theme_qc +
  labs(x = "Analysis-QC RSD (%)", y = "Features", fill = NULL,
       title = "QC precision")
p_missing <- ggplot(plot_meta, aes(injection_order, 100 * final_missing_fraction,
                                   color = qc_kind)) +
  geom_line(aes(group = 1), color = "grey75", linewidth = 0.35) +
  geom_point(size = 1.5) + scale_color_manual(values = pal) + theme_qc +
  labs(x = "Injection order", y = "Missing values after QC filter (%)",
       color = NULL, title = "Sample completeness")
class_palette <- c(
  "M1" = "#0072B2", "M4" = "#D55E00",
  "FL" = "#0072B2", "NS" = "#D55E00",
  "TK" = "#009E73", "TN" = "#CC79A7",
  "Analysis QC" = "#222222", "Conditioning QC" = "#E69F00"
)
class_shapes <- c(
  "M1" = 16, "M4" = 17, "FL" = 16, "NS" = 17,
  "TK" = 15, "TN" = 18, "Analysis QC" = 8, "Conditioning QC" = 4
)
p_pca <- ggplot(pca_scores, aes(PC1, PC2, color = sample_class,
                                shape = sample_class)) +
  geom_point(size = 2.1, alpha = 0.85) +
  scale_color_manual(values = class_palette) +
  scale_shape_manual(values = class_shapes) +
  theme_qc + labs(x = sprintf("PC1 (%.1f%%)", variance[[1L]]),
                  y = sprintf("PC2 (%.1f%%)", variance[[2L]]),
                  color = NULL, shape = NULL, title = "Unsupervised PCA")

overview <- plot_grid(p_peak, p_rt, p_norm, p_rsd, p_missing, p_pca,
                      ncol = 2, labels = LETTERS[1:6], align = "hv")
ggsave(file.path(output_dir, "figures", paste0(mode, "_QC_overview.pdf")),
       overview, width = 12, height = 14, units = "in", device = cairo_pdf)
ggsave(file.path(output_dir, "figures", paste0(mode, "_QC_overview.png")),
       overview, width = 12, height = 14, units = "in", dpi = 300)
ggsave(file.path(output_dir, "figures", paste0(mode, "_PCA.pdf")),
       p_pca, width = 7, height = 5.5, units = "in", device = cairo_pdf)

summary_table <- data.frame(
  metric = c("samples", "biological_samples", "conditioning_QC",
             "analysis_QC", "aligned_features", "presence_pass_features",
             "final_QC_pass_features", "median_QC_RSD_before_percent",
             "median_QC_RSD_after_percent", "PC1_variance_percent",
             "PC2_variance_percent"),
  value = c(nrow(manifest), sum(biological), sum(conditioning_qc),
            sum(analysis_qc), nrow(filled), sum(presence_pass),
            nrow(final_matrix), median(qc_rsd_before, na.rm = TRUE),
            median(qc_rsd_after, na.rm = TRUE), variance[[1L]], variance[[2L]])
)
write.table(summary_table, file.path(output_dir, "qc_summary.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(output_dir, "R_sessionInfo.txt"))
writeLines("0", file.path(output_dir, "validation_status.txt"))
