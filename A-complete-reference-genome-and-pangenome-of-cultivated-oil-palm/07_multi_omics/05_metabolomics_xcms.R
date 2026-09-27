#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) {
  stop("Usage: run_xcms_full.R MODE FILE_LIST MANIFEST OUTPUT_DIR")
}
mode <- args[[1L]]
file_list <- args[[2L]]
manifest_file <- args[[3L]]
output_dir <- args[[4L]]
if (!mode %in% c("neg", "pos")) stop("MODE must be neg or pos")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(xcms)
  library(MSnbase)
  library(BiocParallel)
})

files <- readLines(file_list, warn = FALSE)
files <- files[nzchar(files)]
manifest <- read.delim(manifest_file, check.names = FALSE, stringsAsFactors = FALSE)
sample_name <- sub("\\.mzML$", "", basename(files))
if (length(files) != 128L || nrow(manifest) != 128L || any(!file.exists(files))) {
  stop("Expected 128 existing mzML files and 128 manifest rows")
}
if (!identical(sample_name, manifest$sample_name) || any(manifest$mode != mode) ||
    anyDuplicated(sample_name)) {
  stop("Input files do not exactly match the mode-specific injection-order manifest")
}
sample_group <- ifelse(manifest$sample_type == "QC", "QC", manifest$group_id)
if (any(!nzchar(sample_group)) || length(unique(sample_group)) != 39L) {
  stop("Expected 38 biological groups plus one QC group")
}

writeLines(capture.output(sessionInfo()), file.path(output_dir, "R_sessionInfo.txt"))
write.table(manifest, file.path(output_dir, "sample_manifest.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE, na = "")
timings <- data.frame(stage = character(), elapsed_seconds = numeric())
timed <- function(stage, expr) {
  start <- proc.time()[["elapsed"]]
  value <- force(expr)
  timings <<- rbind(timings, data.frame(
    stage = stage, elapsed_seconds = proc.time()[["elapsed"]] - start
  ))
  write.table(timings, file.path(output_dir, "stage_timings.tsv"), sep = "\t",
              quote = FALSE, row.names = FALSE)
  value
}

pdata <- Biobase::AnnotatedDataFrame(data.frame(
  sample_name = sample_name,
  sample_group = sample_group,
  injection_order = manifest$injection_order,
  sample_type = manifest$sample_type,
  row.names = sample_name
))
raw_data <- timed("read_mzml", readMSData(files = files, pdata = pdata, mode = "onDisk"))
bp <- MulticoreParam(workers = 16L, progressbar = TRUE, stop.on.error = TRUE)
register(bp, default = TRUE)
centwave <- CentWaveParam(
  ppm = 15, peakwidth = c(3, 30), snthresh = 8,
  prefilter = c(5, 1000), mzCenterFun = "wMean", integrate = 1,
  mzdiff = -0.001, fitgauss = FALSE, noise = 1000,
  verboseColumns = TRUE
)
xdata <- timed("find_chrom_peaks",
               findChromPeaks(raw_data, param = centwave, msLevel = 1L,
                              BPPARAM = bp))
saveRDS(xdata, file.path(output_dir, "checkpoint_01_chrom_peaks.rds"),
        compress = FALSE)

peak_sample <- as.integer(chromPeaks(xdata)[, "sample"])
initial_counts <- manifest
initial_counts$chrom_peak_count <- tabulate(peak_sample, nbins = length(files))
write.table(initial_counts, file.path(output_dir, "initial_chrom_peak_counts.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE, na = "")

reference_index <- which(manifest$sample_type == "QC" &
                         manifest$sample_name == paste0("DJ272-ZX01-0301-", mode, "-QC-7"))
if (length(reference_index) != 1L) stop("Expected exactly one QC-7 reference")
obiwarp <- ObiwarpParam(
  binSize = 0.1, centerSample = reference_index,
  response = 1, distFun = "cor_opt", factorDiag = 2, factorGap = 1,
  localAlignment = FALSE, rtimeDifferenceThreshold = 2
)
xdata <- timed("adjust_rtime_obiwarp",
               adjustRtime(xdata, param = obiwarp))
saveRDS(xdata, file.path(output_dir, "checkpoint_02_adjusted_rtime.rds"),
        compress = FALSE)

raw_rt <- rtime(xdata, adjusted = FALSE, bySample = TRUE)
adj_rt <- adjustedRtime(xdata, bySample = TRUE)
rt_summary <- do.call(rbind, lapply(seq_along(files), function(i) {
  z <- adj_rt[[i]] - raw_rt[[i]]
  data.frame(
    sample_name = sample_name[[i]],
    injection_order = manifest$injection_order[[i]],
    sample_type = manifest$sample_type[[i]],
    conditioning_qc = manifest$conditioning_qc[[i]],
    analysis_qc = manifest$analysis_qc[[i]],
    median_shift_seconds = median(z, na.rm = TRUE),
    q05_shift_seconds = unname(quantile(z, 0.05, na.rm = TRUE)),
    q95_shift_seconds = unname(quantile(z, 0.95, na.rm = TRUE)),
    max_abs_shift_seconds = max(abs(z), na.rm = TRUE)
  )
}))
write.table(rt_summary, file.path(output_dir, "retention_time_shift_summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

group_param <- PeakDensityParam(
  sampleGroups = as.integer(factor(sample_group, levels = unique(sample_group))),
  bw = 5, minFraction = 0.5, minSamples = 2,
  binSize = 0.01, ppm = 15, maxFeatures = 100
)
xdata <- timed("group_chrom_peaks",
               groupChromPeaks(xdata, param = group_param))
saveRDS(xdata, file.path(output_dir, "checkpoint_03_grouped_features.rds"),
        compress = FALSE)
xdata <- timed("fill_chrom_peaks",
               fillChromPeaks(xdata, param = FillChromPeaksParam(ppm = 15),
                              BPPARAM = bp))
saveRDS(xdata, file.path(output_dir, "checkpoint_04_filled_features.rds"),
        compress = FALSE)

detected_matrix <- featureValues(xdata, value = "into", method = "medret",
                                 filled = FALSE)
filled_matrix <- featureValues(xdata, value = "into", method = "medret",
                               filled = TRUE)
detected_sample_name <- sub("\\.mzML$", "", colnames(detected_matrix))
filled_sample_name <- sub("\\.mzML$", "", colnames(filled_matrix))
if (anyDuplicated(detected_sample_name) || anyDuplicated(filled_sample_name) ||
    !setequal(detected_sample_name, sample_name) ||
    !setequal(filled_sample_name, sample_name) ||
    !identical(rownames(detected_matrix), rownames(filled_matrix))) {
  stop("Feature matrices do not map one-to-one to feature and sample identities")
}
detected_matrix <- detected_matrix[, match(sample_name, detected_sample_name), drop = FALSE]
filled_matrix <- filled_matrix[, match(sample_name, filled_sample_name), drop = FALSE]
colnames(detected_matrix) <- sample_name
colnames(filled_matrix) <- sample_name
feature_def <- as.data.frame(featureDefinitions(xdata))
feature_def$feature_id <- rownames(feature_def)
feature_def$n_detected_samples <- rowSums(!is.na(detected_matrix) & detected_matrix > 0)
feature_def$n_filled_samples <- rowSums(!is.na(filled_matrix) & filled_matrix > 0)
feature_def$peakidx <- NULL
feature_def <- feature_def[, c("feature_id", setdiff(names(feature_def), "feature_id")),
                           drop = FALSE]
write.table(feature_def, file.path(output_dir, "aligned_feature_definitions.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE, na = "")
write.table(data.frame(feature_id = rownames(detected_matrix), detected_matrix,
                       check.names = FALSE),
            file.path(output_dir, "aligned_feature_areas_detected.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE, na = "")
write.table(data.frame(feature_id = rownames(filled_matrix), filled_matrix,
                       check.names = FALSE),
            file.path(output_dir, "aligned_feature_areas_filled.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE, na = "")

detected_missing <- rowMeans(is.na(detected_matrix) | detected_matrix <= 0)
filled_missing <- rowMeans(is.na(filled_matrix) | filled_matrix <= 0)
feature_summary <- data.frame(
  metric = c("aligned_features", "detected_complete_all_samples",
             "detected_present_at_least_half_samples", "median_detected_missing_fraction",
             "median_filled_missing_fraction", "minimum_sample_chrom_peaks",
             "median_sample_chrom_peaks", "maximum_sample_chrom_peaks"),
  value = c(nrow(filled_matrix), sum(detected_missing == 0),
            sum(detected_missing <= 0.5), median(detected_missing),
            median(filled_missing), min(initial_counts$chrom_peak_count),
            median(initial_counts$chrom_peak_count), max(initial_counts$chrom_peak_count))
)
write.table(feature_summary, file.path(output_dir, "feature_summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
saveRDS(xdata, file.path(output_dir, "xcms_full_processed.rds"), compress = FALSE)

if (nrow(filled_matrix) < 100L || any(initial_counts$chrom_peak_count < 100L) ||
    any(!is.finite(rt_summary$max_abs_shift_seconds))) {
  stop("Full xcms structural validation failed")
}
writeLines("0", file.path(output_dir, "validation_status.txt"))
