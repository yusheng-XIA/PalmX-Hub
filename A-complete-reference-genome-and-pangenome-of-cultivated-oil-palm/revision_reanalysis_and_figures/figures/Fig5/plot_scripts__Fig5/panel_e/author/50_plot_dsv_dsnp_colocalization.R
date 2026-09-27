#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})

args <- commandArgs(trailingOnly = TRUE)
outdir <- if (length(args) >= 1) args[[1]] else "${ANALYSIS_DIR}/21_MS/06_result/dSVs/results-8.9/07_dsv_dsnp_colocalization/all38"
panel_display <- if (length(args) >= 2) args[[2]] else "All38"
fig_dir <- file.path(outdir, "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

read_tsv <- function(name) {
  read.delim(file.path(outdir, name), sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
}

nature_theme <- function(base_size = 7.5) {
  theme_classic(base_size = base_size, base_family = "Arial") +
    theme(
      plot.background = element_rect(fill = "white", colour = NA),
      panel.background = element_rect(fill = "white", colour = NA),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(colour = "#111111", linewidth = 0.34),
      axis.ticks = element_line(colour = "#111111", linewidth = 0.30),
      axis.ticks.length = unit(1.8, "pt"),
      axis.text = element_text(colour = "#111111", size = base_size - 1),
      axis.title = element_text(colour = "#000000", size = base_size),
      axis.title.x = element_text(margin = margin(t = 3)),
      axis.title.y = element_text(margin = margin(r = 3)),
      legend.title = element_blank(),
      legend.text = element_text(size = base_size - 1),
      legend.key.size = unit(3.2, "mm"),
      plot.tag = element_text(family = "Arial", face = "bold", size = 8),
      plot.margin = margin(4, 4, 3, 4)
    )
}

split_at_center <- function(d, group_col) {
  d <- d[order(d[[group_col]], d$Distance_Midpoint_kb), , drop = FALSE]
  side <- ifelse(d$Distance_Midpoint_kb < 0, "left", "right")
  d$Line_Group <- interaction(d[[group_col]], side, drop = TRUE, sep = "_")
  d
}

curve_legend_theme <- function() {
  theme(
    legend.position = c(0.985, 0.955),
    legend.justification = c(1, 1),
    legend.background = element_blank(),
    legend.key = element_blank(),
    legend.spacing.x = unit(1.5, "pt"),
    legend.margin = margin(0, 0, 0, 0)
  )
}

format_p <- function(p) {
  if (is.na(p)) {
    return("NA")
  }
  if (p < 2.2e-16) {
    return("< 2.2e-16")
  }
  if (p < 0.001) {
    return(format(p, scientific = TRUE, digits = 2))
  }
  sprintf("%.3f", p)
}

save_fig <- function(plot, name, width_mm, height_mm) {
  pdf_path <- file.path(fig_dir, paste0(name, ".pdf"))
  png_path <- file.path(fig_dir, paste0(name, ".png"))
  ggsave(pdf_path, plot = plot, width = width_mm, height = height_mm, units = "mm", device = cairo_pdf, bg = "white")
  png(filename = png_path, width = width_mm / 25.4, height = height_mm / 25.4, units = "in", res = 300, type = "cairo", bg = "white")
  print(plot)
  dev.off()
  data.frame(Figure = name, PDF = pdf_path, PNG = png_path, Width_mm = width_mm, Height_mm = height_mm, stringsAsFactors = FALSE)
}

counts <- read_tsv("fig_e_sample_chrom_counts.tsv")
corr <- read_tsv("fig_e_correlation_summary.tsv")
flank <- read_tsv("fig_f_window_summary.tsv")
flank_acc <- read_tsv("fig_f_accumulated_summary.tsv")
flank_acc_raw <- read_tsv("fig_f_accumulated_raw_noncarrier_summary.tsv")
obs_bg <- read_tsv("fig_f_observed_vs_matched_background_summary.tsv")
func_class <- read_tsv("fig_f_functional_class_summary.tsv")
svtype_class <- read_tsv("fig_f_svtype_class_summary.tsv")

counts$dSV_Count <- as.numeric(counts$dSV_Count)
counts$dSNP_Count <- as.numeric(counts$dSNP_Count)
flank$Distance_Midpoint_kb <- as.numeric(flank$Distance_Midpoint_kb)
flank$Mean_Accumulated_dSNP_Count <- as.numeric(flank$Mean_Accumulated_dSNP_Count)
flank$Bootstrap_CI_Low <- as.numeric(flank$Bootstrap_CI_Low)
flank$Bootstrap_CI_High <- as.numeric(flank$Bootstrap_CI_High)
flank$Phase <- factor(flank$Phase, levels = c("coupling", "repulsion"))
flank <- split_at_center(flank, "Phase")
for (obj_name in c("flank_acc", "flank_acc_raw")) {
  obj <- get(obj_name)
  obj$Distance_Midpoint_kb <- as.numeric(obj$Distance_Midpoint_kb)
  obj$Accumulated_dSNP_Count <- as.numeric(obj$Accumulated_dSNP_Count)
  obj$Bootstrap_CI_Low <- as.numeric(obj$Bootstrap_CI_Low)
  obj$Bootstrap_CI_High <- as.numeric(obj$Bootstrap_CI_High)
  obj$Phase <- factor(obj$Phase, levels = c("coupling", "repulsion"))
  assign(obj_name, obj)
}
for (obj_name in c("obs_bg", "func_class", "svtype_class")) {
  obj <- get(obj_name)
  obj$Distance_Midpoint_kb <- as.numeric(obj$Distance_Midpoint_kb)
  obj$Accumulated_dSNP_Count <- as.numeric(obj$Accumulated_dSNP_Count)
  obj$Mean_Per_Focal_DSV <- as.numeric(obj$Mean_Per_Focal_DSV)
  obj$Bootstrap_CI_Low <- as.numeric(obj$Bootstrap_CI_Low)
  obj$Bootstrap_CI_High <- as.numeric(obj$Bootstrap_CI_High)
  obj$Mean_CI_Low <- as.numeric(obj$Mean_CI_Low)
  obj$Mean_CI_High <- as.numeric(obj$Mean_CI_High)
  assign(obj_name, obj)
}

r_value <- as.numeric(corr$Pearson_r[1])
p_value <- as.numeric(corr$P_value[1])
n_value <- as.integer(corr$N[1])
p_text <- format_p(p_value)
p_label <- if (startsWith(p_text, "<")) paste("P", p_text) else paste("P =", p_text)
label_text <- sprintf("r = %.2f\n%s\nN = %d", r_value, p_label, n_value)

phase_cols <- c(coupling = "#E64B35", repulsion = "#6F7378")
phase_fills <- c(coupling = "#E64B35", repulsion = "#6F7378")

p_e <- ggplot(counts, aes(x = dSV_Count, y = dSNP_Count)) +
  geom_point(shape = 21, size = 1.55, stroke = 0.18, colour = "#222222", fill = "#4DBBD5", alpha = 0.76) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 0.55, colour = "#E64B35", fill = "#E64B35", alpha = 0.14) +
  annotate(
    "text",
    x = max(counts$dSV_Count, na.rm = TRUE) * 0.06,
    y = max(counts$dSNP_Count, na.rm = TRUE) * 0.97,
    label = label_text,
    hjust = 0,
    vjust = 1,
    family = "Arial",
    size = 3.0
  ) +
  labs(x = "dSV count per sample-chromosome", y = "dSNP count per sample-chromosome") +
  nature_theme()

p_f <- ggplot(flank, aes(x = Distance_Midpoint_kb, y = Mean_Accumulated_dSNP_Count, colour = Phase, fill = Phase)) +
  geom_line(aes(group = Line_Group), linewidth = 0.42, alpha = 0.92) +
  geom_vline(xintercept = 0, linewidth = 0.28, linetype = "22", colour = "#222222") +
  scale_colour_manual(values = phase_cols, labels = c(coupling = "Coupling", repulsion = "Repulsion")) +
  scale_fill_manual(values = phase_fills, labels = c(coupling = "Coupling", repulsion = "Repulsion")) +
  scale_x_continuous(limits = c(-1000, 1000), breaks = c(-1000, -500, 0, 500, 1000), expand = expansion(mult = c(0.01, 0.01))) +
  coord_cartesian(ylim = c(0, max(flank$Bootstrap_CI_High, na.rm = TRUE) * 1.05)) +
  labs(x = "Distance to dSV (kb)", y = "Mean accumulated dSNPs") +
  nature_theme() +
  curve_legend_theme()

make_accumulated_point_plot <- function(df, legend_labels = c(coupling = "Coupling", repulsion = "Repulsion")) {
  d <- split_at_center(df, "Phase")
  secondary <- d[d$Phase == "repulsion", , drop = FALSE]
  primary <- d[d$Phase != "repulsion", , drop = FALSE]
  ggplot() +
    geom_line(data = secondary, aes(x = Distance_Midpoint_kb, y = Accumulated_dSNP_Count, colour = Phase, group = Line_Group), linewidth = 0.24, alpha = 0.58, lineend = "round") +
    geom_point(data = secondary, aes(x = Distance_Midpoint_kb, y = Accumulated_dSNP_Count, colour = Phase), shape = 16, size = 0.58, alpha = 0.72) +
    geom_line(data = primary, aes(x = Distance_Midpoint_kb, y = Accumulated_dSNP_Count, colour = Phase, group = Line_Group), linewidth = 0.34, alpha = 0.90, lineend = "round") +
    geom_point(data = primary, aes(x = Distance_Midpoint_kb, y = Accumulated_dSNP_Count, colour = Phase), shape = 16, size = 0.82, alpha = 0.95) +
    geom_vline(xintercept = 0, linewidth = 0.28, linetype = "22", colour = "#222222") +
    scale_colour_manual(values = phase_cols, labels = legend_labels) +
    scale_x_continuous(limits = c(-1000, 1000), breaks = c(-1000, -500, 0, 500, 1000), expand = expansion(mult = c(0.01, 0.01))) +
    scale_y_continuous(breaks = function(x) pretty(x, n = 5), expand = expansion(mult = c(0.015, 0.045))) +
    labs(x = "Distance to dSV (kb)", y = "Accumulated dSNP count") +
    guides(colour = guide_legend(override.aes = list(linewidth = 0.45, size = 1.1, alpha = 1))) +
    nature_theme() +
    curve_legend_theme()
}

p_f_acc <- make_accumulated_point_plot(flank_acc)
p_f_acc_raw <- make_accumulated_point_plot(
  flank_acc_raw,
  c(coupling = "Coupling", repulsion = "Repulsion (all non-carriers)")
)

make_series_point_plot <- function(df, y_col, low_col, high_col, y_label, colours, labels, levels) {
  d <- df
  d$Series <- factor(d$Series, levels = levels)
  d$Plot_Y <- d[[y_col]]
  d$Plot_Low <- d[[low_col]]
  d$Plot_High <- d[[high_col]]
  d <- split_at_center(d, "Series")
  secondary <- d[d$Series %in% levels[-1], , drop = FALSE]
  primary <- d[d$Series %in% levels[1], , drop = FALSE]
  line_values <- setNames(c("solid", rep("22", max(0, length(levels) - 1))), levels)
  shape_values <- setNames(c(16, rep(1, max(0, length(levels) - 1))), levels)
  ggplot() +
    geom_line(data = secondary, aes(x = Distance_Midpoint_kb, y = Plot_Y, colour = Series, linetype = Series, group = Line_Group), linewidth = 0.24, alpha = 0.62, lineend = "round") +
    geom_point(data = secondary, aes(x = Distance_Midpoint_kb, y = Plot_Y, colour = Series, shape = Series), size = 0.60, alpha = 0.76) +
    geom_line(data = primary, aes(x = Distance_Midpoint_kb, y = Plot_Y, colour = Series, linetype = Series, group = Line_Group), linewidth = 0.34, alpha = 0.92, lineend = "round") +
    geom_point(data = primary, aes(x = Distance_Midpoint_kb, y = Plot_Y, colour = Series, shape = Series), size = 0.84, alpha = 0.96) +
    geom_vline(xintercept = 0, linewidth = 0.28, linetype = "22", colour = "#222222") +
    scale_colour_manual(values = colours, labels = labels, drop = FALSE) +
    scale_linetype_manual(values = line_values, labels = labels, drop = FALSE) +
    scale_shape_manual(values = shape_values, labels = labels, drop = FALSE) +
    scale_x_continuous(limits = c(-1000, 1000), breaks = c(-1000, -500, 0, 500, 1000), expand = expansion(mult = c(0.01, 0.01))) +
    scale_y_continuous(breaks = function(x) pretty(x, n = 5), expand = expansion(mult = c(0.015, 0.045))) +
    labs(x = "Distance to dSV (kb)", y = y_label) +
    guides(colour = guide_legend(override.aes = list(linewidth = 0.45, size = 1.1, alpha = 1))) +
    nature_theme() +
    curve_legend_theme()
}

obs_bg_levels <- c("Observed_dSV_carrier", "Matched_random_background")
obs_bg_cols <- c(Observed_dSV_carrier = "#D95F45", Matched_random_background = "#8B9096")
obs_bg_labels <- c(Observed_dSV_carrier = "Observed dSV carrier", Matched_random_background = "Matched background")
p_f_obs_bg <- make_series_point_plot(
  obs_bg,
  "Accumulated_dSNP_Count",
  "Bootstrap_CI_Low",
  "Bootstrap_CI_High",
  "Accumulated dSNP count",
  obs_bg_cols,
  obs_bg_labels,
  obs_bg_levels
)

class_levels <- c("CDS_supported", "Conserved_region_proxy_only")
class_n <- tapply(func_class$Effective_Focal_DSV_Count, func_class$Series, max)
class_labels <- c(
  CDS_supported = sprintf("CDS-supported (n=%s)", class_n[["CDS_supported"]]),
  Conserved_region_proxy_only = sprintf("Conserved-only (n=%s)", class_n[["Conserved_region_proxy_only"]])
)
class_cols <- c(CDS_supported = "#D95F45", Conserved_region_proxy_only = "#2E5E8E")
p_f_func <- make_series_point_plot(
  func_class[func_class$Series %in% class_levels, ],
  "Mean_Per_Focal_DSV",
  "Mean_CI_Low",
  "Mean_CI_High",
  "Mean accumulated dSNPs per focal dSV",
  class_cols,
  class_labels,
  class_levels
)

svtype_levels <- c("DEL", "INS")
svtype_n <- tapply(svtype_class$Effective_Focal_DSV_Count, svtype_class$Series, max)
svtype_labels <- c(
  DEL = sprintf("DEL (n=%s)", svtype_n[["DEL"]]),
  INS = sprintf("INS (n=%s)", svtype_n[["INS"]])
)
svtype_cols <- c(DEL = "#2E5E8E", INS = "#D9843B")
p_f_svtype <- make_series_point_plot(
  svtype_class[svtype_class$Series %in% svtype_levels, ],
  "Mean_Per_Focal_DSV",
  "Mean_CI_Low",
  "Mean_CI_High",
  "Mean accumulated dSNPs per focal dSV",
  svtype_cols,
  svtype_labels,
  svtype_levels
)

manifest <- list()
manifest[[length(manifest) + 1]] <- save_fig(p_e, "fig_dsv_dsnp_correlation_e", 89, 76)
manifest[[length(manifest) + 1]] <- save_fig(p_f_obs_bg, "dsv_observed", 183, 72)
manifest[[length(manifest) + 1]] <- save_fig(p_f_func, "dsv_function", 183, 72)
manifest[[length(manifest) + 1]] <- save_fig(p_f_svtype, "dsv_del_ins", 183, 72)

combined <- (p_e | p_f_obs_bg) +
  plot_layout(widths = c(0.95, 1.65)) +
  plot_annotation(tag_levels = list(c("e", "f"))) &
  nature_theme()
manifest[[length(manifest) + 1]] <- save_fig(combined, "fig_dsv_dsnp_ef_combined", 183, 78)

manifest_df <- do.call(rbind, manifest)
manifest_df$Panel <- panel_display
write.table(manifest_df, file.path(fig_dir, "figure_render_manifest.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

cat("[OK] Wrote Fig. e/f plots to ", fig_dir, "\n", sep = "")
