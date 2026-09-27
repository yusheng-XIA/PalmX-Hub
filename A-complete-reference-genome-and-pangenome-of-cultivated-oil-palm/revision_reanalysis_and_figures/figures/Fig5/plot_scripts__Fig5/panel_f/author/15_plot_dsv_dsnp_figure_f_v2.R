#!/usr/bin/env Rscript
# Fig. f v2 in the "1.png" agronomy-journal style, ONE connected figure of two
# tightly-stacked sub-panels:
#   TOP   : horizontal boxplot summary of the per-focal-dSV total accumulated dSNP
#           (integral over +/-1 Mb) for the two series, red/blue, Wilcoxon bracket.
#   BOTTOM: the distance profile (mean +/- SE), small filled markers + line, on a
#           soft light-blue vertical gradient with a full black box frame and a grey
#           dashed reference line at the dSV centre. The bottom y-axis uses a
#           pseudo-log (base 10) transform so the huge central peak and the low
#           background are both visible.
# Standalone; does not modify the v1 figures.

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(grid)
  library(scales)
})

args <- commandArgs(trailingOnly = TRUE)
outdir <- if (length(args) >= 1) args[[1]] else "${ANALYSIS_DIR}/21_MS/06_result/dSVs/results/07_dsv_dsnp_figure_ef"
fig_dir <- file.path(outdir, "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

read_tsv <- function(name) {
  read.delim(file.path(outdir, name), sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
}

grad_grob <- rasterGrob(
  matrix(colorRampPalette(c("#F1F7FC", "#C2D8EE"))(256), ncol = 1),
  width = unit(1, "npc"), height = unit(1, "npc"), interpolate = TRUE
)

sig_stars <- function(p) {
  if (is.na(p)) return("ns")
  if (p < 1e-3) return("***")
  if (p < 1e-2) return("**")
  if (p < 5e-2) return("*")
  "ns"
}

theme_1png <- function(base_size = 9) {
  theme_bw(base_size = base_size, base_family = "Arial") +
    theme(
      panel.grid = element_blank(),
      panel.background = element_blank(),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.7),
      axis.text = element_text(colour = "black", size = base_size - 0.5),
      axis.title = element_text(colour = "black", face = "bold", size = base_size + 1),
      axis.title.x = element_text(margin = margin(t = 2)),
      axis.title.y = element_text(margin = margin(r = 3)),
      axis.ticks = element_line(colour = "black", linewidth = 0.5),
      axis.ticks.length = unit(2.2, "pt"),
      legend.title = element_blank(),
      legend.text = element_text(size = base_size - 0.5, colour = "black"),
      legend.background = element_rect(fill = scales::alpha("white", 0.72), colour = NA),
      legend.key = element_blank(),
      legend.key.size = unit(3.6, "mm"),
      plot.tag = element_text(face = "bold", size = base_size + 4),
      plot.tag.position = c(0.012, 0.86)
    )
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

# ---- TOP: horizontal raincloud (half-violin cloud + boxplot + jittered rain) ----
# group_levels ordered bottom-to-top on the y axis (last element drawn on top).
# Only the extreme top 5% of the top group is trimmed (from the cloud/rain display),
# so most large values are kept; the significance test uses the full uncapped data.
make_raincloud <- function(totals, group_levels, cols, ylabels, xlab, paired, tag = "f") {
  set.seed(42)
  totals$Group <- factor(totals$Group, levels = group_levels)
  g_top <- group_levels[length(group_levels)]; g_bot <- group_levels[1]
  # cap = 95th percentile of the top group -> only its extreme top ~5% is dropped from display
  cap <- max(as.numeric(quantile(totals$Total[totals$Group == g_top], 0.95, na.rm = TRUE)), 1)

  v_top <- totals$Total[totals$Group == g_top]; v_bot <- totals$Total[totals$Group == g_bot]
  if (paired) {
    m <- merge(data.frame(id = totals$ID[totals$Group == g_top], top = v_top),
               data.frame(id = totals$ID[totals$Group == g_bot], bot = v_bot), by = "id")
    pv <- suppressWarnings(wilcox.test(m$top, m$bot, paired = TRUE)$p.value)
  } else {
    pv <- suppressWarnings(wilcox.test(v_top, v_bot)$p.value)
  }
  stars <- sig_stars(pv)

  d <- totals[totals$Total <= cap, , drop = FALSE]        # display data (top 5% of top group dropped)
  d$gy <- as.integer(d$Group)                             # 1..k numeric y positions
  d$ry <- d$gy - 0.30 + runif(nrow(d), -0.10, 0.10)       # rain (jittered points) below the line
  d$by <- d$gy - 0.06                                     # boxplot just below the cloud

  clouds <- do.call(rbind, lapply(levels(d$Group), function(g) {
    v <- d$Total[d$Group == g]; gy <- match(g, levels(d$Group))
    de <- density(v, from = 0, to = cap)
    h <- de$y / max(de$y) * 0.40
    data.frame(x = c(de$x, rev(de$x)), y = c(gy + 0.06 + h, rep(gy + 0.06, length(de$x))), Group = g)
  }))
  clouds$Group <- factor(clouds$Group, levels = group_levels)

  upper <- cap * 1.15; xr <- cap * 1.06; xr_l <- xr - cap * 0.05
  xcand <- c(0, 5, 10, 20, 40, 60, 90, 120, 160, 200)
  xbrks <- xcand[xcand <= upper]
  k <- length(group_levels)
  ggplot() +
    geom_polygon(data = clouds, aes(x = x, y = y, fill = Group, group = Group),
                 colour = "grey45", linewidth = 0.25, alpha = 0.5) +
    geom_point(data = d, aes(x = Total, y = ry, colour = Group), size = 0.35, alpha = 0.22, stroke = 0) +
    stat_boxplot(data = d, aes(x = Total, y = by, group = Group), geom = "errorbar",
                 orientation = "y", width = 0.10, linewidth = 0.38, colour = "grey20") +
    geom_boxplot(data = d, aes(x = Total, y = by, group = Group), orientation = "y",
                 width = 0.16, linewidth = 0.42, outlier.shape = NA, fill = "white", colour = "grey20") +
    annotate("segment", x = xr, xend = xr, y = 1, yend = k, linewidth = 0.5, colour = "black") +
    annotate("segment", x = xr_l, xend = xr, y = 1, yend = 1, linewidth = 0.5, colour = "black") +
    annotate("segment", x = xr_l, xend = xr, y = k, yend = k, linewidth = 0.5, colour = "black") +
    annotate("text", x = xr + cap * 0.05, y = (1 + k) / 2, label = stars, angle = 90, fontface = "bold", size = 3.4) +
    scale_fill_manual(values = cols) +
    scale_colour_manual(values = cols) +
    scale_y_continuous(breaks = seq_len(k), labels = ylabels[group_levels],
                       limits = c(0.55, k + 0.55), expand = expansion(mult = c(0, 0))) +
    scale_x_continuous(trans = "sqrt", breaks = xbrks, expand = expansion(mult = c(0.01, 0.02))) +
    coord_cartesian(xlim = c(0, upper), clip = "on") +
    labs(x = xlab, y = NULL, tag = tag) +
    theme_1png() +
    theme(legend.position = "none",
          axis.text.y = element_text(face = "bold"),
          axis.ticks.y = element_blank(),
          plot.margin = margin(6, 9, 1, 5))
}

# ---- BOTTOM: distance profile, mean +/- SE, small markers, optional pseudo-log y ----
make_profile <- function(df, y_col, low_col, high_col, levels, cols, labels, y_label,
                         log_y = FALSE, y_breaks = waiver()) {
  d <- df
  d$Series <- factor(d$Series, levels = levels)
  d$X <- as.numeric(d$Distance_Midpoint_kb)
  d$Y <- as.numeric(d[[y_col]])
  d$SE <- (as.numeric(d[[high_col]]) - as.numeric(d[[low_col]])) / 3.92  # ~1 SE from 95% bootstrap CI
  d$Lo <- pmax(d$Y - d$SE, 0)
  d$Hi <- d$Y + d$SE
  d <- d[order(d$Series, d$X), , drop = FALSE]
  ymax <- max(d$Hi, na.rm = TRUE) * 1.08

  p <- ggplot(d, aes(x = X, y = Y, colour = Series, fill = Series)) +
    annotation_custom(grad_grob, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) +
    geom_vline(xintercept = 0, linetype = "22", colour = "grey45", linewidth = 0.4) +
    geom_errorbar(aes(ymin = Lo, ymax = Hi), width = 6, linewidth = 0.22, alpha = 0.5) +
    geom_line(aes(group = Series), linewidth = 0.45) +
    geom_point(shape = 21, colour = "grey20", size = 0.85, stroke = 0.22) +
    scale_colour_manual(values = cols, labels = labels) +
    scale_fill_manual(values = cols, labels = labels) +
    scale_x_continuous(limits = c(-1000, 1000), breaks = c(-1000, -500, 0, 500, 1000),
                       expand = expansion(mult = c(0.012, 0.012))) +
    labs(x = "Distance to dSV (kb)", y = y_label) +
    theme_1png() +
    theme(legend.position = c(0.99, 0.97), legend.justification = c(1, 1),
          legend.spacing.y = unit(0.5, "pt"), plot.margin = margin(1, 9, 5, 5))

  if (log_y) {
    p <- p + scale_y_continuous(trans = pseudo_log_trans(sigma = 1, base = 10),
                                breaks = y_breaks, expand = expansion(mult = c(0, 0.05))) +
      coord_cartesian(xlim = c(-1000, 1000), clip = "off")
  } else {
    p <- p + coord_cartesian(ylim = c(0, ymax), xlim = c(-1000, 1000), expand = FALSE, clip = "off")
  }
  p
}

per_focal_totals <- function(file, id_col, group_col, value_col = "Accumulated_dSNP_Count") {
  d <- read_tsv(file)
  d[[value_col]] <- as.numeric(d[[value_col]])
  agg <- aggregate(d[[value_col]], by = list(ID = d[[id_col]], Group = d[[group_col]]), FUN = sum)
  names(agg)[3] <- "Total"
  agg
}

manifest <- list()

# ===== main Fig. f: observed dSV carrier vs matched random background =====
obs_bg <- read_tsv("fig_f_observed_vs_matched_background_summary.tsv")
obs_lv  <- c("Matched_random_background", "Observed_dSV_carrier")   # observed drawn on top
obs_col <- c(Matched_random_background = "#3C6DA8", Observed_dSV_carrier = "#E64B35")
obs_lab <- c(Matched_random_background = "Matched background", Observed_dSV_carrier = "Observed dSV carrier")
obs_ylab <- c(Matched_random_background = "Background", Observed_dSV_carrier = "Observed")

obs_tot <- per_focal_totals("fig_f_observed_vs_matched_background_by_focal_dsv.tsv", "Focal_DSV_ID", "Series")
p_top <- make_raincloud(obs_tot, obs_lv, obs_col, obs_ylab,
                      "Total accumulated dSNP per focal dSV (±1 Mb)", paired = TRUE, tag = "f")
p_bot <- make_profile(obs_bg, "Accumulated_dSNP_Count", "Bootstrap_CI_Low", "Bootstrap_CI_High",
                      rev(obs_lv), obs_col, obs_lab, "Accumulated dSNP count",
                      log_y = TRUE, y_breaks = c(0, 25, 50, 100, 250, 500, 1000))
p_obs <- p_top / p_bot + plot_layout(heights = c(0.42, 1))
manifest[[length(manifest) + 1]] <- save_fig(p_obs, "dsv_observed_v2", 200, 148)

# ===== companion: DEL vs INS in the same style (linear y; range is small) =====
svt <- read_tsv("fig_f_svtype_class_summary.tsv")
svt_n <- tapply(svt$Effective_Focal_DSV_Count, svt$Series, max)
svt_lv  <- c("INS", "DEL")   # DEL drawn on top
svt_col <- c(INS = "#E64B35", DEL = "#3C6DA8")
svt_lab <- c(INS = sprintf("INS (n=%s)", svt_n[["INS"]]), DEL = sprintf("DEL (n=%s)", svt_n[["DEL"]]))
svt_ylab <- c(INS = "INS", DEL = "DEL")

svt_tot <- per_focal_totals("fig_f_svtype_class_by_focal_dsv.tsv", "Focal_DSV_ID", "SVTYPE")
svt_tot <- svt_tot[svt_tot$Group %in% c("DEL", "INS"), , drop = FALSE]
p_top2 <- make_raincloud(svt_tot, svt_lv, svt_col, svt_ylab,
                       "Total accumulated dSNP per focal dSV (±1 Mb)", paired = FALSE, tag = "f")
p_bot2 <- make_profile(svt[svt$Series %in% c("DEL", "INS"), ], "Mean_Per_Focal_DSV", "Mean_CI_Low", "Mean_CI_High",
                       rev(svt_lv), svt_col, svt_lab, "Mean accumulated dSNPs per focal dSV",
                       log_y = FALSE)
p_delins <- p_top2 / p_bot2 + plot_layout(heights = c(0.42, 1))
manifest[[length(manifest) + 1]] <- save_fig(p_delins, "dsv_del_ins_v2", 200, 148)

manifest_df <- do.call(rbind, manifest)
write.table(manifest_df, file.path(fig_dir, "figure_f_v2_render_manifest.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
cat("[OK] Wrote Fig. f v2 (1.png style, connected 2-panel, log y) to ", fig_dir, "\n", sep = "")
