#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ape)
  library(ggtree)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(readr)
  library(ragg)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) stop("usage: render_phylogenies.R <flat_dir> <out_root>")
flat <- normalizePath(args[[1]], mustWork = TRUE)
out_root <- normalizePath(args[[2]], mustWork = FALSE)
out_dir <- file.path(out_root, "01_Extended_Data_Figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

genome_cols <- c(
  American_hap1 = "#009E73", Dura = "#56B4E9", Pisifera = "#0072B2",
  Cocos_nucifera = "#E69F00", Areca_catechu = "#D55E00",
  Phoenix_dactylifera = "#CC79A7", Nypa_fruticans = "#999999",
  Arabidopsis_thaliana = "#000000"
)

draw_tree <- function(prefix, title) {
  tree <- read.tree(file.path(flat, paste0("Fig2e_", prefix, "_phylogeny_ml.contree")))
  meta <- read_tsv(file.path(flat, paste0("Fig2e_", prefix, "_phylogeny_tip_metadata.tsv")),
                   show_col_types = FALSE)
  stopifnot(setequal(tree$tip.label, meta$Tip))
  meta <- meta %>% transmute(label = Tip,
                             display = paste0(Gene_locus, " | ", Genome),
                             Genome,
                             oil_palm_label = ifelse(Oil_palm, "bold", "plain"))
  p <- ggtree(tree, size = 0.28) %<+% meta +
    geom_tiplab(aes(label = display, colour = Genome, fontface = oil_palm_label),
                size = 1.35, align = TRUE, linesize = 0.18, offset = 0.01) +
    scale_colour_manual(values = genome_cols, na.value = "grey30") +
    xlim_tree(max(ggtree(tree)$data$x) * 2.15) +
    labs(title = title, colour = "Genome") +
    theme_tree2() +
    theme(plot.title = element_text(face = "bold", size = 9),
          legend.position = "none", plot.margin = margin(3, 3, 3, 3))
  node_data <- p$data %>% filter(!isTip) %>%
    mutate(support = suppressWarnings(as.numeric(label))) %>%
    filter(!is.na(support), support >= 70)
  if (nrow(node_data) > 0) {
    p <- p + geom_text(data = node_data, aes(x = x, y = y, label = round(support)),
                       inherit.aes = FALSE, size = 1.15, colour = "grey25", nudge_y = 0.18)
  }
  p
}

p1 <- draw_tree("FAD", "FAD2/6 and FAD3/7 maximum-likelihood tree")
p2 <- draw_tree("FAT", "FATA/FATB maximum-likelihood tree")
p3 <- draw_tree("DGAT1", "DGAT1 maximum-likelihood tree")
p4 <- draw_tree("DGAT2", "DGAT2 maximum-likelihood tree")
panel <- ((p1 | p2) / (p3 | p4)) +
  plot_annotation(
    title = "Figure 2 Extended Data 4 | Phylogenetic support for curated lipid-enzyme assignments",
    subtitle = "IQ-TREE consensus topologies; internal values shown when ultrafast-bootstrap support is at least 70. Oil-palm labels are bold.",
    tag_levels = "a"
  )

stem <- file.path(out_dir, "ED_Fig2_04_ML_phylogenies")
pdf(paste0(stem, ".pdf"), width = 16, height = 14, useDingbats = FALSE, family = "Helvetica")
print(panel); dev.off()
svg(paste0(stem, ".svg"), width = 16, height = 14, family = "sans")
print(panel); dev.off()
agg_png(paste0(stem, ".png"), width = 16, height = 14, units = "in", res = 450,
        background = "white")
print(panel); dev.off()

cat("Rendered four IQ-TREE consensus phylogenies\n")
