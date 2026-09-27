#!/usr/bin/env Rscript
# Figure 1a (point version): oil-palm cultivation evidence from the Descals et al.
# 2019 global oil-palm grid (634 detected 100-km cells), plotted as points.
# - No country polygons, no border strokes (previous land-shape version distorted
#   the true cultivation extent because country polygons imply uniform planting).
# - Land is shown as a neutral light-grey silhouette for geographic context only.
# - Point colour: USDA PSD Online Oil, Palm production class (market year 2026).

rm(list = ls())
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(maps)
})

SRC_DIR <- "${ANALYSIS_DIR}/22_answer_reviews/11_map_zone"
OUT_DIR <- "${ANALYSIS_DIR}/22_answer_reviews/00_ms/05_MS/0918_revision"
dir.create(OUT_DIR, showWarnings = FALSE)

EVIDENCE_TSV <- file.path(SRC_DIR, "public_data_colored_tile_evidence.tsv")

COL_MAIN <- "#5365AF"
COL_SECONDARY <- "#8390C9"
COL_MINOR <- "#B99BD6"
COL_LAND <- "#EDEDED"
COL_OCEAN <- "#F7F9FC"
COL_TEXT <- "grey22"

TILES_TSV <- file.path(SRC_DIR, "descals_grid_withOP_tiles.tsv")

pt_tab <- read_tsv(TILES_TSV, show_col_types = FALSE) %>%
  select(long = centroid_long, lat = centroid_lat)

COL_POINT <- "#1F5F8F"

world <- map_data("world") %>%
  mutate(region = if_else(region == "Taiwan", "China", region))

tropic_band <- tibble(xmin = -180, xmax = 180, ymin = -20, ymax = 20)

p <- ggplot() +
  geom_rect(
    data = tropic_band,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    fill = "#EEF5F8",
    alpha = 0.55
  ) +
  geom_polygon(
    data = world,
    aes(x = long, y = lat, group = group),
    fill = COL_LAND,
    color = NA
  ) +
  geom_point(
    data = pt_tab,
    aes(x = long, y = lat),
    shape = 16,
    size = 3.0,
    stroke = 0,
    color = COL_POINT
  ) +
  scale_x_continuous(breaks = seq(-120, 160, 40), expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(-30, 30, 15), expand = c(0, 0)) +
  coord_quickmap(xlim = c(-128, 180), ylim = c(-37, 38)) +
  theme_void(base_family = "sans") +
  theme(
    panel.background = element_rect(fill = COL_OCEAN, color = "grey78", linewidth = 0.22),
    plot.background = element_rect(fill = "white", color = NA),
    axis.text = element_text(size = 10, color = "grey45"),
    axis.ticks = element_line(color = "grey75", linewidth = 0.2),
    plot.margin = margin(4, 5, 4, 5)
  )

ggsave(file.path(OUT_DIR, "Figure1a_oil_palm_cultivation_points.pdf"), p,
       width = 14, height = 4.9, dpi = 600, device = cairo_pdf, bg = "white")
ggsave(file.path(OUT_DIR, "Figure1a_oil_palm_cultivation_points.png"), p,
       width = 14, height = 4.9, dpi = 600, bg = "white")

cat("wrote point map to ", OUT_DIR, "\n", sep = "")
cat("points: ", nrow(pt_tab), "\n", sep = "")