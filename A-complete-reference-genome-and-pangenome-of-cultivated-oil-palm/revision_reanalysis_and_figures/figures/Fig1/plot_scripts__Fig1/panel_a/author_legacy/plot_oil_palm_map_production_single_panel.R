#!/usr/bin/env Rscript

# Figure 1a: compact single-panel map with exact country production labels.
#
# The map displays the 2019 approximately 100-km remote-sensing detection
# cells. Cell colours encode the USDA PSD 2026 palm-oil production class.
# Exact production values for all countries in the USDA table are labelled on
# the same map, avoiding separate regional panels and large unused margins.

rm(list = ls())

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(maps)
  library(grid)
})

OUT_DIR <- "${ANALYSIS_DIR}/22_answer_reviews/11_map_zone"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

EVIDENCE_TSV <- file.path(OUT_DIR, "public_data_colored_tile_evidence.tsv")
USDA_TSV <- file.path(OUT_DIR, "usda_psd_oil_palm_production_2026.tsv")

FONT <- "sans"
COL_OCEAN <- "#F6F8FA"
COL_LAND <- "#E5E7E4"
COL_TEXT <- "#263238"
COL_LINE <- "#647178"
COL_DETECTED_ONLY <- "#B7C1C5"

class_levels <- c(
  "Main producer (>10 Mt)",
  "Secondary producer (0.5-10 Mt)",
  "Minor producer (<0.5 Mt)",
  "Detected plantations (remote sensing only)"
)

class_labels <- c(
  "Main producer (>10 Mt)" = "Main producer >10 Mt",
  "Secondary producer (0.5-10 Mt)" = "Secondary producer 0.5-10 Mt",
  "Minor producer (<0.5 Mt)" = "Minor producer <0.5 Mt",
  "Detected plantations (remote sensing only)" = "Detected cell; no USDA record"
)

class_cols <- c(
  "Main producer (>10 Mt)" = "#0072B2",
  "Secondary producer (0.5-10 Mt)" = "#009E73",
  "Minor producer (<0.5 Mt)" = "#E69F00",
  "Detected plantations (remote sensing only)" = COL_DETECTED_ONLY
)

evidence <- read_tsv(EVIDENCE_TSV, show_col_types = FALSE) %>%
  mutate(
    tile_fill = if_else(
      is.na(map_class),
      "Detected plantations (remote sensing only)",
      map_class
    ),
    tile_fill = factor(tile_fill, levels = class_levels)
  )

usda <- read_tsv(USDA_TSV, show_col_types = FALSE) %>%
  mutate(
    production_mt_million = production_mt / 1e6,
    production_class = factor(production_class, levels = class_levels[1:3])
  )

stopifnot(nrow(evidence) == 634L)
stopifnot(nrow(usda) == 28L)

world <- map_data("world") %>%
  filter(region != "Antarctica")

fmt_mt <- function(x) {
  ifelse(
    x >= 10,
    sprintf("%.1f Mt", x),
    ifelse(x >= 1, sprintf("%.2f Mt", x), sprintf("%.3f Mt", x))
  )
}

# The label positions are arranged around the three main oil-palm belts. They
# are fixed so the figure remains stable and reproducible between runs.
label_df <- tribble(
  ~country_name, ~short_name, ~x, ~y, ~lx, ~ly, ~hjust,
  "Mexico", "Mexico", -94.0, 17.5, -107.0, 22.8, 0,
  "Guatemala", "Guatemala", -90.7, 15.2, -106.0, 17.0, 0,
  "Honduras", "Honduras", -86.2, 15.0, -106.0, 11.2, 0,
  "Costa Rica", "Costa Rica", -84.0, 9.8, -106.0, 6.2, 0,
  "Ecuador", "Ecuador", -78.5, -1.0, -103.5, -2.5, 0,
  "Peru", "Peru", -76.0, -7.0, -101.5, -9.0, 0,
  "Colombia", "Colombia", -73.5, 5.5, -64.0, 14.2, 0.5,
  "Venezuela", "Venezuela", -66.0, 8.0, -59.0, 6.0, 0.5,
  "Dominican Republic", "Dominican Rep.", -70.2, 19.0, -61.0, 23.0, 0.5,
  "Brazil", "Brazil", -48.0, -2.0, -40.5, 3.0, 1,
  "Senegal", "Senegal", -14.5, 14.0, -17.5, 23.0, 0,
  "Guinea", "Guinea", -10.5, 10.5, -7.5, 23.0, 0.5,
  "Sierra Leone", "Sierra Leone", -11.5, 8.2, -28.0, 15.8, 0,
  "Liberia", "Liberia", -9.5, 6.4, -28.0, 6.0, 0,
  "Cote d'Ivoire", "Cote d'Ivoire", -5.5, 6.5, -18.0, 11.2, 0,
  "Ghana", "Ghana", -1.2, 6.5, -2.0, 16.5, 0.5,
  "Togo", "Togo", 1.2, 7.0, 5.0, 12.5, 0.5,
  "Benin", "Benin", 2.3, 9.0, 8.5, 17.0, 0.5,
  "Nigeria", "Nigeria", 7.5, 8.0, 14.0, 14.0, 0.5,
  "Cameroon", "Cameroon", 11.0, 4.5, 18.0, 8.5, 0.5,
  "Congo (Kinshasa)", "DRC", 24.0, -1.5, 31.0, 2.8, 1,
  "Angola", "Angola", 17.5, -9.0, 21.0, -14.8, 0.5,
  "India", "India", 78.0, 11.0, 71.5, 22.0, 0,
  "Thailand", "Thailand", 100.0, 9.0, 96.0, 22.0, 0.5,
  "Malaysia", "Malaysia", 113.5, 4.0, 112.0, 14.0, 0.5,
  "Indonesia", "Indonesia", 117.0, -2.0, 106.0, -14.8, 0,
  "Philippines", "Philippines", 122.0, 11.0, 139.0, 19.0, 1,
  "Papua New Guinea", "Papua New Guinea", 146.0, -6.0, 151.0, -14.8, 1
) %>%
  left_join(
    usda %>% select(country_name_psd, production_mt_million, production_class),
    by = c("country_name" = "country_name_psd")
  ) %>%
  mutate(label = paste0(short_name, "\n", fmt_mt(production_mt_million)))

stopifnot(identical(sort(label_df$country_name), sort(usda$country_name_psd)))

# The limits focus on the oil-palm belt and remove the high-latitude empty
# margins that dominated the previous global-plus-insets composition.
MAP_XLIM <- c(-110, 165)
MAP_YLIM <- c(-22, 26)

cells <- evidence %>%
  filter(
    xmax >= MAP_XLIM[1], xmin <= MAP_XLIM[2],
    ymax >= MAP_YLIM[1], ymin <= MAP_YLIM[2]
  )

legend_df <- tibble(
  x = seq_along(class_levels),
  y = 1,
  tile_fill = factor(class_levels, levels = class_levels)
)

map_plot <- ggplot() +
  geom_polygon(
    data = world,
    aes(x = long, y = lat, group = group),
    fill = COL_LAND, color = NA
  ) +
  geom_rect(
    data = cells,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = tile_fill),
    color = "white", linewidth = 0.06, alpha = 0.98
  ) +
  geom_hline(
    yintercept = c(-20, 0, 20),
    color = c("#D3D8DA", "#BDC5C8", "#D3D8DA"),
    linetype = c("dotted", "solid", "dotted"),
    linewidth = c(0.18, 0.18, 0.18)
  ) +
  geom_segment(
    data = label_df,
    aes(x = x, y = y, xend = lx, yend = ly),
    color = COL_LINE, linewidth = 0.20
  ) +
  geom_point(
    data = label_df,
    aes(x = x, y = y),
    shape = 21, size = 1.25, stroke = 0.22,
    fill = "#59666D", color = "white"
  ) +
  geom_label(
    data = label_df,
    aes(x = lx, y = ly, label = label, hjust = hjust),
    family = FONT, size = 2.05, lineheight = 0.88,
    color = COL_TEXT, fill = "white", alpha = 0.92,
    label.padding = unit(0.10, "lines"), linewidth = 0.12
  ) +
  annotate(
    "text", x = 161.5, y = 24.7, label = "N", size = 3.1,
    family = FONT, color = COL_LINE, fontface = "bold"
  ) +
  annotate(
    "segment", x = 161.5, xend = 161.5, y = 21.7, yend = 23.8,
    linewidth = 0.4, color = COL_LINE,
    arrow = arrow(length = unit(0.11, "cm"))
  ) +
  annotate(
    "segment", x = -104, xend = -54, y = -20.2, yend = -20.2,
    linewidth = 0.5, color = COL_LINE
  ) +
  annotate(
    "segment", x = -104, xend = -104, y = -20.65, yend = -19.75,
    linewidth = 0.5, color = COL_LINE
  ) +
  annotate(
    "segment", x = -54, xend = -54, y = -20.65, yend = -19.75,
    linewidth = 0.5, color = COL_LINE
  ) +
  annotate(
    "text", x = -79, y = -21.2, label = "5,000 km", size = 2.4,
    family = FONT, color = COL_LINE
  ) +
  geom_point(
    data = legend_df,
    aes(x = x, y = y, fill = tile_fill),
    shape = 22, size = 5.4, alpha = 0, show.legend = TRUE
  ) +
  scale_fill_manual(
    values = class_cols,
    breaks = class_levels,
    labels = unname(class_labels[class_levels]),
    drop = FALSE
  ) +
  guides(fill = guide_legend(
    nrow = 1, byrow = TRUE, title = NULL,
    override.aes = list(size = 5.5, alpha = 1)
  )) +
  coord_quickmap(xlim = MAP_XLIM, ylim = MAP_YLIM, expand = FALSE) +
  labs(
    title = "a  Global oil-palm distribution and country-level production",
    subtitle = "Compact view of the main oil-palm belt; all country production values are shown in million tonnes (Mt)"
  ) +
  theme_void(base_family = FONT) +
  theme(
    panel.background = element_rect(fill = COL_OCEAN, color = "#AEB7BC", linewidth = 0.3),
    plot.background = element_rect(fill = "white", color = NA),
    plot.title = element_text(
      family = FONT, face = "bold", size = 14, color = COL_TEXT,
      margin = margin(b = 2)
    ),
    plot.subtitle = element_text(
      family = FONT, size = 8.5, color = "#5C676D",
      margin = margin(b = 3)
    ),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.text = element_text(size = 8.6, color = COL_TEXT, lineheight = 0.88),
    legend.key.width = unit(0.62, "cm"),
    legend.key.height = unit(0.42, "cm"),
    legend.spacing.x = unit(0.24, "cm"),
    legend.margin = margin(t = 2, b = 0),
    plot.margin = margin(3, 4, 2, 4)
  )

png_path <- file.path(OUT_DIR, "Figure1a_oil_palm_production_single_panel.png")
pdf_path <- file.path(OUT_DIR, "Figure1a_oil_palm_production_single_panel.pdf")

pdf_tmp <- tempfile(fileext = ".pdf")
ggsave(pdf_tmp, map_plot, width = 16, height = 4.9, dpi = 600,
       device = cairo_pdf, bg = "white")
stopifnot(file.copy(pdf_tmp, pdf_path, overwrite = TRUE))

png_tmp_base <- tempfile()
pdftoppm_status <- system2(
  "pdftoppm",
  args = c("-png", "-r", "300", "-singlefile", pdf_tmp, png_tmp_base)
)
stopifnot(pdftoppm_status == 0L)
stopifnot(file.copy(paste0(png_tmp_base, ".png"), png_path, overwrite = TRUE))

label_export <- label_df %>%
  select(country_name, short_name, production_mt_million, production_class, label) %>%
  mutate(label = gsub("\\n", " | ", label, fixed = FALSE))
write_tsv(label_export, file.path(OUT_DIR, "Figure1a_production_labels_2026_single_panel.tsv"))

cat("Wrote:\n", png_path, "\n", pdf_path, "\n", sep = "")
cat("Exact production labels included: ", nrow(label_export), " countries\n", sep = "")
