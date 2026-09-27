#!/usr/bin/env Rscript

# Figure 1a: detected oil-palm cells plus country-level palm-oil production.
#
# The global panel shows the 2019 remote-sensing evidence cells. Cell colours
# encode the USDA PSD 2026 production class of the associated country. The
# three lower panels retain the same evidence cells and add exact production
# labels in Mt for all countries present in the USDA production table.

rm(list = ls())

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(maps)
  library(patchwork)
  library(grid)
})

OUT_DIR <- "${ANALYSIS_DIR}/22_answer_reviews/11_map_zone"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

EVIDENCE_TSV <- file.path(OUT_DIR, "public_data_colored_tile_evidence.tsv")
USDA_TSV <- file.path(OUT_DIR, "usda_psd_oil_palm_production_2026.tsv")

FONT <- "sans"
COL_OCEAN <- "#F7F9FB"
COL_LAND <- "#E7E8E5"
COL_GRID <- "#C8CDD1"
COL_TEXT <- "#263238"
COL_DETECTED_ONLY <- "#B7C1C5"

# A high-contrast, colour-blind-friendly production palette.
class_levels <- c(
  "Main producer (>10 Mt)",
  "Secondary producer (0.5-10 Mt)",
  "Minor producer (<0.5 Mt)",
  "Detected plantations (remote sensing only)"
)
class_labels <- c(
  "Main producer (>10 Mt)" = "Main producer\n>10 Mt",
  "Secondary producer (0.5-10 Mt)" = "Secondary producer\n0.5-10 Mt",
  "Minor producer (<0.5 Mt)" = "Minor producer\n<0.5 Mt",
  "Detected plantations (remote sensing only)" = "Detected cell\nno USDA production record"
)
class_cols <- c(
  "Main producer (>10 Mt)" = "#0072B2",
  "Secondary producer (0.5-10 Mt)" = "#009E73",
  "Minor producer (<0.5 Mt)" = "#E69F00",
  "Detected plantations (remote sensing only)" = COL_DETECTED_ONLY
)

evidence <- read_tsv(EVIDENCE_TSV, show_col_types = FALSE) %>%
  mutate(
    map_class = factor(map_class, levels = class_levels),
    tile_fill = if_else(
      is.na(map_class),
      "Detected plantations (remote sensing only)",
      as.character(map_class)
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

# Label anchors and label positions are deliberately fixed so the figure is
# stable between runs and remains readable at manuscript scale.
label_df <- tribble(
  ~panel, ~country_name, ~short_name, ~x, ~y, ~lx, ~ly,
  "Latin America", "Mexico", "Mexico", -94.0, 17.5, -103.0, 22.2,
  "Latin America", "Guatemala", "Guatemala", -90.7, 15.2, -103.0, 13.2,
  "Latin America", "Honduras", "Honduras", -86.2, 15.0, -101.5, 8.7,
  "Latin America", "Costa Rica", "Costa Rica", -84.0, 9.8, -102.0, 4.0,
  "Latin America", "Colombia", "Colombia", -73.5, 5.5, -68.0, 13.8,
  "Latin America", "Venezuela", "Venezuela", -66.0, 8.0, -63.0, 3.2,
  "Latin America", "Ecuador", "Ecuador", -78.5, -1.0, -100.0, -1.0,
  "Latin America", "Peru", "Peru", -76.0, -7.0, -99.0, -8.0,
  "Latin America", "Brazil", "Brazil", -48.0, -2.0, -42.5, 5.8,
  "Latin America", "Dominican Republic", "Dominican Rep.", -70.2, 19.0, -55.5, 22.0,
  "Africa", "Senegal", "Senegal", -14.5, 14.0, -8.0, 16.4,
  "Africa", "Guinea", "Guinea", -10.5, 10.5, -4.5, 16.4,
  "Africa", "Sierra Leone", "Sierra Leone", -11.5, 8.2, -19.0, 9.8,
  "Africa", "Liberia", "Liberia", -9.5, 6.4, -19.0, 4.0,
  "Africa", "Cote d'Ivoire", "Cote d'Ivoire", -5.5, 6.5, -13.0, 12.5,
  "Africa", "Ghana", "Ghana", -1.2, 6.5, -2.0, 12.5,
  "Africa", "Togo", "Togo", 1.2, 7.0, 5.5, 9.5,
  "Africa", "Benin", "Benin", 2.3, 9.0, 8.0, 13.0,
  "Africa", "Nigeria", "Nigeria", 7.5, 8.0, 14.5, 10.8,
  "Africa", "Cameroon", "Cameroon", 11.0, 4.5, 17.0, 5.5,
  "Africa", "Congo (Kinshasa)", "DRC", 24.0, -1.5, 30.0, 1.5,
  "Africa", "Angola", "Angola", 17.5, -9.0, 14.0, -12.0,
  "Asia-Pacific", "India", "India", 78.0, 11.0, 72.0, 21.8,
  "Asia-Pacific", "Thailand", "Thailand", 100.0, 9.0, 94.0, 21.8,
  "Asia-Pacific", "Malaysia", "Malaysia", 113.5, 4.0, 108.0, 14.0,
  "Asia-Pacific", "Indonesia", "Indonesia", 117.0, -2.0, 106.0, -12.0,
  "Asia-Pacific", "Philippines", "Philippines", 122.0, 11.0, 137.0, 19.0,
  "Asia-Pacific", "Papua New Guinea", "Papua New Guinea", 146.0, -6.0, 145.0, -12.0
) %>%
  left_join(
    usda %>% select(country_name_psd, production_mt_million, production_class),
    by = c("country_name" = "country_name_psd")
  ) %>%
  mutate(
    label = paste0(short_name, "\n", fmt_mt(production_mt_million)),
    production_class = factor(production_class, levels = class_levels[1:3]),
    hjust = case_when(
      (panel == "Latin America" & lx <= -90) ~ 0,
      (panel == "Africa" & lx <= -12) ~ 0,
      (panel == "Asia-Pacific" & lx >= 135) ~ 1,
      TRUE ~ 0.5
    )
  )

# Confirm that every USDA country is represented by a production label. The
# map uses the abbreviated display name only; the source table remains exact.
label_countries <- sort(unique(label_df$country_name))
usda_countries <- sort(usda$country_name_psd)
stopifnot(identical(label_countries, usda_countries))

scale_bar <- function(x, y, width, label) {
  list(
    annotate("segment", x = x, xend = x + width, y = y, yend = y,
             linewidth = 0.55, color = "#647078"),
    annotate("segment", x = x, xend = x, y = y - 0.45, yend = y + 0.45,
             linewidth = 0.55, color = "#647078"),
    annotate("segment", x = x + width, xend = x + width,
             y = y - 0.45, yend = y + 0.45,
             linewidth = 0.55, color = "#647078"),
    annotate("text", x = x + width / 2, y = y - 1.25,
             label = label, size = 2.7, color = "#59646B")
  )
}

base_theme <- theme_void(base_family = FONT) +
  theme(
    panel.background = element_rect(fill = COL_OCEAN, color = "#AEB7BC", linewidth = 0.3),
    plot.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(3, 4, 3, 4),
    legend.position = "none"
  )

plot_cells <- function(xlim, ylim, show_labels = FALSE, title = NULL,
                       subtitle = NULL, panel_name = NULL) {
  cells <- evidence %>%
    filter(
      xmax >= xlim[1], xmin <= xlim[2],
      ymax >= ylim[1], ymin <= ylim[2]
    )

  p <- ggplot() +
    geom_polygon(
      data = world,
      aes(x = long, y = lat, group = group),
      fill = COL_LAND, color = NA
    ) +
    geom_rect(
      data = cells,
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = tile_fill),
      color = "white", linewidth = 0.06, alpha = 0.96
    ) +
    geom_hline(
      yintercept = c(-20, 0, 20),
      color = c("#D2D7D9", "#BDC5C8", "#D2D7D9"),
      linetype = c("dotted", "solid", "dotted"),
      linewidth = c(0.18, 0.18, 0.18)
    ) +
    scale_fill_manual(
      values = class_cols,
      breaks = class_levels,
      labels = unname(class_labels[class_levels]),
      drop = FALSE
    ) +
    coord_quickmap(xlim = xlim, ylim = ylim, expand = FALSE) +
    labs(title = title, subtitle = subtitle) +
    base_theme +
    theme(
      plot.title = element_text(
        family = FONT, face = "bold", size = 13, color = COL_TEXT,
        margin = margin(b = 2)
      ),
      plot.subtitle = element_text(
        family = FONT, size = 8.5, color = "#5C676D",
        margin = margin(b = 2)
      )
    )

  if (show_labels) {
    labels <- label_df %>% filter(panel == panel_name)
    p <- p +
      geom_segment(
        data = labels,
        aes(x = x, y = y, xend = lx, yend = ly),
        color = "#6C777C", linewidth = 0.22
      ) +
      geom_point(
        data = labels,
        aes(x = x, y = y),
        shape = 21, size = 1.5, stroke = 0.25,
        fill = "#69757A", color = "white"
      ) +
      geom_label(
        data = labels,
        aes(x = lx, y = ly, label = label, hjust = hjust),
        family = FONT, size = 2.5, lineheight = 0.9,
        color = COL_TEXT, fill = "white", alpha = 0.92,
        label.padding = unit(0.12, "lines"),
        linewidth = 0.15
      )
  }

  p
}

global_map <- plot_cells(
  xlim = c(-130, 180), ylim = c(-38, 38),
  title = "a  Detected closed-canopy oil-palm cells and country-level production",
  subtitle = "Global overview; colour indicates the USDA PSD 2026 palm-oil production class",
  show_labels = FALSE
) +
  annotate(
    "label", x = -124, y = 32.0,
    label = "Each coloured cell is an approximately 100-km\nremote-sensing detection (Descals et al., 2019)",
    hjust = 0, vjust = 1, size = 3.0, lineheight = 0.92,
    family = FONT, color = COL_TEXT, fill = "white", alpha = 0.9,
    label.padding = unit(0.18, "lines"), linewidth = 0.15
  ) +
  annotate(
    "text", x = 171, y = 33.5, label = "N", size = 3.5,
    family = FONT, color = "#5D686D", fontface = "bold"
  ) +
  annotate(
    "segment", x = 171, xend = 171, y = 29.5, yend = 32.0,
    linewidth = 0.45, color = "#5D686D",
    arrow = arrow(length = unit(0.12, "cm"))
  ) +
  scale_bar(-123, -33.2, 50, "5,000 km")

latin_map <- plot_cells(
  xlim = c(-106, -35), ylim = c(-15, 25),
  title = "b  Latin America",
  subtitle = "Country production labels: Mt, USDA PSD 2026",
  show_labels = TRUE, panel_name = "Latin America"
) +
  scale_bar(-103, -12.5, 12, "1,000 km")

africa_map <- plot_cells(
  xlim = c(-21, 36), ylim = c(-15, 20),
  title = "c  Africa",
  subtitle = "Country production labels: Mt, USDA PSD 2026",
  show_labels = TRUE, panel_name = "Africa"
) +
  scale_bar(-18, -12.5, 10, "1,000 km")

asia_map <- plot_cells(
  xlim = c(68, 160), ylim = c(-15, 28),
  title = "d  Asia-Pacific",
  subtitle = "Country production labels: Mt, USDA PSD 2026",
  show_labels = TRUE, panel_name = "Asia-Pacific"
) +
  scale_bar(71, -12.5, 12, "1,000 km")

legend_plot <- ggplot() +
  geom_point(
    data = tibble(x = seq_along(class_levels), y = 1, tile_fill = factor(class_levels, levels = class_levels)),
    aes(x = x, y = y, fill = tile_fill), shape = 22, size = 5.7,
    color = "#69757A", stroke = 0.25, alpha = 0, show.legend = TRUE
  ) +
  scale_fill_manual(
    values = class_cols, breaks = class_levels,
    labels = unname(class_labels[class_levels]), drop = FALSE
  ) +
  guides(fill = guide_legend(
    nrow = 1, byrow = TRUE, title = NULL,
    override.aes = list(size = 6.5, alpha = 1)
  )) +
  theme_void(base_family = FONT) +
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.text = element_text(size = 9.5, color = COL_TEXT, lineheight = 0.9),
    legend.key.width = unit(0.72, "cm"),
    legend.key.height = unit(0.48, "cm"),
    legend.spacing.x = unit(0.32, "cm"),
    plot.margin = margin(0, 0, 0, 0)
  )

figure <- global_map /
  (latin_map | africa_map | asia_map) /
  legend_plot +
  plot_layout(
    heights = c(1.15, 1.0, 0.16),
    guides = "keep"
  )

png_path <- file.path(OUT_DIR, "Figure1a_oil_palm_production_panels.png")
pdf_path <- file.path(OUT_DIR, "Figure1a_oil_palm_production_panels.pdf")

pdf_tmp <- tempfile(fileext = ".pdf")
ggsave(pdf_tmp, figure, width = 16, height = 9.4, dpi = 600, device = cairo_pdf, bg = "white")
stopifnot(file.copy(pdf_tmp, pdf_path, overwrite = TRUE))

# The local R installation has a broken PNG graphics device. Render the
# manuscript PNG from the validated vector PDF instead, preserving text and
# line sharpness while avoiding that device-specific crash.
png_tmp_base <- tempfile()
pdftoppm_status <- system2(
  "pdftoppm",
  args = c("-png", "-r", "300", "-singlefile", pdf_tmp, png_tmp_base)
)
stopifnot(pdftoppm_status == 0L)
stopifnot(file.copy(paste0(png_tmp_base, ".png"), png_path, overwrite = TRUE))

label_export <- label_df %>%
  select(panel, country_name, short_name, production_mt_million, production_class, label) %>%
  mutate(label = gsub("\\n", " | ", label, fixed = FALSE))
write_tsv(label_export, file.path(OUT_DIR, "Figure1a_production_labels_2026.tsv"))

cat("Wrote:\n", png_path, "\n", pdf_path, "\n", sep = "")
cat("Exact production labels included: ", nrow(label_export), " countries\n", sep = "")
