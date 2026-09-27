# plot_figC_v2.R  ── 精致双圆环 + 分层柱状图
library(ggplot2)
library(dplyr)
library(patchwork)
library(scales)

# ── 数据导入 ──────────────────────────────────────
s <- read.csv("sv_summary_for_plot.csv", header = FALSE,
              row.names = 1, col.names = c("key", "value"))

total        <- as.numeric(s["total",      "value"])   # 14009
n_genic      <- as.numeric(s["genic",      "value"])   # 1328
n_intergenic <- as.numeric(s["intergenic", "value"])   # 12681
n_te         <- as.numeric(s["te_overlap", "value"])   # 9228
n_no_te      <- as.numeric(s["no_te",      "value"])   # 4781
n_exon       <- as.numeric(s["exon",       "value"])   # 36
n_intron     <- as.numeric(s["intron",     "value"])   # 692
n_up2k       <- as.numeric(s["up2k",       "value"])   # 166
n_down2k     <- as.numeric(s["down2k",     "value"])   # 95

# ── 颜色定义 ──────────────────────────────────
col_genic      <- "#E57373"   # 红：genic
col_intergenic <- "#64B5F6"   # 蓝：intergenic
col_te         <- "#4DB6AC"   # 青绿：TE overlap
col_no_te      <- "#FFF176"   # 黄：no TE

col_down2k <- "#EF5350"   # 红
col_exon   <- "#66BB6A"   # 绿
col_intron <- "#FFA726"   # 橙
col_up2k   <- "#42A5F5"   # 蓝

# ════════════════════════════════════════════════
# A. 双圆环图 (完美紧贴 + 统一图层)
# ════════════════════════════════════════════════
all_levels <- c(
  "SV overlapping with TE",
  "SV not overlapping with TE",
  "SV in genic region",
  "SV in intergenic region"
)

# 合并数据框，利用 x 轴位置区分内外环
df_combined <- bind_rows(
  data.frame(category = "SV in genic region",         count = n_genic,      ring = 3.5),
  data.frame(category = "SV in intergenic region",    count = n_intergenic, ring = 3.5),
  data.frame(category = "SV overlapping with TE",     count = n_te,         ring = 2.5),
  data.frame(category = "SV not overlapping with TE", count = n_no_te,      ring = 2.5)
) %>%
  mutate(category = factor(category, levels = all_levels))

ring_colors <- c(
  "SV in genic region"         = col_genic,
  "SV in intergenic region"    = col_intergenic,
  "SV overlapping with TE"     = col_te,
  "SV not overlapping with TE" = col_no_te
)

p_donut <- ggplot(df_combined, aes(x = ring, y = count, fill = category)) +
  # width = 1.0 时，内环 [2.0,3.0]，外环 [3.0,4.0]，两环实现0间隙紧贴
  geom_col(width = 1.0, color = "white", linewidth = 0.6) +
  coord_polar(theta = "y", start = 0) +
  # xlim(1.0, 4.05): 1.0 是圆心的空白洞大小，4.05 控制外部边缘留白
  xlim(1.0, 4.05) +
  scale_fill_manual(
    values = ring_colors,
    name   = "Type",
    breaks = all_levels,
    guide  = guide_legend(
      override.aes = list(
        color     = c("grey55", "grey55", "white", "white"),
        linewidth = c(0.35, 0.35, 0, 0)
      )
    )
  ) +
  # ── 中心数字与文字（利用 vjust 垂直分开，永不重叠）──
  annotate("text", x = 1.0, y = 0, label = format(total, big.mark = ","),
           size = 6.5, fontface = "bold", color = "grey15", vjust = -0.2) +
  annotate("text", x = 1.0, y = 0, label = "SVs",
           size = 4.5, color = "grey40", vjust = 1.5) +
  theme_void(base_size = 11) +
  theme(
    legend.position  = "bottom",
    legend.direction = "vertical",
    legend.text      = element_text(size = 10, color = "grey20"),
    legend.key.size  = unit(0.45, "cm"),
    legend.key       = element_rect(color = NA, fill = NA),
    legend.title     = element_text(size = 11, face = "bold", color = "grey10"),
    legend.margin    = margin(t = 0, b = 2),
    plot.margin      = margin(5, 5, 5, 0)
  )

# ════════════════════════════════════════════════
# B. 分层柱状图 (标签置于右侧，颜色完美对应)
# ════════════════════════════════════════════════
sub_total <- n_down2k + n_exon + n_intron + n_up2k  # 989

bar_df <- data.frame(
  region = factor(
    c("Down 2k", "Exon", "Intron", "Up 2k"),
    levels = c("Down 2k", "Exon", "Intron", "Up 2k")
  ),
  count = c(n_down2k, n_exon, n_intron, n_up2k)
) %>%
  mutate(
    pct      = count / sub_total * 100,
    cum_pct  = cumsum(pct),
    midpoint = cum_pct - pct / 2   # 计算 Y 轴堆叠后的中心点
  )

bar_colors <- c(
  "Down 2k" = col_down2k,
  "Exon"    = col_exon,
  "Intron"  = col_intron,
  "Up 2k"   = col_up2k
)

p_bar <- ggplot(bar_df, aes(x = 1, y = pct, fill = region)) +
  # 柱子本身
  geom_col(width = 0.5, color = "white", linewidth = 0.5) +
  # 外部包裹虚线框稍微放宽一点包裹住柱子
  annotate("rect", xmin = 0.72, xmax = 1.28, ymin = 0, ymax = 100,
           fill = NA, color = "grey60", linewidth = 0.5, linetype = "dashed") +
  # 将文字放在柱子的右侧 (x=1.38)，避免窄柱子文字挤压
  geom_text(
    aes(x = 1.38, y = midpoint, label = sprintf("%s (%.1f%%)", region, pct), color = region),
    hjust = 0, size = 3.5, fontface = "bold", show.legend = FALSE
  ) +
  scale_fill_manual(values = bar_colors) +
  scale_color_manual(values = bar_colors) + # 文字颜色与段落色块一致
  # X 轴稍微延伸以容纳右侧长标签文字
  scale_x_continuous(limits = c(0.6, 3.5), expand = c(0, 0)) +
  scale_y_continuous(
    limits = c(0, 100),
    breaks = c(0, 25, 50, 75, 100),
    labels = paste0(c(0, 25, 50, 75, 100), "%"), 
    expand = c(0.01, 0.01) # 上下留极小的白边，避免边框被切
  ) +
  theme_classic(base_size = 10) +
  theme(
    axis.title       = element_blank(),
    axis.text.x      = element_blank(),
    axis.ticks.x     = element_blank(),
    axis.line.x      = element_blank(),
    axis.line.y      = element_line(color = "grey55", linewidth = 0.4),
    axis.ticks.y     = element_line(color = "grey55", linewidth = 0.4),
    axis.text.y      = element_text(size = 9, color = "grey30", margin = margin(r = 5)),
    legend.position  = "none",
    plot.margin      = margin(15, 10, 15, 0)
  )

# ════════════════════════════════════════════════
# C. 组合 & 输出 (调整宽度比例使排版更加紧凑)
# ════════════════════════════════════════════════
# 使用 2.5 : 1 的宽度比例，保证左侧圆环占比且右侧文字不被切除
final <- p_donut + p_bar +
  plot_layout(widths = c(2.5, 1)) &
  theme(plot.background = element_rect(fill = "white", color = NA))

ggsave("figC_v3.pdf", final, width = 9.5, height = 6.5, dpi = 300)
ggsave("figC_v3.png", final, width = 9.5, height = 6.5, dpi = 300)
cat("已保存: figC_v3.pdf / figC_v3.png\n")
