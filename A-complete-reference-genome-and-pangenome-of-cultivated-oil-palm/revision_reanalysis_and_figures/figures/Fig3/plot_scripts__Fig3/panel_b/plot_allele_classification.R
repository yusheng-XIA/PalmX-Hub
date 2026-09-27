#!/usr/bin/env Rscript
# ============================================================
# 等位基因分类可视化 - T2T基因组
# 统一配色 (3个类别各一个颜色), 标签不压线
# ============================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(cowplot)
})

# ============================================================
# 3个类别统一配色 (从用户色板中选取)
# ============================================================
col_biallelic <- rgb(152, 176, 209, maxColorValue = 255)  # 蓝灰
col_samecds   <- rgb(195, 182, 230, maxColorValue = 255)  # 淡紫
col_hapspec   <- rgb(167, 217, 221, maxColorValue = 255)  # 青绿

# ============================================================
# 数据
# ============================================================
plot_data <- data.frame(
  haplotype = c(
    "American (Hap A)", "American (Hap A)", "American (Hap A)",
    "Africa (Hap B)",   "Africa (Hap B)",   "Africa (Hap B)",
    "BK (Hap 1)",       "BK (Hap 1)",       "BK (Hap 1)",
    "BK (Hap 2)",       "BK (Hap 2)",       "BK (Hap 2)"
  ),
  category = rep(c("Biallelic genes", "Allele with same CDS", "Haplotype-specific genes"), 4),
  count = c(27106,264,3214, 27106,264,7430, 18439,7747,7720, 18439,7747,8486),
  stringsAsFactors = FALSE
)

plot_data <- plot_data %>%
  group_by(haplotype) %>%
  mutate(total = sum(count), pct = count/total*100) %>%
  ungroup()

plot_data$haplotype <- factor(plot_data$haplotype,
  levels = rev(c("American (Hap A)","Africa (Hap B)","BK (Hap 1)","BK (Hap 2)")))

# 堆叠顺序: Biallelic在左(先堆), Same CDS中, Hap-specific在右
plot_data$category <- factor(plot_data$category,
  levels = c("Haplotype-specific genes","Allele with same CDS","Biallelic genes"))

cat_cols <- c(
  "Biallelic genes"           = col_biallelic,
  "Allele with same CDS"     = col_samecds,
  "Haplotype-specific genes"  = col_hapspec
)

# ============================================================
# 标签: 避开边框
# ============================================================
label_data <- plot_data %>%
  arrange(haplotype, desc(category)) %>%
  group_by(haplotype) %>%
  mutate(
    cum_pct   = cumsum(pct),
    seg_left  = cum_pct - pct,       # 段左边界
    seg_right = cum_pct,             # 段右边界
    label_pos = seg_left + pct/2,    # 段中心
    # 两行: 段宽 > 12%
    label_txt = ifelse(pct > 12,
      paste0(round(pct,1), "%\n(n=", format(count, big.mark=","), ")"), ""),
    # 一行: 段宽 5~12%
    label_small = ifelse(pct > 5 & pct <= 12,
      paste0(round(pct,1), "%\n(n=", format(count, big.mark=","), ")"), ""),
    # 很窄 1~5%: 只放百分比
    label_tiny = ifelse(pct > 1.5 & pct <= 5,
      paste0(round(pct,1), "%"), "")
  ) %>%
  ungroup()

# 对于小标签, 如果文字中心太靠近左右边框(< 4%), 向中间微调
label_data <- label_data %>%
  mutate(
    # 估算文字占的宽度百分比(两行标签约占8%, 一行约5%, tiny约3%)
    est_w = case_when(
      nchar(label_txt)   > 0 ~ 8,
      nchar(label_small) > 0 ~ 6,
      nchar(label_tiny)  > 0 ~ 3,
      TRUE ~ 0
    ),
    # 如果标签左沿会超出段左边界, 右移
    label_pos = ifelse(label_pos - est_w/2 < seg_left + 1.5,
                       seg_left + est_w/2 + 1.5, label_pos),
    # 如果标签右沿会超出段右边界, 左移
    label_pos = ifelse(label_pos + est_w/2 > seg_right - 1.5,
                       seg_right - est_w/2 - 1.5, label_pos)
  )

# ============================================================
# 主图
# ============================================================
p_main <- ggplot(plot_data, aes(x = pct, y = haplotype, fill = category)) +
  geom_bar(stat = "identity", width = 0.55, color = "black", linewidth = 0.9) +
  # 大标签
  geom_text(
    data = label_data %>% filter(nchar(label_txt) > 0),
    aes(x = label_pos, label = label_txt),
    family = "Arial", size = 4.5, color = "black", lineheight = 0.85
  ) +
  # 中标签
  geom_text(
    data = label_data %>% filter(nchar(label_small) > 0),
    aes(x = label_pos, label = label_small),
    family = "Arial", size = 4, color = "black", lineheight = 0.85
  ) +
  # 小标签
  geom_text(
    data = label_data %>% filter(nchar(label_tiny) > 0),
    aes(x = label_pos, label = label_tiny),
    family = "Arial", size = 3.5, color = "black"
  ) +
  scale_fill_manual(
    values = cat_cols, name = NULL,
    guide = guide_legend(reverse = TRUE, nrow = 1,
      override.aes = list(color = "black", linewidth = 0.5))
  ) +
  scale_x_continuous(expand = c(0,0), breaks = seq(0,100,25),
                     labels = paste0(seq(0,100,25))) +
  geom_hline(yintercept = 2.5, linetype = "dotted", color = "gray50", linewidth = 0.5) +
  # 右侧总数
  geom_text(
    data = plot_data %>% distinct(haplotype, total),
    aes(x = 103, y = haplotype, label = paste0("n=", format(total, big.mark=",")), fill = NULL),
    family = "Arial", size = 4.5, hjust = 0, color = "black"
  ) +
  # 分组标注
  annotate("text", x = 120, y = 3.5, label = "American\n\u00d7 Africa",
    family = "Arial", size = 4.5, hjust = 0.5, color = "black", lineheight = 0.85) +
  annotate("text", x = 120, y = 1.5, label = "BK\nhap1 \u00d7 hap2",
    family = "Arial", size = 4.5, hjust = 0.5, color = "black", lineheight = 0.85) +
  coord_cartesian(clip = "off", xlim = c(0,100)) +
  labs(x = "Proportion of genes (%)", y = NULL) +
  theme_classic(base_size = 16, base_family = "Arial") +
  theme(
    legend.position    = "bottom",
    legend.text        = element_text(size = 14, color = "black", family = "Arial"),
    legend.key.size    = unit(0.55, "cm"),
    legend.key.width   = unit(0.8, "cm"),
    legend.spacing.x   = unit(0.5, "cm"),
    legend.margin      = margin(t = 8),
    axis.text.y  = element_text(size = 16, color = "black"),
    axis.text.x  = element_text(size = 16, color = "black"),
    axis.title.x = element_text(size = 16, color = "black"),
    axis.line    = element_line(color = "black", linewidth = 0.6),
    axis.ticks   = element_line(color = "black", linewidth = 0.6),
    plot.margin  = margin(15, 80, 5, 10)
  )

# ============================================================
# 输出
# ============================================================
outdir <- "${ANALYSIS_DIR}/21_MS/02_result/01_figure"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

ggsave(file.path(outdir, "Fig_allele_classification_T2T.pdf"),
  p_main, width = 14, height = 6.5, dpi = 300, device = cairo_pdf)

cat("\nDone! Output:\n")
cat("  ", file.path(outdir, "Fig_allele_classification_T2T.pdf"), "\n")