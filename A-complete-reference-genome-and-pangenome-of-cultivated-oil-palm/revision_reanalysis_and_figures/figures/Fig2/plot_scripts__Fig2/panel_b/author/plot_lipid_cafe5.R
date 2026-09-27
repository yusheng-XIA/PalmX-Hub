#!/usr/bin/env Rscript
#=====================================================
# 油脂基因家族扩张收缩可视化
# Fig: 气泡图 + 全局vs油脂对比 + 通路热图
#=====================================================

for (pkg in c("ggplot2", "dplyr", "tidyr", "RColorBrewer", "cowplot", "ggrepel"))
    if (!requireNamespace(pkg, quietly=TRUE))
        install.packages(pkg, repos="https://cloud.r-project.org")

library(ggplot2)
library(dplyr)
library(tidyr)
library(cowplot)

OUT <- "${ANALYSIS_DIR}/20_results/Figure2/07_new_figure/02_comparative_orthofinder/phylo_divtime_cafe/04_figures"
dir.create(OUT, showWarnings=FALSE, recursive=TRUE)

# ============================================================
# 数据准备
# ============================================================

# 油脂相关扩张基因 (三层汇总)
lipid_expand <- data.frame(
    OG = c("OG0012737","OG0014386","OG0002452","OG0006523","OG0001300",
           "OG0014259","OG0013166","OG0019450","OG0000040","OG0000409",
           "OG0022011","OG0013873","OG0016696","OG0016322","OG0007476",
           "OG0003800","OG0005194","OG0000011","OG0001323","OG0005223",
           "OG0004640","OG0006259","OG0013137","OG0017678","OG0000635",
           "OG0004665","OG0003737","OG0005835","OG0010136","OG0006144",
           "OG0017646","OG0000292"),
    Gene = c("WRI1","DGAT","Desaturase-1","Desaturase-2","Oleosin","Caleosin",
             "Enoyl reductase-1","Enoyl reductase-2","Thioesterase-1","Thioesterase-2",
             "GPAT","KCS11","CER3","Desaturase-3","Desaturase-4",
             "Desaturase-5","FAD-1","FAD-2","FAD-3","FAD-4",
             "FA metabolism-1","FA metabolism-2","FA metabolism-3","FA metabolism-4",
             "FA related-1","FA related-2","OLE18","TAG related-1","TAG related-2","ACS7",
             "MFP1","TAG related-3"),
    Pathway = c("Transcription factor","TAG synthesis","FA modification","FA modification",
                "Oil body","Oil body","FA synthesis","FA synthesis",
                "FA termination","FA termination","Kennedy pathway","FA elongation",
                "Wax/cutin","FA modification","FA modification",
                "FA modification","FA modification","FA modification","FA modification","FA modification",
                "FA metabolism","FA metabolism","FA metabolism","FA metabolism",
                "FA related","FA related","Oil body","TAG synthesis","TAG synthesis","TAG synthesis",
                "beta-oxidation","TAG synthesis"),
    Layer1 = c(5,5,7,6,2,0,3,2,1,1,1,1,1,1,2,4,0,3,1,2,3,2,3,1,1,1,1,1,1,1,3,2),
    Layer2 = c(0,0,5,5,2,1,2,0,1,0,0,1,1,1,2,0,0,1,1,0,3,2,0,0,1,1,1,1,1,0,3,0),
    Layer3 = c(3,2,3,3,1,0,1,1,0,0,1,1,0,0,0,0,1,0,0,1,1,1,1,1,0,1,0,0,0,1,2,0),
    stringsAsFactors = FALSE
)

# 油脂相关收缩基因 (主要的)
lipid_contract <- data.frame(
    OG = c("OG0001638","OG0000913","OG0001496","OG0001325"),
    Gene = c("SAD","Lipase family","BASS2/FAD","FA beta-oxidation"),
    Pathway = c("FA modification","Lipid degradation","FA modification","beta-oxidation"),
    Layer1 = c(0,-15,0,-2),
    Layer2 = c(0,0,-2,-2),
    Layer3 = c(-3,-3,0,0),
    stringsAsFactors = FALSE
)

# ============================================================
# Figure 1: 全局收缩 vs 油脂扩张 对比图
# ============================================================
cat("绘制Figure 1: 全局vs油脂对比...\n")

contrast_data <- data.frame(
    Layer = rep(c("Layer 1\n(30 species)", "Layer 2\n(Arecaceae)", "Layer 3\n(Oil crops)"), each=4),
    Category = rep(c("All families\n(expand)", "All families\n(contract)", 
                     "Lipid families\n(expand)", "Lipid families\n(contract)"), 3),
    Count = c(770, -2373, 114, -168,   # Layer1
              478, -1780, 75, -157,     # Layer2
              300, -2695, 52, -217),    # Layer3
    stringsAsFactors = FALSE
)

contrast_data$Layer <- factor(contrast_data$Layer, 
    levels=c("Layer 1\n(30 species)", "Layer 2\n(Arecaceae)", "Layer 3\n(Oil crops)"))
contrast_data$Category <- factor(contrast_data$Category,
    levels=c("All families\n(expand)", "All families\n(contract)", 
             "Lipid families\n(expand)", "Lipid families\n(contract)"))

# 计算百分比
pct_data <- data.frame(
    Layer = c("Layer 1\n(30 species)", "Layer 2\n(Arecaceae)", "Layer 3\n(Oil crops)"),
    pct_expand = c(114/770*100, 75/478*100, 52/300*100),
    pct_contract = c(168/2373*100, 157/1780*100, 217/2695*100),
    stringsAsFactors = FALSE
)

p1 <- ggplot(contrast_data, aes(x=Layer, y=Count, fill=Category)) +
    geom_col(position=position_dodge(width=0.8), width=0.7) +
    scale_fill_manual(values=c("All families\n(expand)"="#FFCDD2",
                                "All families\n(contract)"="#BBDEFB",
                                "Lipid families\n(expand)"="#D32F2F",
                                "Lipid families\n(contract)"="#1565C0")) +
    geom_hline(yintercept=0, color="grey30", size=0.5) +
    labs(x="", y="Number of gene families",
         title="Oil palm ancestor: genome-wide contraction vs lipid gene expansion",
         fill="") +
    theme_minimal(base_size=12) +
    theme(
        plot.title = element_text(size=13, face="bold"),
        legend.position = "bottom",
        panel.grid.major.x = element_blank(),
        axis.text.x = element_text(size=11)
    ) +
    # 添加比例注释
    annotate("text", x=1, y=900, label="14.8% lipid", size=3, color="#D32F2F", fontface="bold") +
    annotate("text", x=2, y=600, label="15.7% lipid", size=3, color="#D32F2F", fontface="bold") +
    annotate("text", x=3, y=450, label="17.3% lipid", size=3, color="#D32F2F", fontface="bold")

ggsave(file.path(OUT, "Fig_lipid_contrast.pdf"), p1, width=10, height=6, dpi=300)
cat("  ✓ Fig_lipid_contrast.pdf\n")


# ============================================================
# Figure 2: 油脂通路基因热图 (核心基因 × 3层)
# ============================================================
cat("绘制Figure 2: 油脂核心基因热图...\n")

# 选核心基因
core_genes <- data.frame(
    Gene = c("WRI1","DGAT","GPAT","Oleosin","Caleosin","OLE18",
             "KCS11","CER3","Thioesterase",
             "Desaturase (×5)","Enoyl reductase (×2)","FAD (×4)",
             "Lipase","ACS7","MFP1",
             "SAD"),
    Pathway = c("Transcription\nfactor","TAG\nsynthesis","Kennedy\npathway",
                "Oil body","Oil body","Oil body",
                "FA elongation","Wax/cutin","FA termination",
                "FA modification","FA synthesis","FA modification",
                "Lipid\ndegradation","TAG\nsynthesis","beta-oxidation",
                "FA modification"),
    Category = c(rep("Regulation", 1), rep("TAG/Oil body", 5), 
                 rep("FA synthesis/\nmodification", 6), rep("Degradation/\nOther", 4)),
    Layer1 = c(5, 5, 1, 2, 0, 1, 1, 1, 2, 20, 5, 6, -15, 1, 3, 0),
    Layer2 = c(0, 0, 0, 2, 1, 1, 1, 1, 1, 11, 2, 2, 0, 0, 3, 0),
    Layer3 = c(3, 2, 1, 1, 0, 0, 1, 0, 0, 6, 2, 2, -3, 1, 2, -3),
    stringsAsFactors = FALSE
)

core_long <- core_genes %>%
    pivot_longer(cols=c(Layer1, Layer2, Layer3), names_to="Layer", values_to="Change") %>%
    mutate(Layer = recode(Layer, "Layer1"="Global\n(30 spp)", 
                                 "Layer2"="Palm\n(12 spp)", 
                                 "Layer3"="Oil crops\n(16 spp)"))

# 基因排序
core_long$Gene <- factor(core_long$Gene, levels=rev(core_genes$Gene))
core_long$Layer <- factor(core_long$Layer, levels=c("Global\n(30 spp)", "Palm\n(12 spp)", "Oil crops\n(16 spp)"))

# 颜色: 红=扩张, 蓝=收缩, 白=0
p2 <- ggplot(core_long, aes(x=Layer, y=Gene, fill=Change)) +
    geom_tile(color="white", size=1) +
    geom_text(aes(label=ifelse(Change > 0, paste0("+", Change),
                               ifelse(Change < 0, as.character(Change), ""))),
              size=3.5, fontface="bold",
              color=ifelse(core_long$Change > 5 | core_long$Change < -5, "white", "grey20")) +
    scale_fill_gradient2(
        low="#1565C0", mid="white", high="#D32F2F",
        midpoint=0, limits=c(-15, 20),
        name="Gene family\nsize change"
    ) +
    # 分类标注 (左侧)
    facet_grid(Category ~ ., scales="free_y", space="free_y", switch="y") +
    labs(x="", y="",
         title="Lipid-related gene family changes on oil palm ancestral branch") +
    theme_minimal(base_size=12) +
    theme(
        plot.title = element_text(size=13, face="bold"),
        strip.placement = "outside",
        strip.text.y.left = element_text(angle=0, hjust=1, size=9, face="bold"),
        panel.grid = element_blank(),
        axis.text.y = element_text(size=10),
        axis.text.x = element_text(size=10),
        legend.position = "right"
    )

ggsave(file.path(OUT, "Fig_lipid_heatmap.pdf"), p2, width=8, height=8, dpi=300)
cat("  ✓ Fig_lipid_heatmap.pdf\n")


# ============================================================
# Figure 3: 气泡图 - 所有油脂扩张家族
# ============================================================
cat("绘制Figure 3: 气泡图...\n")

# 选择三层至少有一层 change >= 2 的基因
bubble_data <- lipid_expand %>%
    filter(pmax(Layer1, Layer2, Layer3) >= 2) %>%
    pivot_longer(cols=c(Layer1, Layer2, Layer3), names_to="Layer", values_to="Change") %>%
    filter(Change > 0) %>%
    mutate(Layer = recode(Layer, "Layer1"="Global (30 spp)", 
                                 "Layer2"="Palm (12 spp)", 
                                 "Layer3"="Oil crops (16 spp)"))

bubble_data$Gene <- factor(bubble_data$Gene, 
    levels=rev(unique(lipid_expand$Gene[order(-lipid_expand$Layer1)])))
bubble_data$Layer <- factor(bubble_data$Layer, 
    levels=c("Global (30 spp)", "Palm (12 spp)", "Oil crops (16 spp)"))

# 通路颜色
pathway_colors <- c(
    "Transcription factor"="#E91E63",
    "TAG synthesis"="#F44336",
    "Oil body"="#FF9800",
    "FA modification"="#4CAF50",
    "FA synthesis"="#2196F3",
    "FA termination"="#9C27B0",
    "FA elongation"="#00BCD4",
    "Kennedy pathway"="#FF5722",
    "Wax/cutin"="#795548",
    "FA metabolism"="#607D8B",
    "FA related"="#9E9E9E",
    "beta-oxidation"="#3F51B5",
    "Lipid degradation"="#455A64"
)

p3 <- ggplot(bubble_data, aes(x=Layer, y=Gene)) +
    geom_point(aes(size=Change, fill=Pathway), shape=21, color="grey30", stroke=0.3, alpha=0.85) +
    scale_size_continuous(range=c(2, 12), breaks=c(1,3,5,7,10), name="Expansion\n(# families)") +
    scale_fill_manual(values=pathway_colors, name="Pathway") +
    labs(x="", y="",
         title="Lipid gene family expansion on oil palm ancestral branch") +
    theme_minimal(base_size=11) +
    theme(
        plot.title = element_text(size=13, face="bold"),
        panel.grid.major = element_line(color="grey90"),
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(size=9),
        axis.text.x = element_text(size=10),
        legend.position = "right"
    )

ggsave(file.path(OUT, "Fig_lipid_bubble.pdf"), p3, width=11, height=9, dpi=300)
cat("  ✓ Fig_lipid_bubble.pdf\n")


# ============================================================
# Figure 4: Kennedy pathway 示意 + 扩缩标注
# ============================================================
cat("绘制Figure 4: Kennedy pathway...\n")

# 用ggplot手绘简化版Kennedy pathway
pathway_steps <- data.frame(
    step = c("Acetyl-CoA", "Malonyl-CoA", "C16:0-ACP", "C18:0-ACP", "C18:1-ACP",
             "Free FA", "G3P", "LPA", "PA", "DAG", "TAG", "Oil body"),
    x = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12),
    y = c(5, 5, 5, 5, 5, 5, 3, 3, 3, 3, 3, 3),
    stringsAsFactors = FALSE
)

enzymes <- data.frame(
    enzyme = c("ACC","KAS","KAS","SAD","Thioesterase","GPAT","LPAAT","PAP","DGAT"),
    Gene_change = c("", "", "", "-3 (SAD)", "+1", "+1 (GPAT)", "", "", "+5 (DGAT)"),
    x_start = c(1.5, 2.5, 3.5, 4.5, 5.5, 7.5, 8.5, 9.5, 10.5),
    y_pos = c(5.5, 5.5, 5.5, 5.5, 5.5, 3.5, 3.5, 3.5, 3.5),
    color = c("grey50","grey50","grey50","#1565C0","#D32F2F","#D32F2F","grey50","grey50","#D32F2F"),
    stringsAsFactors = FALSE
)

regulators <- data.frame(
    name = c("WRI1", "Oleosin", "Caleosin", "Desaturase\n(×5 expanded)", "FAD\n(×4 expanded)"),
    x = c(6, 12.5, 12.5, 4, 5),
    y = c(6.5, 3.5, 2.5, 6.5, 6.5),
    change = c("+5", "+2", "+1", "+20", "+6"),
    stringsAsFactors = FALSE
)

p4 <- ggplot() +
    # 步骤节点
    geom_point(data=pathway_steps, aes(x=x, y=y), 
               shape=21, size=8, fill="lightyellow", color="grey40") +
    geom_text(data=pathway_steps, aes(x=x, y=y, label=step),
              size=2.2, fontface="bold") +
    
    # 箭头 (FA synthesis chain)
    geom_segment(data=data.frame(x=1:5, xend=2:6), 
                 aes(x=x+0.3, xend=xend-0.3, y=5, yend=5),
                 arrow=arrow(length=unit(0.15,"cm")), color="grey50", size=0.5) +
    # 箭头 (Kennedy pathway)
    geom_segment(data=data.frame(x=7:10, xend=8:11),
                 aes(x=x+0.3, xend=xend-0.3, y=3, yend=3),
                 arrow=arrow(length=unit(0.15,"cm")), color="grey50", size=0.5) +
    # Free FA → G3P (转折)
    geom_segment(aes(x=6, xend=7, y=4.6, yend=3.4),
                 arrow=arrow(length=unit(0.15,"cm")), color="grey50", 
                 size=0.5, linetype="dashed") +
    # TAG → Oil body
    geom_segment(aes(x=11.3, xend=11.7, y=3, yend=3),
                 arrow=arrow(length=unit(0.15,"cm")), color="grey50", size=0.5) +
    
    # 酶标注
    geom_label(data=enzymes, aes(x=x_start, y=y_pos, label=enzyme),
               size=2.5, fill="white", label.size=0.3, color=enzymes$color,
               fontface="bold") +
    # 扩缩数字
    geom_text(data=enzymes %>% filter(Gene_change != ""), 
              aes(x=x_start, y=y_pos+0.4, label=Gene_change),
              size=2.5, color=enzymes$color[enzymes$Gene_change != ""], fontface="bold") +
    
    # 调控因子
    geom_label(data=regulators, aes(x=x, y=y, label=paste0(name, "\n", change)),
               size=2.5, fill="#FFF3E0", color="#E65100", label.size=0.5, fontface="bold") +
    # WRI1 → pathway 调控箭头
    geom_curve(aes(x=6, y=6.1, xend=3, yend=5.3),
               arrow=arrow(length=unit(0.15,"cm")), color="#E65100",
               curvature=-0.3, size=0.4, linetype="dashed") +
    
    # 标题和主题
    labs(title="Fatty acid biosynthesis & TAG assembly pathway",
         subtitle="Red = expanded on oil palm ancestral branch | Blue = contracted") +
    xlim(0.5, 13.5) + ylim(1.5, 7.5) +
    theme_void() +
    theme(
        plot.title = element_text(size=13, face="bold"),
        plot.subtitle = element_text(size=10, color="grey40")
    )

ggsave(file.path(OUT, "Fig_lipid_pathway.pdf"), p4, width=14, height=6, dpi=300)
cat("  ✓ Fig_lipid_pathway.pdf\n")


# ============================================================
# Figure 5: 饼图对比 - expand中油脂占比
# ============================================================
cat("绘制Figure 5: 占比对比...\n")

ratio_data <- data.frame(
    Layer = rep(c("Global (30 spp)", "Palm (12 spp)", "Oil crops (16 spp)"), each=2),
    Type = rep(c("Lipid-related", "Other"), 3),
    Expand = c(114, 770-114, 75, 478-75, 52, 300-52),
    stringsAsFactors = FALSE
)

ratio_data$Layer <- factor(ratio_data$Layer, 
    levels=c("Global (30 spp)", "Palm (12 spp)", "Oil crops (16 spp)"))

p5 <- ggplot(ratio_data, aes(x=Layer, y=Expand, fill=Type)) +
    geom_col(position="fill", width=0.6) +
    scale_fill_manual(values=c("Lipid-related"="#D32F2F", "Other"="#E0E0E0")) +
    scale_y_continuous(labels=scales::percent) +
    geom_text(data=ratio_data %>% filter(Type=="Lipid-related"),
              aes(label=sprintf("%d (%.1f%%)", Expand, Expand/(Expand + 
                  ratio_data$Expand[ratio_data$Layer==Layer & ratio_data$Type=="Other"])*100)),
              y=0.08, color="white", fontface="bold", size=3.5) +
    labs(x="", y="Proportion of expanded gene families",
         title="Lipid gene enrichment in expanded families",
         fill="") +
    theme_minimal(base_size=12) +
    theme(
        plot.title = element_text(size=13, face="bold"),
        legend.position = "bottom",
        panel.grid.major.x = element_blank()
    )

ggsave(file.path(OUT, "Fig_lipid_ratio.pdf"), p5, width=8, height=5, dpi=300)
cat("  ✓ Fig_lipid_ratio.pdf\n")


cat("\n========================================\n")
cat("完成! 输出:\n")
cat("  ", file.path(OUT, "Fig_lipid_contrast.pdf"), " - 全局收缩vs油脂扩张\n")
cat("  ", file.path(OUT, "Fig_lipid_heatmap.pdf"), "  - 核心基因热图\n")
cat("  ", file.path(OUT, "Fig_lipid_bubble.pdf"), "   - 气泡图\n")
cat("  ", file.path(OUT, "Fig_lipid_pathway.pdf"), "  - Kennedy通路\n")
cat("  ", file.path(OUT, "Fig_lipid_ratio.pdf"), "    - 占比图\n")
cat("========================================\n")
