#!/usr/bin/env Rscript
#=====================================================
# 修正版: 棕榈科 Layer2 进化树
# 修复 x 轴标签偏移问题
#=====================================================

library(ggtree)
library(ggplot2)
library(ape)
library(dplyr)

BASE <- "${ANALYSIS_DIR}/20_results/Figure2/07_new_figure/02_comparative_orthofinder/phylo_divtime_cafe"
CAFE <- file.path(BASE, "03_cafe5")
OUT  <- file.path(BASE, "04_figures")

# ============================================================
# 1. 读树 (MCMCTree × 100)
# ============================================================
cat("1. 读树...\n")
tree <- read.tree(file.path(BASE, "02_mcmctree/mcmctree_clean.nwk"))
tree$edge.length <- tree$edge.length * 100
root_age <- max(node.depth.edgelength(tree))
cat("  全局 Root:", round(root_age, 2), "Mya\n")

# ============================================================
# 2. 读 Layer2 CAFE5 数据
# ============================================================
cat("2. 读 CAFE5 Layer2...\n")
read_clade <- function(path) {
    df <- read.delim(path, header=TRUE, comment.char="#",
                     col.names=c("Taxon_ID", "Increase", "Decrease"))
    df$cafe_id <- as.integer(gsub(".*<(\\d+)>", "\\1", df$Taxon_ID))
    df$name <- trimws(gsub("<.*", "", df$Taxon_ID))
    df$net <- df$Increase - df$Decrease
    df
}
cafe_l2 <- read_clade(file.path(CAFE, "layer2_palm/gamma_results/Gamma_clade_results.txt"))

# ============================================================
# 3. 节点映射函数
# ============================================================
build_map <- function(tr, cafe_df, anc_defs) {
    rows <- list()
    for (i in 1:nrow(cafe_df)) {
        nm <- cafe_df$name[i]
        if (nm != "" && nm %in% tr$tip.label) {
            gn <- which(tr$tip.label == nm)
            rows[[length(rows)+1]] <- data.frame(
                gn=gn, cid=cafe_df$cafe_id[i],
                inc=cafe_df$Increase[i], dec=cafe_df$Decrease[i],
                net=cafe_df$net[i], tip=TRUE, stringsAsFactors=FALSE)
        }
    }
    for (ad in anc_defs) {
        cid <- ad[[1]]; t1 <- ad[[2]]; t2 <- ad[[3]]
        if (!(t1 %in% tr$tip.label) || !(t2 %in% tr$tip.label)) next
        gn <- getMRCA(tr, c(which(tr$tip.label==t1), which(tr$tip.label==t2)))
        row <- cafe_df[cafe_df$cafe_id == cid, ]
        if (nrow(row) == 1) {
            rows[[length(rows)+1]] <- data.frame(
                gn=gn, cid=cid, inc=row$Increase, dec=row$Decrease,
                net=row$net, tip=FALSE, stringsAsFactors=FALSE)
        }
    }
    do.call(rbind, rows)
}

# ============================================================
# 4. 构建棕榈科子树
# ============================================================
cat("4. 构建棕榈科子树...\n")

oil_palms <- c("American_hap1", "Dura", "Pisifera")
palms <- c("Calamus", "Daemonorops", "Nypa_fruticans", "Phoenix_dactylifera",
           "Areca_catechu", "Cocos_nucifera", oil_palms)
palm_spp <- c(palms, "Oryza_sativa", "Musa_acuminata", "Musa_balbisiana")

palm_tree <- drop.tip(tree, tree$tip.label[!tree$tip.label %in% palm_spp])
palm_root <- max(node.depth.edgelength(palm_tree))
cat("  棕榈科 Root:", round(palm_root, 2), "Mya\n")

# Layer2 节点定义
L2_anc <- list(
    list(2,"Dura","Pisifera"), list(4,"American_hap1","Dura"),
    list(6,"Cocos_nucifera","Dura"), list(8,"Areca_catechu","Cocos_nucifera"),
    list(10,"Phoenix_dactylifera","Cocos_nucifera"),
    list(14,"Nypa_fruticans","Cocos_nucifera"),
    list(15,"Calamus","Cocos_nucifera"), list(16,"Calamus","Daemonorops"),
    list(18,"Musa_acuminata","Calamus"),
    list(19,"Musa_acuminata","Musa_balbisiana"),
    list(20,"Oryza_sativa","Musa_acuminata")
)

map2 <- build_map(palm_tree, cafe_l2, L2_anc)
cat("  Layer2 映射:", nrow(map2), "nodes\n")

# ============================================================
# 5. 绘图 (修正 x 轴)
# ============================================================
cat("5. 绘图...\n")

# ★ 关键修正: x 轴标签基于实际 palm_root
# ggtree: x=0 是根, x=palm_root 是 tips (0 Mya)
# 标签公式: label = palm_root - x
# 为了让标签为整十数, 用 round(palm_root, -1) = 110 作为最大标签
max_lbl <- round(palm_root, -1)  # 110
x_max <- palm_root + 40          # tip标签后面留空间
brk_at <- seq(0, max_lbl, 10)    # x 位置: 0, 10, 20, ..., 110
lbl <- max_lbl - brk_at           # 标签: 110, 100, 90, ..., 10, 0

cat("  palm_root =", round(palm_root, 2), "\n")
cat("  breaks =", paste(brk_at, collapse=", "), "\n")
cat("  labels =", paste(lbl, collapse=", "), "\n")

p2 <- ggtree(palm_tree, size=1.2, ladderize=TRUE, layout="roundrect") +
    theme_tree2() +
    scale_x_continuous(
        breaks = brk_at,
        labels = lbl,            # ★ 修正: 基于 palm_root 的正确标签
        limits = c(0, x_max)
    ) +
    xlab("Divergence time (Mya)") +
    ggtitle("Gene family evolution in Arecaceae (Layer 2)")

# 高亮油棕
op2 <- getMRCA(palm_tree, which(palm_tree$tip.label %in% oil_palms))
p2 <- p2 + geom_hilight(node=op2, fill="#FFCDD2", alpha=0.3, extend=90)

# Tip 标签
p2 <- p2 + geom_tiplab(fontface="italic", size=5, offset=3)

# 获取 xy 坐标
td2 <- p2$data

# 内部节点标注
int2 <- map2[!map2$tip, ]
for (i in 1:nrow(int2)) {
    r <- int2[i, ]
    nd <- td2[td2$node == r$gn, ]
    if (nrow(nd) == 0) next
    xp <- nd$x[1]; yp <- nd$y[1]
    bg <- ifelse(r$net > 0, "#FFCDD2", ifelse(r$net < 0, "#BBDEFB", "#E0E0E0"))
    bd <- ifelse(r$net > 0, "#C62828", ifelse(r$net < 0, "#1565C0", "grey50"))
    sz <- min(6, max(2.5, log10(abs(r$net)+1)*2.2))

    p2 <- p2 +
        annotate("point", x=xp, y=yp, shape=21, size=sz, fill=bg, color=bd, stroke=0.5) +
        annotate("text", x=xp+1, y=yp+0.25, label=sprintf("+%d", r$inc),
                 size=3.5, color="#C62828", fontface="bold", hjust=0) +
        annotate("text", x=xp+1, y=yp-0.25, label=sprintf("-%d", r$dec),
                 size=3.5, color="#1565C0", fontface="bold", hjust=0)
}

# Tip 节点标注
tip2 <- map2[map2$tip, ]
for (i in 1:nrow(tip2)) {
    r <- tip2[i, ]
    nd <- td2[td2$node == r$gn, ]
    if (nrow(nd) == 0) next
    xp <- palm_root + 20
    yp <- nd$y[1]
    p2 <- p2 +
        annotate("text", x=xp, y=yp,
                 label=sprintf("+%d/-%d", r$inc, r$dec),
                 size=3.2, color="grey35", hjust=0)
}

# 图例
p2 <- p2 +
    annotate("point", x=x_max-25, y=2.5, shape=21, size=5, fill="#FFCDD2", color="#C62828") +
    annotate("text", x=x_max-22, y=2.5, label="Expansion", size=3.5, hjust=0) +
    annotate("point", x=x_max-25, y=1.5, shape=21, size=5, fill="#BBDEFB", color="#1565C0") +
    annotate("text", x=x_max-22, y=1.5, label="Contraction", size=3.5, hjust=0) +
    theme(plot.title=element_text(size=16, face="bold"),
          axis.text.x=element_text(size=10))

# ============================================================
# 6. 保存
# ============================================================
outfile <- file.path(OUT, "Fig2b_palm_v6_fixed.pdf")
ggsave(outfile, p2, width=14, height=8, dpi=300)
cat("\n  ✓ 已保存:", outfile, "\n")

# 验证: 打印关键节点的实际分化时间
cat("\n===== 关键节点分化时间验证 =====\n")
key_nodes <- list(
    c("Cocos_nucifera", "Dura"),
    c("Areca_catechu", "Cocos_nucifera"),
    c("Calamus", "Daemonorops"),
    c("Dura", "Pisifera")
)
for (kn in key_nodes) {
    gn <- getMRCA(palm_tree, c(which(palm_tree$tip.label==kn[1]),
                                which(palm_tree$tip.label==kn[2])))
    nd <- td2[td2$node == gn, ]
    if (nrow(nd) > 0) {
        actual_time <- palm_root - nd$x[1]
        cat(sprintf("  %s / %s → x=%.1f, 真实时间=%.2f Mya\n",
                    kn[1], kn[2], nd$x[1], actual_time))
    }
}

cat("\n完成!\n")
