#!/usr/bin/env Rscript
#============================================
# ggtree绘制: 分化时间树 + 所有节点扩缩标注
# 参考教程: 用treeio/tidytree导入外部数据
#============================================

if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager", repos="https://cloud.r-project.org")
for (pkg in c("ggtree", "treeio"))
    if (!requireNamespace(pkg, quietly=TRUE)) BiocManager::install(pkg)
for (pkg in c("ggplot2", "ape", "dplyr", "tidytree"))
    if (!requireNamespace(pkg, quietly=TRUE))
        install.packages(pkg, repos="https://cloud.r-project.org")

library(ggtree)
library(treeio)
library(ggplot2)
library(ape)
library(dplyr)
library(tidytree)

OUT <- "${ANALYSIS_DIR}/20_results/Figure2/07_new_figure/02_comparative_orthofinder/phylo_divtime_cafe/04_figures"

# =============================================
# 辅助函数: 从id_tree匹配节点
# =============================================

# 从id_tree.nwk提取 Name<ID> 映射
parse_id_tree <- function(id_tree_file) {
    nwk <- readLines(id_tree_file)[1]
    # tip: Name<ID>
    tips <- regmatches(nwk, gregexpr("[A-Za-z][A-Za-z_]*\\d*<\\d+>", nwk))[[1]]
    tip_df <- data.frame(
        name = gsub("<.*", "", tips),
        cafe_id = as.integer(gsub(".*<(\\d+)>", "\\1", tips)),
        stringsAsFactors = FALSE
    )
    # internal: )<ID>
    ints <- regmatches(nwk, gregexpr("\\)<\\d+>", nwk))[[1]]
    int_df <- data.frame(
        name = "",
        cafe_id = as.integer(gsub(".*<(\\d+)>", "\\1", ints)),
        stringsAsFactors = FALSE
    )
    rbind(tip_df, int_df)
}

# 通过tip descendants匹配内部节点
get_tips_under <- function(tree, node) {
    if (node <= Ntip(tree)) return(tree$tip.label[node])
    stack <- c(node)
    tips <- c()
    while(length(stack) > 0) {
        cur <- stack[length(stack)]; stack <- stack[-length(stack)]
        ch <- tree$edge[tree$edge[,1]==cur, 2]
        for (c in ch) {
            if (c <= Ntip(tree)) tips <- c(tips, tree$tip.label[c])
            else stack <- c(stack, c)
        }
    }
    sort(tips)
}

# 从CAFE5 asr.tre解析每个内部节点的tip后代
parse_asr_topology <- function(asr_file) {
    lines <- readLines(asr_file)
    tree_line <- grep("^  TREE ", lines, value=TRUE)[1]
    nwk <- sub(".*=\\s*", "", tree_line)
    nwk <- gsub("\\*", "", nwk)
    nwk <- gsub("(<\\d+>)_\\d+", "\\1", nwk)
    
    # 提取所有节点
    all_nodes <- list()
    
    # Tips: Name<ID>
    tip_matches <- gregexpr("[A-Za-z][A-Za-z_]*\\d*<\\d+>", nwk)
    for (t in regmatches(nwk, tip_matches)[[1]]) {
        name <- gsub("<.*", "", t)
        id <- as.integer(gsub(".*<(\\d+)>", "\\1", t))
        all_nodes[[as.character(id)]] <- list(name=name, is_tip=TRUE, tips=name)
    }
    
    # 用栈解析内部节点的子节点
    # 简化: 直接解析括号结构
    pos <- 1
    stack <- list()  # 每层的tip集合
    stack[[1]] <- c()
    level <- 1
    
    i <- 1
    chars <- strsplit(nwk, "")[[1]]
    n <- length(chars)
    
    while (i <= n) {
        ch <- chars[i]
        if (ch == "(") {
            level <- level + 1
            stack[[level]] <- c()
        } else if (ch == ")") {
            # 找紧跟的 <ID>
            rest <- substr(nwk, i+1, min(i+20, n))
            m <- regexpr("<(\\d+)>", rest)
            if (m > 0) {
                id_str <- regmatches(rest, m)
                id <- as.integer(gsub("[<>]", "", id_str))
                current_tips <- stack[[level]]
                all_nodes[[as.character(id)]] <- list(name="", is_tip=FALSE, tips=current_tips)
                # 传递tips到上一层
                level <- level - 1
                stack[[level]] <- c(stack[[level]], current_tips)
            } else {
                level <- level - 1
            }
        } else {
            # 检查是否是tip开头
            rest <- substr(nwk, i, min(i+100, n))
            m <- regexpr("^[A-Za-z][A-Za-z_]*\\d*<\\d+>", rest)
            if (m > 0) {
                matched <- regmatches(rest, m)
                name <- gsub("<.*", "", matched)
                stack[[level]] <- c(stack[[level]], name)
                i <- i + nchar(matched) - 1
            }
        }
        i <- i + 1
    }
    
    all_nodes
}

# 将cafe_id映射到ggtree node号
map_cafe_to_ggtree <- function(tree, asr_nodes) {
    mapping <- data.frame(cafe_id=integer(), ggtree_node=integer(), stringsAsFactors=FALSE)
    
    for (id_str in names(asr_nodes)) {
        info <- asr_nodes[[id_str]]
        cafe_id <- as.integer(id_str)
        
        if (info$is_tip) {
            gn <- which(tree$tip.label == info$name)
            if (length(gn) == 1) {
                mapping <- rbind(mapping, data.frame(cafe_id=cafe_id, ggtree_node=gn))
            }
        } else {
            # 内部节点: 通过tips匹配
            cafe_tips <- sort(info$tips)
            # 找ggtree中具有相同tip集合的节点
            for (node in (Ntip(tree)+1):(Ntip(tree)+Nnode(tree))) {
                gt_tips <- sort(get_tips_under(tree, node))
                if (identical(cafe_tips, gt_tips)) {
                    mapping <- rbind(mapping, data.frame(cafe_id=cafe_id, ggtree_node=node))
                    break
                }
            }
        }
    }
    mapping
}


# =============================================
# 绘图函数
# =============================================

plot_cafe_tree <- function(tree, node_data, cafe_mapping, title, outfile,
                           width=14, height=11, use_sig=FALSE) {
    
    # 合并数据
    plot_data <- merge(cafe_mapping, node_data, by.x="cafe_id", by.y="node_id")
    
    if (use_sig) {
        plot_data$expand <- plot_data$sig_expand
        plot_data$contract <- plot_data$sig_contract
        plot_data$net <- plot_data$sig_net
    } else {
        plot_data$expand <- plot_data$all_expand
        plot_data$contract <- plot_data$all_contract
        plot_data$net <- plot_data$all_net
    }
    
    # 分类
    oil_palms <- c("American_hap1", "Dura", "Pisifera")
    palms <- c("Calamus", "Daemonorops", "Nypa_fruticans", "Phoenix_dactylifera",
               "Areca_catechu", "Cocos_nucifera", oil_palms)
    monocots <- c(palms, "Oryza_sativa", "Musa_acuminata", "Musa_balbisiana", "Acorus_calamus")
    
    grp <- ifelse(tree$tip.label %in% oil_palms, "Oil palm",
           ifelse(tree$tip.label %in% palms, "Arecaceae",
           ifelse(tree$tip.label %in% monocots, "Monocot", "Eudicot")))
    names(grp) <- tree$tip.label
    grp_col <- c("Oil palm"="#D32F2F", "Arecaceae"="#2E7D32",
                  "Monocot"="#1565C0", "Eudicot"="#455A64")
    
    # 构建tidytree
    tree_tbl <- as_tibble(tree)
    tree_tbl$display_label <- ifelse(is.na(tree_tbl$label), "", gsub("_", " ", tree_tbl$label))
    
    # 合并扩缩数据
    tree_tbl <- tree_tbl %>%
        left_join(plot_data %>% select(ggtree_node, expand, contract, net, is_tip),
                  by=c("node"="ggtree_node"))
    
    # 标注文字
    tree_tbl$ec_label <- ifelse(!is.na(tree_tbl$expand),
                                paste0("+", tree_tbl$expand, " / -", tree_tbl$contract),
                                "")
    tree_tbl$net_color <- ifelse(is.na(tree_tbl$net), "grey50",
                          ifelse(tree_tbl$net > 0, "#C62828", "#1565C0"))
    tree_tbl$node_size <- ifelse(is.na(tree_tbl$net), 0,
                          pmin(6, pmax(1.5, log10(abs(tree_tbl$net)+1)*2)))
    
    # 转回treedata
    tree_td <- as.treedata(tree_tbl)
    
    root_age <- max(node.depth.edgelength(tree))
    x_max <- root_age * 1.45
    
    # 基础树
    p <- ggtree(tree_td, size=0.7, ladderize=TRUE) +
        theme_tree2() +
        xlab("Divergence time (Mya)")
    
    # 时间轴 (反向)
    breaks <- seq(0, floor(root_age/20)*20, 20)
    p <- p + scale_x_continuous(
        breaks=breaks,
        labels=rev(breaks),
        limits=c(0, x_max)
    )
    
    # 地质年代背景
    if (root_age > 100) {
        p <- p +
            annotate("rect", xmin=root_age-23.03, xmax=root_age, 
                     ymin=-Inf, ymax=Inf, fill="#FFF9C4", alpha=0.3) +
            annotate("rect", xmin=root_age-66, xmax=root_age-23.03, 
                     ymin=-Inf, ymax=Inf, fill="#F0F4C3", alpha=0.25) +
            annotate("rect", xmin=max(0, root_age-145.5), xmax=root_age-66, 
                     ymin=-Inf, ymax=Inf, fill="#DCEDC8", alpha=0.2)
    }
    
    # 高亮油棕clade
    op_tips <- intersect(oil_palms, tree$tip.label)
    pm_tips <- intersect(palms, tree$tip.label)
    if (length(op_tips) >= 2) {
        op_mrca <- getMRCA(tree, which(tree$tip.label %in% op_tips))
        p <- p + geom_hilight(node=op_mrca, fill="#FFCDD2", alpha=0.3, extend=x_max*0.3)
    }
    if (length(pm_tips) >= 2) {
        pm_mrca <- getMRCA(tree, which(tree$tip.label %in% pm_tips))
        p <- p + geom_hilight(node=pm_mrca, fill="#C8E6C9", alpha=0.2, extend=x_max*0.3)
    }
    
    # Tip标签 (斜体+颜色)
    tip_df <- data.frame(
        node = 1:Ntip(tree),
        label = tree$tip.label,
        group = grp[tree$tip.label],
        stringsAsFactors = FALSE
    )
    p <- p + geom_tiplab(
        data=tip_df,
        aes(label=gsub("_", " ", label), color=group),
        fontface="italic", size=2.8, offset=root_age*0.02
    ) + scale_color_manual(values=grp_col, guide="none")
    
    # === 所有节点标注扩缩 ===
    
    # 内部节点: 圆点 + 数字
    int_nodes <- plot_data %>% filter(is_tip == "False" | is_tip == FALSE)
    for (i in 1:nrow(int_nodes)) {
        r <- int_nodes[i, ]
        clr <- ifelse(r$net > 0, "#C62828", ifelse(r$net < 0, "#1565C0", "grey50"))
        sz <- min(5.5, max(1.8, log10(abs(r$net)+1)*2))
        
        p <- p +
            geom_point2(aes(subset=(node==r$ggtree_node)),
                        shape=21, size=sz, fill=clr, color="white", 
                        stroke=0.3, alpha=0.85) +
            # 扩张数 (红)
            geom_text2(aes(subset=(node==r$ggtree_node)),
                       label=paste0("+", r$expand),
                       hjust=-0.1, vjust=-0.8, size=1.8, 
                       color="#C62828", fontface="bold") +
            # 收缩数 (蓝)
            geom_text2(aes(subset=(node==r$ggtree_node)),
                       label=paste0("-", r$contract),
                       hjust=-0.1, vjust=1.5, size=1.8,
                       color="#1565C0", fontface="bold")
    }
    
    # Tip节点: 在标签后面标注
    tip_nodes <- plot_data %>% filter(is_tip == "True" | is_tip == TRUE)
    for (i in 1:nrow(tip_nodes)) {
        r <- tip_nodes[i, ]
        clr <- ifelse(r$net > 0, "#C62828", "#1565C0")
        offset_x <- root_age * 0.28  # tip标签后面
        p <- p +
            geom_text2(aes(subset=(node==r$ggtree_node)),
                       label=sprintf("+%d/-%d", r$expand, r$contract),
                       hjust=-3.5, size=1.6, color="grey30")
    }
    
    # 图例
    p <- p +
        annotate("point", x=root_age*0.95, y=3, shape=21, size=4, fill="#C62828", color="white") +
        annotate("text", x=root_age*0.97, y=3, label="Net expansion", size=2.3, hjust=0) +
        annotate("point", x=root_age*0.95, y=1.5, shape=21, size=4, fill="#1565C0", color="white") +
        annotate("text", x=root_age*0.97, y=1.5, label="Net contraction", size=2.3, hjust=0) +
        ggtitle(title) +
        theme(plot.title=element_text(size=12, face="bold"),
              axis.text.x=element_text(size=8))
    
    ggsave(outfile, p, width=width, height=height, dpi=300)
    cat("  ✓", outfile, "\n")
}


# =============================================
# 主流程
# =============================================
cat("开始绘图...\n\n")

# === Layer1: 全局30物种 ===
cat("--- Layer1 Global ---\n")
tree1 <- read.tree(file.path(OUT, "mcmctree_clean.nwk"))
tree1$edge.length <- tree1$edge.length * 100

node_data1 <- read.delim(file.path(OUT, "layer1_global_node_data.tsv"))
asr_file1 <- "${ANALYSIS_DIR}/20_results/Figure2/07_new_figure/02_comparative_orthofinder/phylo_divtime_cafe/03_cafe5/layer1_global/gamma_results/Gamma_asr.tre"
asr_nodes1 <- parse_asr_topology(asr_file1)
mapping1 <- map_cafe_to_ggtree(tree1, asr_nodes1)
cat("  映射:", nrow(mapping1), "nodes\n")

# 全部家族版本
plot_cafe_tree(tree1, node_data1, mapping1,
               "Gene family evolution across 30 species (all families)",
               file.path(OUT, "Fig2a_global_all.pdf"),
               width=16, height=12, use_sig=FALSE)

# 显著家族版本
plot_cafe_tree(tree1, node_data1, mapping1,
               "Gene family evolution across 30 species (significant families, p<0.05)",
               file.path(OUT, "Fig2a_global_sig.pdf"),
               width=16, height=12, use_sig=TRUE)


# === Layer2: 棕榈科 ===
cat("\n--- Layer2 Palm ---\n")
palm_spp <- c("Calamus", "Daemonorops", "Nypa_fruticans", "Phoenix_dactylifera",
              "Areca_catechu", "Cocos_nucifera", "American_hap1", "Dura", "Pisifera",
              "Oryza_sativa", "Musa_acuminata", "Musa_balbisiana")
tree2 <- drop.tip(tree1, tree1$tip.label[!tree1$tip.label %in% palm_spp])

node_data2 <- read.delim(file.path(OUT, "layer2_palm_node_data.tsv"))
asr_file2 <- "${ANALYSIS_DIR}/20_results/Figure2/07_new_figure/02_comparative_orthofinder/phylo_divtime_cafe/03_cafe5/layer2_palm/gamma_results/Gamma_asr.tre"
asr_nodes2 <- parse_asr_topology(asr_file2)
mapping2 <- map_cafe_to_ggtree(tree2, asr_nodes2)
cat("  映射:", nrow(mapping2), "nodes\n")

plot_cafe_tree(tree2, node_data2, mapping2,
               "Gene family evolution in Arecaceae (all families)",
               file.path(OUT, "Fig2b_palm_all.pdf"),
               width=13, height=8, use_sig=FALSE)

plot_cafe_tree(tree2, node_data2, mapping2,
               "Gene family evolution in Arecaceae (significant, p<0.05)",
               file.path(OUT, "Fig2b_palm_sig.pdf"),
               width=13, height=8, use_sig=TRUE)


# === Layer3: 产油作物 ===
cat("\n--- Layer3 OilCrop ---\n")
# Layer3有不同的树拓扑，需要读chronos树
layer3_tree_file <- "${ANALYSIS_DIR}/20_results/Figure2/07_new_figure/02_comparative_orthofinder/phylo_divtime_cafe/03_cafe5/layer3_oilcrop/cafe_tree.nwk"
tree3 <- read.tree(layer3_tree_file)

node_data3 <- read.delim(file.path(OUT, "layer3_oilcrop_node_data.tsv"))
asr_file3 <- "${ANALYSIS_DIR}/20_results/Figure2/07_new_figure/02_comparative_orthofinder/phylo_divtime_cafe/03_cafe5/layer3_oilcrop/gamma_results/Gamma_asr.tre"
asr_nodes3 <- parse_asr_topology(asr_file3)
mapping3 <- map_cafe_to_ggtree(tree3, asr_nodes3)
cat("  映射:", nrow(mapping3), "nodes\n")

plot_cafe_tree(tree3, node_data3, mapping3,
               "Gene family evolution in oil crops (all families)",
               file.path(OUT, "Fig2c_oilcrop_all.pdf"),
               width=14, height=9, use_sig=FALSE)

cat("\n========================================\n")
cat("完成! 输出目录:", OUT, "\n")
cat("  Fig2a_global_all.pdf   - 全局30物种 (全部家族)\n")
cat("  Fig2a_global_sig.pdf   - 全局30物种 (显著家族)\n") 
cat("  Fig2b_palm_all.pdf     - 棕榈科 (全部家族)\n")
cat("  Fig2b_palm_sig.pdf     - 棕榈科 (显著家族)\n")
cat("  Fig2c_oilcrop_all.pdf  - 产油作物 (全部家族)\n")
cat("========================================\n")
