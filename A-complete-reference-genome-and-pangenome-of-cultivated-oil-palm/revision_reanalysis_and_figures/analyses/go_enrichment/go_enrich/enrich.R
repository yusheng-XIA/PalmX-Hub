suppressPackageStartupMessages({library(clusterProfiler); library(GO.db)})
args <- commandArgs(TRUE); mode <- args[1]
t2g <- read.delim("term2gene.tsv")
tt <- AnnotationDbi::select(GO.db, keys=unique(t2g$term), columns=c("TERM","ONTOLOGY"), keytype="GOID")
t2g <- t2g[t2g$term %in% tt$GOID[!is.na(tt$TERM)], ]
t2n <- tt[, c("GOID","TERM")]
run <- function(sets, universe, tag) {
  universe <- intersect(universe, t2g$gene)
  out <- list()
  for (s in names(sets)) {
    g <- intersect(sets[[s]], universe)
    e <- enricher(g, universe=universe, TERM2GENE=t2g, TERM2NAME=t2n, pAdjustMethod="BH",
                  pvalueCutoff=0.05, qvalueCutoff=1, minGSSize=10, maxGSSize=500)
    r <- if (is.null(e)) NULL else as.data.frame(e)
    cat(tag, s, "genes", length(sets[[s]]), "with GO", length(g), "sig", if (is.null(r)) 0 else sum(r$p.adjust < 0.05), "\n")
    if (!is.null(r) && nrow(r)) { r <- r[r$p.adjust < 0.05, ]; if (nrow(r)) { r$set <- s; out[[s]] <- r } }
  }
  res <- do.call(rbind, out)
  if (!is.null(res)) { res$ontology <- tt$ONTOLOGY[match(res$ID, tt$GOID)]; cat(tag, "universe", length(universe), "\n") }
  res
}
if (mode == "wgcna") {
  m <- read.delim("wgcna_modules.tsv")
  sets <- split(m$gene, m$module); sets <- sets[names(sets) != "grey"]
  res <- run(sets, m$gene, "wgcna")
  write.table(res, "go_wgcna_modules.tsv", sep="\t", quote=FALSE, row.names=FALSE)
} else {
  mk <- read.delim("sn_markers_all.tsv"); det <- read.delim("sn_gene_detection.tsv")
  sig <- mk[mk$p_val_adj < 0.05 & mk$avg_log2FC >= 1 & mk$pct.1 >= 0.1, ]
  sets <- split(sig$gene, paste0("C", sig$cluster))
  universe <- unique(mk$gene)  # genes tested in the marker analysis
  res <- run(sets, universe, "sn")
  write.table(res, "go_sn_clusters.tsv", sep="\t", quote=FALSE, row.names=FALSE)
}
cat("clusterProfiler", as.character(packageVersion("clusterProfiler")), "\n")
