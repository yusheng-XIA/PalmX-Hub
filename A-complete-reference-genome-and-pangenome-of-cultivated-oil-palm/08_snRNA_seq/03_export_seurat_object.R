# Export counts, Harmony embedding and SNN graph from the snRNA Seurat object for 04_/05_ (Python)
suppressPackageStartupMessages(library(Matrix))
args <- commandArgs(trailingOnly = TRUE)   # usage: Rscript 03_export_seurat_object.R sce.all_int.rds export_dir
out <- args[2]; dir.create(out, showWarnings=FALSE, recursive=TRUE)
obj <- readRDS(args[1])   # merged Harmony-integrated Seurat object from 02_seurat_analysis.R
meta <- slot(obj, "meta.data")
cat("graphs:", names(slot(obj, "graphs")), "\n"); cat("reductions:", names(slot(obj, "reductions")), "\n")
cmd <- slot(obj, "commands"); cat("commands:", names(cmd), "\n")
for (k in names(cmd)) { p <- slot(cmd[[k]], "params"); cat("==", k, "\n"); str(p[setdiff(names(p), c("object"))]) }
rna <- slot(obj, "assays")[["RNA"]]; layers <- slot(rna, "layers")
fn <- attributes(slot(rna, "features"))$dimnames[[1]]; cn <- attributes(slot(rna, "cells"))$dimnames[[1]]
cnt <- layers[["counts"]]
cnt2 <- sparseMatrix(i = slot(cnt, "i") + 1, p = slot(cnt, "p"), x = slot(cnt, "x"), dims = slot(cnt, "Dim"))
writeMM(cnt2, file.path(out, "counts.mtx")); writeLines(fn, file.path(out, "genes.txt")); writeLines(cn, file.path(out, "cells.txt"))
hv <- slot(rna, "meta.data"); cat("rna meta.data cols:", colnames(hv), "\n")
for (r in names(slot(obj, "reductions"))) {
  e <- slot(slot(obj, "reductions")[[r]], "cell.embeddings")
  write.table(e, file.path(out, paste0("emb_", r, ".tsv")), sep="\t", quote=FALSE, col.names=NA)
}
for (g in names(slot(obj, "graphs"))) {
  G <- slot(obj, "graphs")[[g]]
  T <- summary(sparseMatrix(i = slot(G, "i") + 1, p = slot(G, "p"), x = slot(G, "x"), dims = slot(G, "Dim")))
  cat(g, "nnz", nrow(T), "dimnames match cells:", identical(slot(G, "Dimnames")[[1]], rownames(meta)), "\n")
  write.table(T, gzfile(file.path(out, paste0("graph_", g, ".tsv.gz"))), sep="\t", quote=FALSE, row.names=FALSE)
  writeLines(slot(G, "Dimnames")[[1]], file.path(out, paste0("graph_", g, "_names.txt")))
}
write.table(data.frame(cell=rownames(meta), meta, check.names=FALSE), file.path(out, "meta.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
cat("identical cells order:", identical(cn, rownames(meta)), "\n"); cat("done\n")
