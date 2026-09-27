# GO enrichment used throughout the study (clusterProfiler v4.10.1; v4.20.0 for the WGCNA-module and
# single-nucleus-cluster enrichments). GO terms from eggNOG-mapper; the background is restricted to genes tested
# in the corresponding analysis that carry GO annotations; BH-adjusted P < 0.05.
library(clusterProfiler)
args <- commandArgs(trailingOnly = TRUE)   # gene_list.txt background.txt gene2go.tsv go_terms.tsv out.tsv
genes    <- readLines(args[1])
universe <- readLines(args[2])
gene2go  <- read.delim(args[3], header = FALSE, col.names = c("gene", "go"))   # one gene-GO pair per line
term2name <- read.delim(args[4], header = FALSE, col.names = c("go", "name"))
term2gene <- gene2go[gene2go$gene %in% universe, c("go", "gene")]

ego <- enricher(gene = intersect(genes, term2gene$gene), universe = unique(term2gene$gene),
                TERM2GENE = term2gene, TERM2NAME = term2name,
                pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 1)
write.table(as.data.frame(ego), args[5], sep = "\t", quote = FALSE, row.names = FALSE)
