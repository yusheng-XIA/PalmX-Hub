B=${ANALYSIS_DIR}/05_GWAS/00_analysis/06_GWAS/shared_genotype/GP_maf_allChr.bim
O=${CLUSTER_WORK}/snp_repair_r308/out
awk "{n[\$1]++; if(\$4>m[\$1])m[\$1]=\$4; b=int(\$4/5e6); h[\$1\" \"b]++} END{for(c in n) print c, n[c], m[c] > \"$O/gwas_bim_chrom.tsv\"; for(k in h) print k, h[k] > \"$O/gwas_bim_5Mb.tsv\"}" $B
