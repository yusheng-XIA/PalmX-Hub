G=${ANALYSIS_DIR}/05_GWAS/00_analysis/02_genomeDB/05_snp
awk "{n[\$1]++; if(\$4>m[\$1])m[\$1]=\$4} END{for(c in n) print c, n[c], m[c]}" $G/LD_pruned_bed.bim | sort -V
cat $G/run_evolution_pipeline.sh
