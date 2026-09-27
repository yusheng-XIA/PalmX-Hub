#!/bin/bash
OF=${ANALYSIS_DIR}/20_results/Figure2/07_new_figure/02_comparative_orthofinder/OrthoFinder_Results/Results_Feb05
ls -la $OF/Orthogroups/
head -1 $OF/Orthogroups/Orthogroups.GeneCount.tsv | tr '\t' '\n' | nl
wc -l $OF/Orthogroups/Orthogroups.GeneCount.tsv
python3 - <<'PY'
import csv,sys
OF="${ANALYSIS_DIR}/20_results/Figure2/07_new_figure/02_comparative_orthofinder/OrthoFinder_Results/Results_Feb05"
f=open(OF+"/Orthogroups/Orthogroups.tsv"); h=f.readline().rstrip("\n").split("\t")
ex={}
for line in f:
    r=line.rstrip("\n").split("\t")
    for i,c in enumerate(r[1:],1):
        if c and h[i] not in ex: ex[h[i]]=c.split(", ")[:3]
    if len(ex)==len(h)-1: break
for k,v in ex.items(): print(k,v)
PY
ls $OF/Orthogroups/ | head; ls ${ANALYSIS_DIR}/20_results/Figure2/07_new_figure/02_comparative_orthofinder/data_clean/
