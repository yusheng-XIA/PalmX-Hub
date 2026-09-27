#!/bin/bash
# ${COMPUTE_HOST}; reads user (read-only), writes work/enh_fig2stats/data
W=${CLUSTER_WORK}/enh_fig2stats; mkdir -p $W/data
OF=${ANALYSIS_DIR}/20_results/Figure2/07_new_figure/02_comparative_orthofinder/OrthoFinder_Results/Results_Feb05
cp $OF/Orthogroups/Orthogroups.GeneCount.tsv $W/data/
md5sum $OF/Orthogroups/Orthogroups.GeneCount.tsv $OF/Orthogroups/Orthogroups.tsv > $W/data/input_md5.txt
python3 - <<'PY'
import re
OF="${ANALYSIS_DIR}/20_results/Figure2/07_new_figure/02_comparative_orthofinder/OrthoFinder_Results/Results_Feb05"
W="${CLUSTER_WORK}/enh_fig2stats/data"
f=open(OF+"/Orthogroups/Orthogroups.tsv"); h=f.readline().rstrip("\n").split("\t")
sp=h[1:]
out=open(W+"/og_locus_counts.tsv","w")
out.write("Orthogroup\t"+"\t".join(sp)+"\n")
def gene(s,p):
    p=p.split("|",1)[1]
    if s=="Cocos_nucifera": return re.sub(r"\.\d+$","",p)
    return p
for line in f:
    r=line.rstrip("\n").split("\t")
    r+= [""]*(len(h)-len(r))
    cnt=[]
    for s,c in zip(sp,r[1:]):
        cnt.append(str(len({gene(s,x) for x in c.split(", ")}) if c else 0))
    out.write(r[0]+"\t"+"\t".join(cnt)+"\n")
out.close()
PY
RUN=${ANALYSIS_DIR}/22_answer_reviews/00_ms/03_V3/03_figure3/05_multiomics_integration/runs/RUN-MULTIOMICS-INTEGRATION-20260721-001/outputs
cp $RUN/stage26_full_chemical_phenotype_atlas_attempt002/molecular_phenotype_scores_by_sample.tsv $W/data/
md5sum $RUN/stage26_full_chemical_phenotype_atlas_attempt002/molecular_phenotype_scores_by_sample.tsv >> $W/data/input_md5.txt
ls -la $W/data; echo DONE
