import pandas as pd, numpy as np
W="${CLUSTER_WORK}/snp_repair_rerun/gwas"
d=pd.read_csv(W+"/shell/out/finemap_shell_new.tsv",sep="\t")
cols=["trait","n","n_snp_window","n_sv_window","snp_lead","snp_lead_P_SNPmodel","snp_lead_MAF","CS95_size_SNPmodel_f0.2","CS95_nSV_SNPmodel_f0.2","PP_snplead_SNPmodel_f0.2","top_variant_SNPmodel_f0.2","top_PP_SNPmodel_f0.2","sv_lead","sv_lead_P_SVmodel","condA_P_snp_given_sv","r2_snp_sv"]
pd.set_option("display.width",300); print(d[cols].to_string())
s=pd.read_csv(W+"/shell/out/shell_snp_lead_conditioned_on_each_SV_new.tsv",sep="\t")
n=int(open(W+"/geno/n_tests.txt").read())
for tr in ["Nut_weight_g"]:
    nw=s[s.trait==tr]; print(tr,"SVs conditioned",len(nw),"max P",nw.P_snp_lead_given_sv.max(),"all sig",bool((nw.P_snp_lead_given_sv<0.05/n).all()))
o=pd.read_csv("${CLUSTER_WORK}/enh_sv_finemap/out/finemap_all.tsv",sep="\t")
o=o[(o.trait=="Nut_weight_g")&(o.chrom=="chr01B")]
print("PUBLISHED:",o[[c for c in cols if c in o.columns]].to_string())
