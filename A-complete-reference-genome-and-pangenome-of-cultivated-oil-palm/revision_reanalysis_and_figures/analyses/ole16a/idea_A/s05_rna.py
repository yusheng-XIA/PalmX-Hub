"""Gate 3a: bulk RNA-seq of oleosin / caleosin / LDAP(REF) genes and seed-tissue marker programmes (114 libraries)."""
import os, re, numpy as np, pandas as pd
B="${ANALYSIS_DIR}"
RUN=B+"/22_answer_reviews/00_ms/03_V3/03_figure3/05_multiomics_integration/runs/RUN-MULTIOMICS-INTEGRATION-20260721-001/outputs"
O="out/rna/"; os.makedirs(O,exist_ok=True)
STAGES=["0d","15d","35d","50d","65d","80d","95d","110d","125d","140d","155d","170d","185d","12h","24h","36h","48h","60h","72h"]
cw=pd.read_csv(RUN+"/stage1_identity/sample_crosswalk_114.tsv",sep="\t",dtype=str)
counts=pd.read_csv(B+"/22_answer_reviews/00_ms/03_V3/03_figure3/00_minipan/03_rnaseq_mapping/count_matrices/pangraphrna_hisat2_graph.gene_counts.tsv",sep="\t",index_col=0)
S=cw.rna_sample.tolist(); C=counts[S].astype(float)
sfu=pd.read_csv(B+"/22_answer_reviews/00_ms/05_MS/MS_revision_3/Final_20260823_VectorRevision/8月29日/10_Table1_ST8_AI_Cleanup_20260902/ST9_latest_reanalysis_audit/RNA_114_joint_size_factors.tsv",sep="\t").set_index("sample").DESeq2_size_factor
N=C.div(sfu.loc[S],axis=1)
# gene length for TPM-like within-sample rank: use counts per million instead (rank within library)
cpm=C.div(C.sum(0),axis=1)*1e6
ann=pd.read_csv(B+"/20_results/Figure2/07_new_figure/05_omic/1.final_counts/GO_annotation/Africa_hap2/Africa_hap2.emapper.annotations",sep="\t",comment="#",header=None,dtype=str,low_memory=False)
ann=ann.set_index(0); desc=(ann[7].fillna("")+" | "+ann[8].fillna("")+" | "+ann[20].fillna(""))
LD=pd.read_csv("out/dom/ld_candidates_domains.tsv",sep="\t",dtype=str)
g2f={}
for r in LD.itertuples():
    for g in r.genes.split(";"):
        if g.startswith("Africa_hap2:"): g2f[g.split(":")[1]]=r.family
# seed / embryo programme markers from annotation (defined before looking at data)
pat={"seed_storage_globulin":r"vicilin|legumin|11S globulin|7S globulin|Cupin_1|Cupin_2|glutelin|2S albumin|seed storage",
     "LEA_dehydrin":r"late embryogenesis|\bLEA\b|LEA_|Dehydrin|dehydrin",
     "seed_master_TF":r"\bABI3\b|\bLEC1\b|\bLEC2\b|\bFUS3\b|\bVP1\b|B3 domain.*ABI3"}
mk={}
for k,p in pat.items():
    ids=desc.index[desc.str.contains(p,regex=True,case=False)]
    ids=[i for i in ids if i in C.index and i not in g2f]; mk[k]=ids
    print(k,len(ids))
rows=[]
genes=[g for g in g2f if g in C.index]
print("LD genes in counts",len(genes),"of",len(g2f))
cwi=cw.set_index("rna_sample")
long=[]
for g in genes:
    for s in S:
        long.append(dict(gene=g,family=g2f[g],name=desc.get(g,"").split(" | ")[1] if g in desc.index else "",sample=s,genotype=cwi.loc[s,"genotype"],stage=cwi.loc[s,"stage"],norm=N.loc[g,s],cpm=cpm.loc[g,s],raw=C.loc[g,s]))
L=pd.DataFrame(long); L.to_csv(O+"ld_gene_sample_long.tsv.gz",sep="\t",index=False)
M=L.groupby(["family","gene","name","genotype","stage"]).norm.mean().unstack("stage")[STAGES]
M.to_csv(O+"ld_gene_stage_mean_norm.tsv",sep="\t")
# within-library rank of each OLE16 gene (percentile among all genes by CPM)
pct=cpm.rank(pct=True)
rk=[]
for g in ["evm.TU.chr11B.1497","evm.TU.chr04B.697","evm.TU.chr03B.2346"]:
    for s in S:
        rk.append(dict(gene=g,sample=s,genotype=cwi.loc[s,"genotype"],stage=cwi.loc[s,"stage"],cpm=cpm.loc[g,s],pctile=pct.loc[g,s]))
pd.DataFrame(rk).groupby(["gene","genotype","stage"])[["cpm","pctile"]].mean().unstack("genotype").to_csv(O+"ole16_cpm_percentile.tsv",sep="\t")
# marker programme sums
mrows=[]
for k,ids in mk.items():
    tot=N.loc[ids].sum(0)
    for s in S: mrows.append(dict(programme=k,n_genes=len(ids),sample=s,genotype=cwi.loc[s,"genotype"],stage=cwi.loc[s,"stage"],norm_sum=tot[s]))
MR=pd.DataFrame(mrows); MR.groupby(["programme","genotype","stage"]).norm_sum.mean().unstack("stage")[STAGES].to_csv(O+"seed_programme_stage_means.tsv",sep="\t")
# top individual seed-storage genes at 170-185 d FL
for k,ids in mk.items():
    sub=N.loc[ids]; fl=[s for s in S if cwi.loc[s,"genotype"]=="FL" and cwi.loc[s,"stage"] in ("170d","185d")]
    tn=[s for s in S if cwi.loc[s,"genotype"]=="TN" and cwi.loc[s,"stage"] in ("170d","185d")]
    t=pd.DataFrame({"FL_170_185":sub[fl].mean(1),"TN_170_185":sub[tn].mean(1),"max_any":sub.max(1),"desc":[desc.get(i,"")[:80] for i in ids]}).sort_values("max_any",ascending=False).head(15)
    t.to_csv(O+f"top_{k}.tsv",sep="\t")
# Pearson/Spearman of OLE16 with other LD genes across FL samples
print(M.round(1).to_string())
