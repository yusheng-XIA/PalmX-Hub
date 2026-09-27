"""Gate 3b: OLE16 / LD-coat genes in snRNA nuclei; ambient and seed-contamination checks."""
import os, numpy as np, pandas as pd, scipy.io, scipy.sparse as sp
X="${CLUSTER_WORK}/enh_sn_pop/sn/export/"; SO="${CLUSTER_WORK}/enh_sn_pop/sn/out/"
O="out/sn/"; os.makedirs(O,exist_ok=True)
B="${ANALYSIS_DIR}"
meta=pd.read_csv(X+"meta.tsv",sep="\t",index_col=0)
genes=open(X+"genes.txt").read().split(); gi={g:i for i,g in enumerate(genes)}
C=scipy.io.mmread(X+"counts.mtx").tocsr()   # genes x cells
LAB={"E_1_95_2":"FL_95","E_1_125_4":"FL_125","E_1_185_9":"FL_185","E_4_95_5":"TN_95","E_4_125_4":"TN_125","E_4_185_5":"TN_185"}
meta["lib"]=meta["orig.ident"].map(LAB); meta["cl"]=meta.seurat_clusters.astype(int)
dbl=pd.read_csv(SO+"scrublet_per_nucleus.tsv.gz",sep="\t",index_col=0).reindex(meta.index)
meta["dbl_score"]=dbl["dbl_score"].values
ann=pd.read_csv(B+"/20_results/Figure2/07_new_figure/05_omic/1.final_counts/GO_annotation/Africa_hap2/Africa_hap2.emapper.annotations",sep="\t",comment="#",header=None,dtype=str,low_memory=False).set_index(0)
desc=(ann[7].fillna("")+" | "+ann[8].fillna("")+" | "+ann[20].fillna(""))
LD=pd.read_csv("out/dom/ld_candidates_domains.tsv",sep="\t",dtype=str)
g2f={}
for r in LD.itertuples():
    for g in r.genes.split(";"):
        if g.startswith("Africa_hap2:"): g2f[g.split(":")[1]]=r.family
storage=[g for g in desc.index[desc.str.contains(r"vicilin|legumin|globulin|glutelin|2S albumin|seed storage|cupincin",case=False,regex=True)] if g in gi]
lea=[g for g in desc.index[desc.str.contains(r"late embryogenesis abundant|\bLEA\b|Dehydrin",case=False,regex=True)] if g in gi]
print("storage genes in sn",len(storage),"LEA",len(lea))
tot=np.asarray(C.sum(0)).ravel()
def row(g): return np.asarray(C[gi[g]].todense()).ravel()
OLEa="evm.TU.chr11B.1497"; OLEb="evm.TU.chr04B.697"
print("OLE16a in sn genes:",OLEa in gi,"OLE16b:",OLEb in gi)
rows=[]
for g,f in sorted(g2f.items()):
    if g not in gi: rows.append(dict(gene=g,family=f,in_sn=False)); continue
    x=row(g)
    for L,idx in meta.groupby("lib").indices.items():
        rows.append(dict(gene=g,family=f,in_sn=True,lib=L,nuclei=len(idx),pct_pos=100*(x[idx]>0).mean(),umi_per_1e4=1e4*x[idx].sum()/tot[idx].sum(),max_umi=int(x[idx].max())))
pd.DataFrame(rows).to_csv(O+"ld_genes_by_library.tsv",sep="\t",index=False)
# per-cluster in each library for OLE16a and main LDAP
res=[]
for g in [OLEa,OLEb,"evm.TU.chr10B.1088","evm.TU.chr05B.1024","evm.TU.chr14B.647"]:
    if g not in gi: continue
    x=row(g)
    for (L,c),idx in meta.groupby(["lib","cl"]).indices.items():
        if len(idx)<30: continue
        lam=x[meta.index.get_indexer(meta.index[meta.lib==L])].sum()/tot[meta.lib==L].sum()*tot[idx]  # ambient-like expectation if uniform
        res.append(dict(gene=g,lib=L,cl=c,n=len(idx),pct_pos=100*(x[idx]>0).mean(),exp_pct_uniform=100*np.mean(1-np.exp(-lam)),umi_per_1e4=1e4*x[idx].sum()/tot[idx].sum()))
R=pd.DataFrame(res); R.to_csv(O+"ld_genes_by_cluster.tsv",sep="\t",index=False)
# OLE16a-positive nuclei characterization
x=row(OLEa); pos=x>0
meta["OLEa"]=x; meta["storage_umi"]=np.asarray(C[[gi[g] for g in storage]].sum(0)).ravel() if storage else 0
meta["lea_umi"]=np.asarray(C[[gi[g] for g in lea]].sum(0)).ravel() if lea else 0
meta["ldap_umi"]=row("evm.TU.chr10B.1088")+row("evm.TU.chr05B.1024")+row("evm.TU.chr14B.647")
summ=[]
for L,d in meta.groupby("lib"):
    p=d.OLEa>0
    summ.append(dict(lib=L,n=len(d),OLEa_pos=int(p.sum()),pct=100*p.mean(),OLEa_umi_total=int(d.OLEa.sum()),frac_lib_umi=d.OLEa.sum()/d.nCount_RNA.sum(),
        median_nUMI_pos=d.nCount_RNA[p].median() if p.any() else np.nan,median_nUMI_neg=d.nCount_RNA[~p].median(),
        dbl_pos=d.dbl_score[p].median() if p.any() else np.nan,dbl_neg=d.dbl_score[~p].median(),
        storage_per1e4_pos=1e4*d.storage_umi[p].sum()/d.nCount_RNA[p].sum() if p.any() else np.nan,storage_per1e4_neg=1e4*d.storage_umi[~p].sum()/d.nCount_RNA[~p].sum(),
        lea_per1e4_pos=1e4*d.lea_umi[p].sum()/d.nCount_RNA[p].sum() if p.any() else np.nan,lea_per1e4_neg=1e4*d.lea_umi[~p].sum()/d.nCount_RNA[~p].sum(),
        ldap_pct_pos=100*(d.ldap_umi[p]>0).mean() if p.any() else np.nan,ldap_pct_neg=100*(d.ldap_umi[~p]>0).mean(),
        top_clusters_pos=";".join(f"C{k}:{v}" for k,v in d.cl[p].value_counts().head(6).items())))
S=pd.DataFrame(summ); S.to_csv(O+"ole16a_positive_nuclei_summary.tsv",sep="\t",index=False); print(S.to_string())
# cluster composition of OLE16a+ vs all nuclei, FL_185
d=meta[meta.lib=="FL_185"]; comp=pd.DataFrame({"all":d.cl.value_counts(normalize=True),"OLEa_pos":d.cl[d.OLEa>0].value_counts(normalize=True),"n_all":d.cl.value_counts(),"n_pos":d.cl[d.OLEa>0].value_counts()}).fillna(0)
comp["pct_pos_in_cluster"]=100*comp.n_pos/comp.n_all; comp["enrich"]=comp.OLEa_pos/comp["all"]
comp.sort_values("n_all",ascending=False).to_csv(O+"FL185_cluster_OLEa.tsv",sep="\t"); print(comp.sort_values("n_all",ascending=False).round(3).to_string())
# which storage genes are expressed at all in sn (top)
st=pd.DataFrame({"gene":storage,"desc":[desc[g][:70] for g in storage],"umi_total":[int(row(g).sum()) for g in storage]}).sort_values("umi_total",ascending=False)
st.head(20).to_csv(O+"storage_genes_sn_top.tsv",sep="\t",index=False); print(st.head(10).to_string())
meta[["lib","cl","nCount_RNA","OLEa","storage_umi","lea_umi","ldap_umi","dbl_score"]].to_csv(O+"per_nucleus_ld.tsv.gz",sep="\t")
