import numpy as np, pandas as pd, scipy.io
X="${CLUSTER_WORK}/enh_sn_pop/sn/export/"
B="${ANALYSIS_DIR}"
meta=pd.read_csv(X+"meta.tsv",sep="\t",index_col=0); genes=np.array(open(X+"genes.txt").read().split())
C=scipy.io.mmread(X+"counts.mtx").tocsc()
LAB={"E_1_95_2":"FL_95","E_1_125_4":"FL_125","E_1_185_9":"FL_185","E_4_95_5":"TN_95","E_4_125_4":"TN_125","E_4_185_5":"TN_185"}
meta["lib"]=meta["orig.ident"].map(LAB); cl=meta.seurat_clusters.astype(int).values
print(pd.crosstab(meta.seurat_clusters,meta.lib).to_string())
ann=pd.read_csv(B+"/20_results/Figure2/07_new_figure/05_omic/1.final_counts/GO_annotation/Africa_hap2/Africa_hap2.emapper.annotations",sep="\t",comment="#",header=None,dtype=str,low_memory=False).set_index(0)
desc=(ann[7].fillna("")+" | "+ann[8].fillna(""))
tot=np.asarray(C.sum(0)).ravel()
Cn=C.multiply(1e4/tot).tocsr()
lg=Cn.copy(); lg.data=np.log1p(lg.data)
det=(C>0).tocsr()
out=[]
for k in [18,12,7,9,1]:
    a=cl==k; b=~a
    ma=np.asarray(lg[:,a].mean(1)).ravel(); mb=np.asarray(lg[:,b].mean(1)).ravel()
    pa=np.asarray(det[:,a].mean(1)).ravel(); pb=np.asarray(det[:,b].mean(1)).ravel()
    fc=(ma-mb)/np.log(2); sel=np.where((pa>0.25))[0]; o=sel[np.argsort(-fc[sel])][:20]
    for i in o: out.append(dict(cl=k,gene=genes[i],lfc=round(fc[i],2),pct_in=round(pa[i],3),pct_out=round(pb[i],3),desc=desc.get(genes[i],"")[:90]))
D=pd.DataFrame(out); D.to_csv("out/sn/cluster_markers_C18_C12_C7_C9_C1.tsv",sep="\t",index=False)
pd.set_option("display.width",250); pd.set_option("display.max_colwidth",90); print(D.to_string())
