"""Protein-level deployment of LD-coat families and seed-storage proteins (directLFQ matrix, 114 samples)."""
import os, numpy as np, pandas as pd
B="${ANALYSIS_DIR}"
P=B+"/22_answer_reviews/00_ms/03_V3/03_figure3/07_proteomics_multiomics_20260724"
R=P+"/01_reference_database"; O="out/prot/"; os.makedirs(O,exist_ok=True)
mat=pd.read_csv(P+"/02_current_Astral_114/runs/RUN-PROT-ASTRAL-DIRECTLFQ-20260725-001/outputs/current114_directlfq_protein_abundance_integration_keys.tsv",sep="\t",low_memory=False)
sc=[c for c in mat.columns if "|" in c]
mp=pd.read_csv(R+"/FL_TN_unified_exact_sequence_nr_map.tsv",sep="\t",dtype=str).set_index("unified_protein_id")
ann=pd.read_csv(B+"/20_results/Figure2/07_new_figure/05_omic/1.final_counts/GO_annotation/Africa_hap2/Africa_hap2.emapper.annotations",sep="\t",comment="#",header=None,dtype=str,low_memory=False).set_index(0)
desc=(ann[7].fillna("")+" | "+ann[8].fillna("")+" | "+ann[20].fillna(""))
LD=pd.read_csv("out/dom/ld_candidates_domains.tsv",sep="\t",dtype=str); fam=dict(zip(LD.unified_id,LD.family))
def ahap2(u):
    gs=[g.split(":")[1] for g in str(mp.loc[u,"source_gene_ids"]).split(";") if g.startswith("Africa_hap2:")] if u in mp.index else []
    return gs
recs=[]
for r in mat.itertuples(index=False):
    grp=str(r[0]); us=grp.split(";")
    fams={fam[u] for u in us if u in fam}
    gs=sorted({g for u in us for g in ahap2(u)})
    d=" || ".join(desc.get(g,"")[:60] for g in gs[:2])
    st=any(pd.Series([desc.get(g,"") for g in gs]).str.contains(r"vicilin|legumin|globulin|glutelin|2S albumin|seed storage|cupincin",case=False,regex=True)) if gs else False
    if fams or st:
        recs.append(dict(group=grp,family=";".join(sorted(fams)) if fams else "seed_storage",genes=";".join(gs),desc=d))
G=pd.DataFrame(recs); print(G.family.value_counts())
M=mat.set_index(mat.columns[0]).loc[G.group,sc].apply(pd.to_numeric,errors="coerce")
M=M.where(M>0)
out=G.set_index("group").join(M); out.to_csv(O+"ld_storage_groups_by_sample.tsv",sep="\t")
