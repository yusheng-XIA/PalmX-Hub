"""Gate 2: peptide-level requantification of lipid-droplet coat proteins from DIA-NN report.parquet."""
import re, numpy as np, pandas as pd, pyarrow.parquet as pq
from collections import defaultdict
P="${ANALYSIS_DIR}/22_answer_reviews/00_ms/03_V3/03_figure3/07_proteomics_multiomics_20260724"
F=P+"/02_current_Astral_114/runs/RUN-PROT-ASTRAL-FULL-20260724-006/final_pass/"
R=P+"/01_reference_database"
O="out/pep/"; import os; os.makedirs(O,exist_ok=True)
def readfa(p):
    d={};k=None
    for l in open(p):
        l=l.strip()
        if l.startswith('>'): k=l[1:].split()[0]; d[k]=[]
        elif k: d[k].append(l)
    return {k:''.join(v).rstrip('*') for k,v in d.items()}
seq=readfa(R+"/FL_TN_unified_exact_sequence_nr.fasta")
mp=pd.read_csv(R+"/FL_TN_unified_exact_sequence_nr_map.tsv",sep="\t",dtype=str).set_index("unified_protein_id")
C=pd.read_csv("out/dom/ld_candidates_domains.tsv",sep="\t",dtype=str)
cand=set(C.unified_id)
fam=dict(zip(C.unified_id,C.family))
mani=pd.read_csv(P+"/00_contract_manifest/current_sample_manifest_114.tsv",sep="\t",dtype=str)
mani["run"]=mani.raw_basename.str.replace(".raw","",regex=False)
r2k=dict(zip(mani.run,mani.integration_key))
avail=pq.ParquetFile(F+"report.parquet").schema.names
print(avail)
want=["Run","Protein.Group","Protein.Ids","Protein.Names","Genes","Stripped.Sequence","Modified.Sequence","Precursor.Id","Proteotypic",
      "Precursor.Quantity","Precursor.Normalised","Q.Value","Global.Q.Value","PG.Q.Value","Lib.Q.Value","Lib.PG.Q.Value","Global.PG.Q.Value","PG.MaxLFQ","RT","Precursor.Charge"]
cols=[c for c in want if c in avail]
t=pq.read_table(F+"report.parquet",columns=cols).to_pandas()
print("rows",len(t))
t["run"]=t.Run.astype(str).str.replace(".raw","",regex=False).map(lambda x: os.path.basename(x))
t["key"]=t.run.map(r2k); print("unmapped runs",t.key.isna().sum())
# acceptance identical to directLFQ input: all q <= 0.01, quantity > 0
m=np.ones(len(t),bool)
for q in ["Q.Value","Global.Q.Value","PG.Q.Value","Lib.PG.Q.Value","Global.PG.Q.Value"]:
    if q in t: m&=t[q].values<=0.01
m&=t["Precursor.Quantity"].values>0
t["accepted"]=m
# global run-level intensity distribution for detection-limit context
acc=t[m]
q=acc.groupby("key")["Precursor.Quantity"].quantile([0.01,0.05,0.25,0.5]).unstack(); q.columns=["p01","p05","p25","p50"]
q.to_csv(O+"run_precursor_quantiles.tsv",sep="\t")
# LD rows (any row whose Protein.Ids touch a candidate)
ids=t["Protein.Ids"].astype(str)
hit=ids.map(lambda s: any(x in cand for x in s.split(";")))
L=t[hit].copy(); print("LD rows",len(L),"accepted",L.accepted.sum())
# in-silico peptide mapping (I/L-equivalent) to all unified proteins -> loci & materials
allp=list(seq.items())
def norm(s): return s.replace("I","L")
nseq={k:norm(v) for k,v in allp}
peps=L["Stripped.Sequence"].unique()
pm={}
for p in peps:
    np_=norm(p); hits=[k for k,v in nseq.items() if np_ in v]
    genes=set(); mats=set()
    for h in hits:
        for g in str(mp.loc[h,"source_gene_ids"]).split(";"): genes.add(g)
        for x in str(mp.loc[h,"materials"]).split(";"): mats.add(x)
    pm[p]=dict(n_prot=len(hits),prots=";".join(sorted(hits)),genes=";".join(sorted(genes)),materials=";".join(sorted(mats)),
               fams=";".join(sorted({fam.get(h,"other") for h in hits})))
PM=pd.DataFrame(pm).T; PM.index.name="peptide"; PM.to_csv(O+"ld_peptide_mapping.tsv",sep="\t")
L=L.join(PM,on="Stripped.Sequence")
L.to_csv(O+"ld_precursor_rows.tsv.gz",sep="\t",index=False)
print(PM.to_string())
