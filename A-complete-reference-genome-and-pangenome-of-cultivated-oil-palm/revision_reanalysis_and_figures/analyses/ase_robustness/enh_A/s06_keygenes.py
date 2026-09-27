#!/usr/bin/env python3
"""Part 3: ASE direction/magnitude (DNA-corrected) of oleic / unsaturated-FA key genes in FL."""
import numpy as np, pandas as pd
from pathlib import Path
W = Path("${CLUSTER_WORK}/enh_A"); O = W/"out"
A = Path("${ANALYSIS_DIR}/22_answer_reviews")
FAT = A/"10_fatty_acid_pathway_reconstruction/tables/fa_pathway_gene_catalog_integrated.tsv"
SL = A/"00_ms/03_V3/03_figure3/01_ASE/02_seedless/ase/config/validated_one_to_one_gene_pairs.tsv"
KEY = ["SAD","FAD2","FAD3","FAD6","FAD7/FAD8","FATA/B","KAS I/II","KASIII","DGAT","PDAT","LPAT"]
cat = pd.read_csv(FAT, sep="\t")[["Enzyme","GeneID","Preferred_name","RNA_FL_mean_TPM"]].drop_duplicates("GeneID")
cat = cat[cat.Enzyme.isin(KEY)]
pairs = pd.read_csv(SL, sep="\t").rename(columns={"gene_a":"GeneID","gene_b":"gene_B_FLHap1"})[["GeneID","gene_B_FLHap1"]]
G = pd.read_csv(O/"gene_DNA_ratios.tsv.gz", sep="\t", index_col=0)
from scipy import stats
UNI = A/"00_ms/03_V3/03_figure3/05_multiomics_integration/runs/RUN-MULTIOMICS-ALLELE-CNS-V4-20260723-001/outputs/stage3_existing_ASE_unification_attempt001/gene_stage_ASE_unified.tsv.gz"
use = ["analysis","gene_id","stage","eligible","robust_ase","log2_allele_ratio","ref_fragments","informative_fragments","variance_weight"]
u = pd.read_csv(UNI, sep="\t", usecols=use, low_memory=False); u = u[u.analysis=="FL"]
for c in ["eligible","robust_ase"]: u[c] = u[c].astype(str).str.lower().eq("true")
u = u[u.eligible].merge(G[["l2_sr_sym"]], left_on="gene_id", right_index=True, how="inner").dropna(subset=["l2_sr_sym"]).reset_index(drop=True)
def bh(p):
    p=np.asarray(p,float); n=len(p); o=np.argsort(p); q=p[o]*n/np.arange(1,n+1); q=np.minimum.accumulate(q[::-1])[::-1]; out=np.empty(n); out[o]=np.minimum(q,1); return out
p0 = np.clip(2**u.l2_sr_sym/(1+2**u.l2_sr_sym), .05, .95)
z = (u.ref_fragments - p0*u.informative_fragments)/np.sqrt(p0*(1-p0)*u.variance_weight)
u["p_c"] = 2*stats.norm.sf(np.abs(z)); u["padj_c"] = np.nan
for st, idx in u.groupby("stage").groups.items(): u.loc[idx,"padj_c"] = bh(u.loc[idx,"p_c"].values)
u["l2c"] = u.log2_allele_ratio - u.l2_sr_sym; u["robust_c"] = (u.padj_c<0.05) & (u.l2c.abs()>=0.5)
u.to_csv(O/"FL_gene_stage_DNAcorrected_shortread_sym.tsv.gz", sep="\t", index=False, compression="gzip")
x = u
r = pd.read_csv(O/"gene_window_ASE_unified.tsv.gz", sep="\t"); r = r[(r.analysis=="FL") & (r.stage_group=="All stages")]
rows = []
for c in cat.itertuples():
    g = c.GeneID; xs = x[x.gene_id == g]; rr = r[r.gene_id == g]
    d = G.loc[g] if g in G.index else None
    row = dict(enzyme=c.Enzyme, gene_A_FLHap2=g, name=c.Preferred_name, FL_TPM=round(c.RNA_FL_mean_TPM, 1),
               n_eligible_stages=int(rr.n_elig.iloc[0]) if len(rr) else 0)
    if d is not None:
        row.update(DNA_sr_Aref=d.l2_srA, DNA_sr_Bref=d.l2_srB, DNA_sr_recip=d.l2_sr_sym, DNA_hifi_recip=d.l2_hf_sym, DNA_sites=d.nsite_srA)
    if len(rr):
        row.update(RNA_median_log2AB=rr.med_all.iloc[0], RNA_Aref_sum=int(rr.refsum.iloc[0]), RNA_Bref_sum=int(rr.altsum.iloc[0]))
    if len(xs):
        row.update(corr_median_log2AB=xs.l2c.median(), stages_robustA_orig=int((xs.robust_ase & (xs.log2_allele_ratio>0)).sum()),
                   stages_robustB_orig=int((xs.robust_ase & (xs.log2_allele_ratio<0)).sum()),
                   stages_robustA_corr=int((xs.robust_c & (xs.l2c>0)).sum()), stages_robustB_corr=int((xs.robust_c & (xs.l2c<0)).sum()))
    rows.append(row)
t = pd.DataFrame(rows).merge(pairs.rename(columns={"GeneID":"gene_A_FLHap2"}), how="left")
def verdict(r):
    if pd.isna(r.get("corr_median_log2AB")): return "not testable (no DNA-corrected ASE)"
    a, b = r.stages_robustA_corr, r.stages_robustB_corr
    if b >= 2 and a == 0: return "E. oleifera (B) allele higher"
    if a >= 2 and b == 0: return "African (A) allele higher"
    if a + b == 0: return "balanced"
    return "stage-dependent / mixed"
t["verdict_corrected"] = t.apply(verdict, axis=1)
t = t.sort_values(["enzyme","FL_TPM"], ascending=[True, False])
t.to_csv(O/"part3_key_FA_genes.tsv", sep="\t", index=False, float_format="%.3f")
pd.set_option("display.width", 300); pd.set_option("display.max_columns", 30)
print(t.to_string())
print(t.groupby("enzyme").verdict_corrected.value_counts().to_string())
