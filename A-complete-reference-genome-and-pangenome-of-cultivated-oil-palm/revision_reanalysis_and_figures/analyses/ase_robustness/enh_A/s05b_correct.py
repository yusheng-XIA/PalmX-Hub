#!/usr/bin/env python3
"""Part 2b: DNA-corrected RNA ASE, genome vs modules (reads gene_DNA_ratios.tsv.gz from s05_dna.py)."""
import sys, resource, numpy as np, pandas as pd
from pathlib import Path
from scipy import stats
W = Path("${CLUSTER_WORK}/enh_A"); O = W/"out"
A = Path("${ANALYSIS_DIR}/22_answer_reviews")
UNI = A/"00_ms/03_V3/03_figure3/05_multiomics_integration/runs/RUN-MULTIOMICS-ALLELE-CNS-V4-20260723-001/outputs/stage3_existing_ASE_unification_attempt001/gene_stage_ASE_unified.tsv.gz"
rng = np.random.default_rng(20260924); NPERM = 10000
L = open(O/"part2b_log.txt", "w")
def P(*a):
    s = " ".join(map(str, a)); print(s, flush=True); L.write(s + "\n"); L.flush()
G = pd.read_csv(O/"gene_DNA_ratios.tsv.gz", sep="\t", index_col=0)
# ---- RNA ----
use = ["analysis","gene_id","stage","stage_group","eligible","robust_ase","log2_allele_ratio","ref_fragments","alt_fragments","informative_fragments","variance_weight","padj"]
r = pd.read_csv(UNI, sep="\t", usecols=use, low_memory=False)
r = r[(r.analysis == "FL")]
for c in ["eligible","robust_ase"]: r[c] = r[c].astype(str).str.lower().eq("true")
r = r[r.eligible].merge(G[["l2_hf_sym","l2_hfA","l2_sr_sym","l2_srA","nsite_hfA","nsite_hfB"]], left_on="gene_id", right_index=True, how="left")
P("FL eligible gene-stage rows", len(r), "genes", r.gene_id.nunique(), "with HiFi-sym DNA ratio", r.dropna(subset=["l2_hf_sym"]).gene_id.nunique())

def bh(p):
    p = np.asarray(p, float); n = len(p); o = np.argsort(p); q = p[o] * n / np.arange(1, n + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]; out = np.empty(n); out[o] = np.minimum(q, 1); return out

def recall(df, dcol):
    """shift log2 by DNA log2 and re-test with a gene-specific null p0 = DNA A-fraction."""
    x = df.dropna(subset=[dcol]).copy()
    p0 = np.clip(2**x[dcol] / (1 + 2**x[dcol]), 0.05, 0.95)
    z = (x.ref_fragments - p0 * x.informative_fragments) / np.sqrt(p0 * (1 - p0) * x.variance_weight)
    x["p_c"] = 2 * stats.norm.sf(np.abs(z)); x["padj_c"] = np.nan
    for st, idx in x.groupby("stage").groups.items(): x.loc[idx, "padj_c"] = bh(x.loc[idx, "p_c"].values)
    x["l2c"] = x.log2_allele_ratio - x[dcol]
    x["robust_c"] = (x.padj_c < 0.05) & (x.l2c.abs() >= 0.5)
    return x

mem = pd.read_csv(O/"module_membership.tsv", sep="\t"); mem = mem[mem.analysis == "FL"]
mods = mem.groupby("trait_module").gene_id.apply(set)
WIN = ["Days 0–65", "Days 80–140", "Days 155–185", "Hours 12–72", "All stages"]

def gene_table(x, lcol, rcol):
    x2 = x.copy(); x3 = x.copy(); x3["stage_group"] = "All stages"; x2 = pd.concat([x2, x3])
    g = x2.groupby(["stage_group", "gene_id"]).agg(robust_any=(rcol, "any"), mean_depth=("informative_fragments", "mean")).reset_index()
    rb = x2[x2[rcol]].groupby(["stage_group", "gene_id"])[lcol].median().rename("med").reset_index()
    return g.merge(rb, how="left")

def compare(gt, label):
    rows = []
    for w in WIN:
        x = gt[gt.stage_group == w]; rr = x[x.robust_any].reset_index(drop=True); vals = rr.med.values
        rows.append(dict(version=label, window=w, set="Genome-wide", n_eligible=len(x), n_robust=len(rr),
                         pct_B=100*(vals < 0).mean(), median_log2AB=np.median(vals)))
        for m, gs in mods.items():
            inm = rr.gene_id.isin(gs).values; k = inm.sum()
            if k < 3: continue
            obs = np.median(vals[inm])
            null = np.median(vals[rng.integers(0, len(vals), (NPERM, k))], axis=1)
            rows.append(dict(version=label, window=w, set=m, n_eligible=int(x.gene_id.isin(gs).sum()), n_robust=int(k),
                             pct_B=100*(vals[inm] < 0).mean(), median_log2AB=obs,
                             p_perm_lower=(1 + (null <= obs).sum())/(NPERM + 1),
                             p_wilcox=stats.mannwhitneyu(vals[inm], vals[~inm], method="asymptotic").pvalue))
    return rows

rows = []
# uncorrected on the subset with DNA (same genes as corrected) to isolate the effect of correction
for dcol, lab in [("l2_sr_sym", "Short-read reciprocal"), ("l2_srA", "Short-read A-ref only"), ("l2_hf_sym", "HiFi reciprocal")]:
    P("version", lab, "maxrss GB", resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1e6)
    x = recall(r, dcol); P(" recall done", resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1e6, x.padj_c.isna().sum())
    gt = gene_table(x, "log2_allele_ratio", "robust_ase"); P(" gt done", len(gt), resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1e6)
    rows += compare(gt, f"RNA uncorrected (genes with {lab} DNA)"); P(" cmp done", resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1e6)
    rows += compare(gene_table(x, "l2c", "robust_c"), f"RNA - DNA ({lab})")
    if dcol == "l2_sr_sym":
        x.to_csv(O/"FL_gene_stage_DNAcorrected_sr_sym.tsv.gz", sep="\t", index=False, compression="gzip")
        P("Short-read reciprocal corrected: robust gene-stage rows", int(x.robust_c.sum()), "vs original", int(x.robust_ase.sum()),
          "; B-biased share original", round(100*(x[x.robust_ase].log2_allele_ratio < 0).mean(), 2), "corrected", round(100*(x[x.robust_c].l2c < 0).mean(), 2))
# DNA-balanced genes only (both HiFi A-ref and B-ref |log2|<0.25, >=2 sites)
Gb = G[(G.l2_srA.abs() < 0.2) & (G.l2_srB.abs() < 0.2) & (G.nsite_srA >= 2)].index
P("DNA-balanced genes (short-read |log2|<0.2 on both references):", len(Gb))
xb = r[r.gene_id.isin(Gb)]
rows += compare(gene_table(xb, "log2_allele_ratio", "robust_ase"), "RNA uncorrected, DNA-balanced genes only")
res = pd.DataFrame(rows); res.to_csv(O/"part2_corrected_genome_vs_module.tsv", sep="\t", index=False, float_format="%.5g")
pd.set_option("display.width", 250); pd.set_option("display.max_rows", 500)
P(res.to_string())
P("maxrss GB", resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1e6)
