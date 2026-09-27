#!/usr/bin/env python3
"""Part 2/3: DNA allele ratios at FL diagnostic sites; DNA-corrected RNA ASE; genome vs modules; key FA genes."""
import sys, numpy as np, pandas as pd
from pathlib import Path
from scipy import stats
W = Path("${CLUSTER_WORK}/enh_A"); O = W/"out"
A = Path("${ANALYSIS_DIR}/22_answer_reviews")
UNI = A/"00_ms/03_V3/03_figure3/05_multiomics_integration/runs/RUN-MULTIOMICS-ALLELE-CNS-V4-20260723-001/outputs/stage3_existing_ASE_unification_attempt001/gene_stage_ASE_unified.tsv.gz"
SL = A/"00_ms/03_V3/03_figure3/01_ASE/02_seedless/ase/config/validated_one_to_one_gene_pairs.tsv"
rng = np.random.default_rng(20260924); NPERM = 5000
L = open(O/"part2_log.txt", "w")
def P(*a):
    s = " ".join(map(str, a)); print(s); L.write(s + "\n"); L.flush()

s = pd.read_csv(W/"sites/FL_sites_liftB.tsv.gz", sep="\t")
def load(side, src):
    d = pd.read_csv(W/f"pileup/{side}.{src}.counts.tsv.gz", sep="\t")
    return d
for src, tag in [("seedless3-10", "sr"), ("hifi_remap", "hf")]:
    a = load("A", src).rename(columns={"chrom": "chrom", "pos": "pos"})
    m = s.merge(a, on=["chrom", "pos"], how="left").fillna({"A": 0, "C": 0, "G": 0, "T": 0})
    tot = m[["A", "C", "G", "T"]].sum(1)
    s[f"{tag}A_a"] = [r[x] for r, x in zip(m[["A","C","G","T"]].to_dict("records"), m.ref)]
    s[f"{tag}A_b"] = [r[x] for r, x in zip(m[["A","C","G","T"]].to_dict("records"), m.alt)]
    s[f"{tag}A_o"] = tot - s[f"{tag}A_a"] - s[f"{tag}A_b"]
    b = load("B", src).rename(columns={"chrom": "chrom_B", "pos": "pos_B"})
    mb = s[["chrom_B", "pos_B", "ref_onB", "alt_onB", "mapq_B"]].merge(b, on=["chrom_B", "pos_B"], how="left").fillna({"A": 0, "C": 0, "G": 0, "T": 0})
    rec = mb[["A","C","G","T"]].to_dict("records")
    ok = (mb.mapq_B >= 20).values
    s[f"{tag}B_a"] = np.where(ok, [r[x] if isinstance(x, str) else 0 for r, x in zip(rec, mb.ref_onB)], np.nan)
    s[f"{tag}B_b"] = np.where(ok, [r[x] if isinstance(x, str) else 0 for r, x in zip(rec, mb.alt_onB)], np.nan)
s.to_csv(O/"site_DNA_counts.tsv.gz", sep="\t", index=False, compression="gzip")

# ---- site-level QC/distribution ----
for tag, nm in [("sr", "short-read"), ("hf", "HiFi")]:
    for side in ["A", "B"]:
        a, b = s[f"{tag}{side}_a"], s[f"{tag}{side}_b"]; d = a + b
        ok = d >= 10
        fr = (a / d)[ok]
        P(f"[{nm} -> {side}-ref] sites depth>=10: {ok.sum()} / {s[f'{tag}{side}_a'].notna().sum()}; median depth {d[ok].median():.0f};"
          f" mean A-fraction {fr.mean():.4f}; median {fr.median():.4f}; A-only(<5% B) {100*(fr>0.95).mean():.2f}%; B-only(<5% A) {100*(fr<0.05).mean():.2f}%;"
          f" 0.3-0.7 {100*((fr>=0.3)&(fr<=0.7)).mean():.2f}%")

# ---- gene-level DNA ratios (sum over sites; sites with DNA depth within 5..3x median) ----
def gene_ratio(tag, side, retained_only=False):
    a, b = s[f"{tag}{side}_a"], s[f"{tag}{side}_b"]; d = a + b
    med = d[d > 0].median()
    keep = (d >= 5) & (d <= 3 * med)
    if retained_only: keep &= s.retained.astype(str).str.upper().eq("TRUE")
    g = pd.DataFrame({"gene_id": s.gene_id[keep], "a": a[keep], "b": b[keep]}).groupby("gene_id").agg(a=("a", "sum"), b=("b", "sum"), n=("a", "size"))
    g[f"l2_{tag}{side}"] = np.log2((g.a + 0.5) / (g.b + 0.5))
    return g.rename(columns={"a": f"a_{tag}{side}", "b": f"b_{tag}{side}", "n": f"nsite_{tag}{side}"})
G = None
for tag in ["hf", "sr"]:
    for side in ["A", "B"]:
        g = gene_ratio(tag, side); G = g if G is None else G.join(g, how="outer")
for tag in ["hf", "sr"]:
    G[f"l2_{tag}_sym"] = (G[f"l2_{tag}A"] + G[f"l2_{tag}B"]) / 2      # reciprocal-reference average cancels reference bias
    G[f"refbias_{tag}"] = (G[f"l2_{tag}A"] - G[f"l2_{tag}B"]) / 2
G.to_csv(O/"gene_DNA_ratios.tsv.gz", sep="\t", compression="gzip")
for c in ["l2_hfA", "l2_hfB", "l2_hf_sym", "refbias_hf", "l2_srA", "l2_srB", "l2_sr_sym", "refbias_sr"]:
    v = G[c].dropna(); P(f"gene DNA {c}: n={len(v)} median={v.median():.3f} mean={v.mean():.3f} IQR=[{v.quantile(.25):.3f},{v.quantile(.75):.3f}] |x|>0.5: {100*(v.abs()>0.5).mean():.1f}%  <-0.5: {100*(v< -0.5).mean():.1f}%  >0.5: {100*(v>0.5).mean():.1f}%")

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
            null = np.array([np.median(vals[rng.choice(len(vals), k, replace=False)]) for _ in range(NPERM)])
            rows.append(dict(version=label, window=w, set=m, n_eligible=int(x.gene_id.isin(gs).sum()), n_robust=int(k),
                             pct_B=100*(vals[inm] < 0).mean(), median_log2AB=obs,
                             p_perm_lower=(1 + (null <= obs).sum())/(NPERM + 1),
                             p_wilcox=stats.mannwhitneyu(vals[inm], vals[~inm]).pvalue))
    return rows

rows = []
# uncorrected on the subset with DNA (same genes as corrected) to isolate the effect of correction
for dcol, lab in [("l2_hf_sym", "HiFi reciprocal"), ("l2_hfA", "HiFi A-ref only"), ("l2_sr_sym", "Short-read reciprocal"), ("l2_srA", "Short-read A-ref only")]:
    x = recall(r, dcol)
    rows += compare(gene_table(x, "log2_allele_ratio", "robust_ase"), f"RNA uncorrected (genes with {lab} DNA)")
    rows += compare(gene_table(x, "l2c", "robust_c"), f"RNA - DNA ({lab})")
    if dcol == "l2_hf_sym":
        x.to_csv(O/"FL_gene_stage_DNAcorrected_hifi_sym.tsv.gz", sep="\t", index=False, compression="gzip")
        P("HiFi-sym corrected: robust gene-stage rows", int(x.robust_c.sum()), "vs original", int(x.robust_ase.sum()),
          "; B-biased share original", round(100*(x[x.robust_ase].log2_allele_ratio < 0).mean(), 2), "corrected", round(100*(x[x.robust_c].l2c < 0).mean(), 2))
# DNA-balanced genes only (both HiFi A-ref and B-ref |log2|<0.25, >=2 sites)
Gb = G[(G.l2_hfA.abs() < 0.25) & (G.l2_hfB.abs() < 0.25) & (G.nsite_hfA >= 2)].index
xb = r[r.gene_id.isin(Gb)]
rows += compare(gene_table(xb, "log2_allele_ratio", "robust_ase"), "RNA uncorrected, DNA-balanced genes only")
res = pd.DataFrame(rows); res.to_csv(O/"part2_corrected_genome_vs_module.tsv", sep="\t", index=False, float_format="%.5g")
pd.set_option("display.width", 250); pd.set_option("display.max_rows", 500)
P(res.to_string())
