#!/usr/bin/env python3
"""DNA-corrected FL ASE (Extended Data Fig. 5).

1. Allele counts at each diagnostic site are read from the FL-Hap2 (A) and FL-Hap1 (B) short-read pileups
   (05b_fl_dna_allele_counts.sh). Sites with DNA depth 5..3x the median depth are summed per gene and a gene
   DNA log2(A/B) is computed on each reference; the two values are averaged (reciprocal-reference average).
2. For every eligible FL gene-stage, the RNA log2(A/B) is corrected by subtracting the gene's DNA log2(A/B) and
   re-tested with the score test of 04_ase_test.py using the DNA-derived allelic proportion as the null
   (clipped to 0.05-0.95); BH within stage; robust = FDR < 0.05 and |corrected log2(A/B)| >= 0.5.
3. Genome-wide share of B-biased genes and trait-module medians are compared before and after correction
   (optional --modules: columns trait_module, gene_id), with 10,000 random gene sets of equal size.
"""
import argparse

import numpy as np
import pandas as pd
from scipy import stats

rng = np.random.default_rng(20260924)
NPERM = 10000
BASES = ["A", "C", "G", "T"]
WINDOWS = {**{i: "Days 0-65" for i in range(1, 6)}, **{i: "Days 80-140" for i in range(6, 11)},
           **{i: "Days 155-185" for i in range(11, 14)}, **{i: "Hours 12-72" for i in range(14, 20)}}


def bh(p):
    p = np.asarray(p, float); n = len(p); o = np.argsort(p)
    q = np.minimum.accumulate((p[o] * n / np.arange(1, n + 1))[::-1])[::-1]
    out = np.empty(n); out[o] = np.minimum(q, 1); return out


def site_counts(sites, pileup_dir):
    a = pd.read_csv(f"{pileup_dir}/A.short_read.counts.tsv.gz", sep="\t")
    m = sites.merge(a, on=["chrom", "pos"], how="left").fillna({b: 0 for b in BASES})
    rec = m[BASES].to_dict("records")
    sites["A_a"] = [r[x] for r, x in zip(rec, m.ref)]
    sites["A_b"] = [r[x] for r, x in zip(rec, m.alt)]
    b = pd.read_csv(f"{pileup_dir}/B.short_read.counts.tsv.gz", sep="\t").rename(columns={"chrom": "chrom_B", "pos": "pos_B"})
    mb = sites[["chrom_B", "pos_B", "ref_onB", "alt_onB", "mapq_B"]].merge(b, on=["chrom_B", "pos_B"], how="left").fillna({x: 0 for x in BASES})
    rec = mb[BASES].to_dict("records")
    ok = (mb.mapq_B >= 20).values
    sites["B_a"] = np.where(ok, [r[x] if isinstance(x, str) else 0 for r, x in zip(rec, mb.ref_onB)], np.nan)
    sites["B_b"] = np.where(ok, [r[x] if isinstance(x, str) else 0 for r, x in zip(rec, mb.alt_onB)], np.nan)
    return sites


def gene_ratio(s, side):
    a, b = s[f"{side}_a"], s[f"{side}_b"]; d = a + b
    keep = (d >= 5) & (d <= 3 * d[d > 0].median())
    g = pd.DataFrame({"gene_id": s.gene_id[keep], "a": a[keep], "b": b[keep]}).groupby("gene_id").agg(
        a=("a", "sum"), b=("b", "sum"), n=("a", "size"))
    g[f"l2_{side}"] = np.log2((g.a + 0.5) / (g.b + 0.5))
    return g.rename(columns={"a": f"a_{side}", "b": f"b_{side}", "n": f"nsite_{side}"})


def recall(df, dcol):
    x = df.dropna(subset=[dcol]).copy()
    p0 = np.clip(2 ** x[dcol] / (1 + 2 ** x[dcol]), 0.05, 0.95)
    z = (x.ref_fragments - p0 * x.informative_fragments) / np.sqrt(p0 * (1 - p0) * x.variance_weight)
    x["p_c"] = 2 * stats.norm.sf(np.abs(z)); x["padj_c"] = np.nan
    for _, idx in x.groupby("stage").groups.items():
        x.loc[idx, "padj_c"] = bh(x.loc[idx, "p_c"].values)
    x["l2c"] = x.log2_allele_ratio - x[dcol]
    x["robust_c"] = (x.padj_c < 0.05) & (x.l2c.abs() >= 0.5)
    return x


def compare(x, lcol, rcol, modules, label):
    """Per developmental window (and all stages): share of B-biased genes and trait-module medians."""
    rows = []
    x = pd.concat([x.assign(window=x.stage_index.map(WINDOWS)), x.assign(window="All stages")])
    for w, xw in x.groupby("window"):
        g = xw.groupby("gene_id").agg(robust_any=(rcol, "any")).reset_index()
        med = xw[xw[rcol]].groupby("gene_id")[lcol].median().rename("med").reset_index()
        rr = g.merge(med)[lambda d: d.robust_any].reset_index(drop=True)
        vals = rr.med.values
        if len(vals) == 0:
            continue
        n_b = int((vals < 0).sum())
        rows.append(dict(version=label, window=w, set="Genome-wide", n_robust=len(rr), pct_B=100 * n_b / len(vals),
                         p_binom=stats.binomtest(n_b, len(vals), 0.5).pvalue, median_log2AB=np.median(vals)))
        for m, gs in (modules or {}).items():
            inm = rr.gene_id.isin(gs).values; k = int(inm.sum())
            if k < 3:
                continue
            obs = np.median(vals[inm])
            null = np.median(vals[rng.integers(0, len(vals), (NPERM, k))], axis=1)
            p_lo = (1 + (null <= obs).sum()) / (NPERM + 1); p_hi = (1 + (null >= obs).sum()) / (NPERM + 1)
            rows.append(dict(version=label, window=w, set=m, n_robust=k, pct_B=100 * (vals[inm] < 0).mean(),
                             median_log2AB=obs, p_perm_two_sided=min(1.0, 2 * min(p_lo, p_hi))))
    out = pd.DataFrame(rows)
    mod = out.set.ne("Genome-wide")
    if mod.any():   # BH across module-window combinations
        out.loc[mod, "padj_perm"] = bh(out.loc[mod, "p_perm_two_sided"].values)
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sites", required=True, help="FL_sites_liftB.tsv.gz from 05a_lift_sites_to_hap1.py")
    ap.add_argument("--pileup-dir", required=True)
    ap.add_argument("--ase", required=True, help="gene_stage_ASE.tsv.gz from 04_ase_test.py")
    ap.add_argument("--modules", help="optional trait-module membership (trait_module, gene_id)")
    ap.add_argument("--out-prefix", required=True)
    a = ap.parse_args()

    s = site_counts(pd.read_csv(a.sites, sep="\t"), a.pileup_dir)
    G = gene_ratio(s, "A").join(gene_ratio(s, "B"), how="outer")
    G["l2_sym"] = (G.l2_A + G.l2_B) / 2          # reciprocal-reference average
    G["refbias"] = (G.l2_A - G.l2_B) / 2
    G.to_csv(f"{a.out_prefix}.gene_DNA_ratios.tsv.gz", sep="\t")

    r = pd.read_csv(a.ase, sep="\t", low_memory=False)
    r = r[(r.analysis == "FL") & r.eligible.astype(str).str.lower().eq("true")]
    r["robust_ase"] = r.robust_ase.astype(str).str.lower().eq("true")
    r = r.merge(G[["l2_sym"]], left_on="gene_id", right_index=True, how="left")
    x = recall(r, "l2_sym")
    x.to_csv(f"{a.out_prefix}.gene_stage_DNA_corrected.tsv.gz", sep="\t", index=False)

    modules = None
    if a.modules:
        mem = pd.read_csv(a.modules, sep="\t")
        modules = mem.groupby("trait_module").gene_id.apply(set).to_dict()
    res = pd.concat([compare(x, "log2_allele_ratio", "robust_ase", modules, "RNA uncorrected"),
                     compare(x, "l2c", "robust_c", modules, "RNA - DNA")])
    res.to_csv(f"{a.out_prefix}.genome_vs_module.tsv", sep="\t", index=False, float_format="%.5g")


if __name__ == "__main__":
    main()
