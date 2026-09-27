#!/usr/bin/env python3
"""Replicate-aware ASE test and temporal ASE classes for the 19-stage FL and TN series.

Input
  --counts   per-library gene allele counts from 03_count_fragments.py (concatenated)
             columns: analysis, sample, gene_id, ref_fragments, alt_fragments, informative_fragments
  --samples  sample sheet: sample, analysis (FL/TN), stage, stage_index (1-19), replicate
Allele A = backbone (REF) allele: FL-Hap2 for FL, dura/TK-like for TN; allele B = alternative path.

Per replicate, at least 10 informative fragments; a gene-stage is tested when >= 2 replicates qualify
and the pooled depth is >= 30 fragments. Over-dispersion (rho) is profiled once per material and stage
by beta-binomial likelihood; the score test uses variance p0(1-p0) * sum_r n_r[1 + (n_r - 1) rho].
BH correction within each material and stage. Robust ASE: |log2((A+0.5)/(B+0.5))| >= 0.5 and FDR < 0.05;
for TN the call must also hold under p0 = 0.505 (maximum residual reference bias).

Temporal classes across stages (per gene, counting robust A- and B-biased stages):
  NoASE  no robust ASE at any qualifying stage
  HapDom >= 2 robust stages, all towards the same allele
  Sub    >= 2 A-biased and >= 2 B-biased stages
  NoDiff all other genes
"""
import argparse
import math

import numpy as np
import pandas as pd
from scipy.optimize import minimize_scalar
from scipy.special import betaln, gammaln
from scipy.stats import norm

MIN_REP_FRAGMENTS, MIN_REPLICATES, MIN_POOLED = 10, 2, 30
LOG2_MIN, FDR = 0.5, 0.05
P0 = {"FL": 0.5, "TN": 0.5}
P0_BIAS = {"FL": 0.5, "TN": 0.505}


def bh_adjust(values: pd.Series) -> pd.Series:
    arr = values.to_numpy(dtype=float)
    out = np.full(len(arr), np.nan)
    ok = np.isfinite(arr)
    p = arr[ok]
    if len(p):
        order = np.argsort(p)
        adj = np.minimum.accumulate((p[order] * len(p) / np.arange(1, len(p) + 1))[::-1])[::-1]
        restored = np.empty_like(adj)
        restored[order] = np.minimum(adj, 1.0)
        out[np.where(ok)[0]] = restored
    return pd.Series(out, index=values.index)


def beta_binomial_logpmf(k, n, p, rho):
    rho = float(np.clip(rho, 1e-7, 0.5))
    p = np.clip(p, 1e-5, 1 - 1e-5)
    conc = 1.0 / rho - 1.0
    a, b = p * conc, (1 - p) * conc
    return (gammaln(n + 1) - gammaln(k + 1) - gammaln(n - k + 1)
            + betaln(k + a, n - k + b) - betaln(a, b))


def estimate_rho(part: pd.DataFrame) -> float:
    """Profile a common stage-level rho using gene-specific pooled allele fractions."""
    p_gene = part.groupby("gene_id").apply(
        lambda x: x["ref_fragments"].sum() / x["informative_fragments"].sum(), include_groups=False)
    p = part["gene_id"].map(p_gene).to_numpy(float)
    k = part["ref_fragments"].to_numpy(float)
    n = part["informative_fragments"].to_numpy(float)
    use = (n >= MIN_REP_FRAGMENTS) & (p > 0) & (p < 1)
    if use.sum() < 100:
        return 0.01
    k, n, p = k[use], n[use], p[use]

    def objective(logit_rho):
        rho = 0.30 / (1 + math.exp(-logit_rho))
        return float(-np.sum(beta_binomial_logpmf(k, n, p, rho)))

    fit = minimize_scalar(objective, bounds=(-10, 5), method="bounded", options={"xatol": 1e-4})
    rho = 0.30 / (1 + math.exp(-float(fit.x))) if fit.success else 0.01
    return float(np.clip(rho, 1e-6, 0.30))


def compute_ase(counts: pd.DataFrame, meta: pd.DataFrame) -> pd.DataFrame:
    counts = counts.merge(meta, on=["sample", "analysis"], how="left", validate="many_to_one")
    if counts["stage"].isna().any():
        raise SystemExit("count sample absent from the sample sheet")
    keys = ["analysis", "gene_id", "stage", "stage_index"]
    all_groups = counts.groupby(keys).agg(observed_replicates=("sample", "nunique")).reset_index()

    qual = counts[counts["informative_fragments"] >= MIN_REP_FRAGMENTS].copy()
    rho = {k: estimate_rho(part) for k, part in qual.groupby(["analysis", "stage_index"])}
    qual["rho"] = [rho[k] for k in zip(qual["analysis"], qual["stage_index"])]
    qual["variance_weight"] = qual["informative_fragments"] * (1 + (qual["informative_fragments"] - 1) * qual["rho"])
    qsum = qual.groupby(keys).agg(
        qualifying_replicates=("sample", "nunique"),
        ref_fragments=("ref_fragments", "sum"),
        alt_fragments=("alt_fragments", "sum"),
        informative_fragments=("informative_fragments", "sum"),
        variance_weight=("variance_weight", "sum"),
        rho=("rho", "first"),
    ).reset_index()
    ase = all_groups.merge(qsum, on=keys, how="left")
    for c in ("qualifying_replicates", "ref_fragments", "alt_fragments", "informative_fragments"):
        ase[c] = ase[c].fillna(0).astype(np.int64)
    ase["eligible"] = (ase["qualifying_replicates"] >= MIN_REPLICATES) & (ase["informative_fragments"] >= MIN_POOLED)
    ase["log2_allele_ratio"] = np.log2((ase["ref_fragments"] + 0.5) / (ase["alt_fragments"] + 0.5))

    for suffix, p0_map in (("", P0), ("_bias_sensitivity", P0_BIAS)):
        p0 = ase["analysis"].map(p0_map).astype(float)
        z = (ase["ref_fragments"] - p0 * ase["informative_fragments"]) / np.sqrt(p0 * (1 - p0) * ase["variance_weight"])
        ase[f"z{suffix}"] = z
        ase[f"pvalue{suffix}"] = np.where(ase["eligible"], 2 * norm.sf(np.abs(z)), np.nan)
        ase[f"padj{suffix}"] = ase.groupby(["analysis", "stage_index"])[f"pvalue{suffix}"].transform(bh_adjust)

    ase["robust_ase"] = (ase["eligible"] & (ase["log2_allele_ratio"].abs() >= LOG2_MIN)
                         & (ase["padj"] < FDR) & (ase["padj_bias_sensitivity"] < FDR))
    ase["ase_call"] = np.where(~ase["eligible"], "Insufficient", "Balanced")
    ase.loc[ase["robust_ase"] & (ase["log2_allele_ratio"] > 0), "ase_call"] = "Allele_A_biased"
    ase.loc[ase["robust_ase"] & (ase["log2_allele_ratio"] < 0), "ase_call"] = "Allele_B_biased"
    return ase.sort_values(["analysis", "stage_index", "gene_id"])


def classify(n_a: int, n_b: int) -> str:
    if n_a + n_b == 0:
        return "NoASE"
    if n_a + n_b >= 2 and (n_a == 0 or n_b == 0):
        return "HapDom"
    if n_a >= 2 and n_b >= 2:
        return "Sub"
    return "NoDiff"


def temporal_classes(ase: pd.DataFrame) -> pd.DataFrame:
    tested = ase[ase["eligible"]]
    n_a = tested[tested["ase_call"] == "Allele_A_biased"].groupby(["analysis", "gene_id"]).size()
    n_b = tested[tested["ase_call"] == "Allele_B_biased"].groupby(["analysis", "gene_id"]).size()
    g = tested[["analysis", "gene_id"]].drop_duplicates().set_index(["analysis", "gene_id"])
    g["A_biased_stages"] = n_a.reindex(g.index).fillna(0).astype(int)
    g["B_biased_stages"] = n_b.reindex(g.index).fillna(0).astype(int)
    g["temporal_class"] = [classify(a, b) for a, b in zip(g["A_biased_stages"], g["B_biased_stages"])]
    return g.reset_index()


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--counts", required=True)
    ap.add_argument("--samples", required=True)
    ap.add_argument("--out-prefix", required=True)
    a = ap.parse_args()
    counts = pd.read_csv(a.counts, sep="\t", dtype={"sample": str, "gene_id": str})
    meta = pd.read_csv(a.samples, sep="\t", dtype={"sample": str})[["sample", "analysis", "stage", "stage_index", "replicate"]]
    ase = compute_ase(counts, meta)
    ase.to_csv(f"{a.out_prefix}.gene_stage_ASE.tsv.gz", sep="\t", index=False)
    cls = temporal_classes(ase)
    cls.to_csv(f"{a.out_prefix}.gene_temporal_class.tsv", sep="\t", index=False)
    print(cls.groupby(["analysis", "temporal_class"]).size().to_string())


if __name__ == "__main__":
    main()
