#!/usr/bin/env python3
"""Full-scan QQ and Manhattan display tables for Supplementary Fig. 9.

Input: full EMMAX SV .ps files (370,706 tests each; columns id, beta, se, P), copied from
03_V3/05_figure/07_GWAS_zero_extreme_excluded_20260723/sv/results/yield/<trait>/<trait>.ps
QQ: expected = -log10(i/N) over all N tests (same convention as ED Fig. 6b); every test with
P < 1e-3 kept at its exact rank; the bulk is thinned (every 20th rank) for display.
Manhattan: existing Source Data rule (every 13th test in file order + all tests above the
Bonferroni threshold); P stored as scientific-notation strings.
"""
import gzip
from pathlib import Path
from statistics import NormalDist

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
N_EXP = 370706
BONF = 0.05 / N_EXP
out = {}
for trait, qqname, mname, lead in [("Shell_weight_g", "SF9b_shell_weight_QQ.tsv", "SF9a_manhattan.tsv",
                                    "SV_chr01B_7154140_2e12f281"),
                                   ("Flesh_thickness_mm", "SF9e_flesh_thickness_QQ.tsv", "SF9d_manhattan.tsv",
                                    "SV_chr01B_6487575_305019c5")]:
    d = pd.read_csv(HERE / f"full/{trait}.ps.gz", sep="\t", header=None, names=["variant_id", "beta", "se", "p"],
                    dtype={"p": str})
    d["p_float"] = d.p.astype(float)
    N = len(d)
    assert N == N_EXP, N
    lam = NormalDist().inv_cdf(1 - np.median(d.p_float) / 2) ** 2 / 0.454936
    # QQ
    p = np.sort(d.p_float.values)
    rank = np.arange(1, N + 1)
    keep = (p < 1e-3) | ((rank - 1) % 20 == 0) | (rank == N)
    qq = pd.DataFrame({"rank": rank[keep], "expected_neglog10P": -np.log10(rank[keep] / N),
                       "observed_neglog10P": -np.log10(np.clip(p[keep], 1e-300, None)),
                       "display": np.where(p[keep] < 1e-3, "exact_tail", "thinned_every_20th")})
    qq["full_test_count"] = N
    qq["lambda_GC_full"] = round(lam, 6)
    qq.to_csv(HERE / qqname, sep="\t", index=False, float_format="%.6g")
    # Manhattan
    parts = d.variant_id.str.extract(r"SV_(chr\d+B)_(\d+)_")
    d["chrom"], d["pos"] = parts[0], parts[1].astype(int)
    d["neglog10_p"] = -np.log10(d.p_float)
    d["significant"] = d.p_float < BONF
    idx = np.arange(N)
    sel = (idx % 13 == 0) | d.significant.values
    mm = d.loc[sel, ["chrom", "pos", "variant_id", "p", "neglog10_p", "beta", "se", "significant"]].copy()
    mm["p"] = mm.p.map(lambda s: f"{float(s):.3e}")
    mm.to_csv(HERE / mname, sep="\t", index=False, float_format="%.6g")
    out[trait] = dict(N=N, lam=lam, n_sig=int(d.significant.sum()), min_p=p[0], top=d.loc[d.p_float.idxmin(), "variant_id"],
                      lead_p=float(d.loc[d.variant_id == lead, "p_float"].iloc[0]) if (d.variant_id == lead).any() else None,
                      qq_rows=len(qq), qq_tail=int((qq.display == "exact_tail").sum()), max_exp=qq.expected_neglog10P.max(),
                      max_obs=qq.observed_neglog10P.max(), manh_rows=len(mm))
for k, v in out.items():
    print(k, v)
