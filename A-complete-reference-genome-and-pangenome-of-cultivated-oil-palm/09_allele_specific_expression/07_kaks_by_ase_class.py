#!/usr/bin/env python3
"""Ka/Ks and Ks of allelic pairs among the four temporal ASE classes (Fig. 3f; Extended Data Fig. 4).

Gene pairs with 0 < Ks <= 0.10 and 0 <= Ka/Ks <= 3 are retained. Classes are compared with Kruskal-Wallis
tests and pairwise two-sided Mann-Whitney U tests (Holm correction). Pairwise differences are estimated as
Hodges-Lehmann shifts with distribution-free (Moses) confidence intervals and tested for equivalence with two
one-sided Mann-Whitney tests (TOST), using a margin of +/-10% of the median of all one-to-one allelic pairs of the
same hybrid passing the filters; Cliff's delta is reported as a standardized effect size.

usage: 07_kaks_by_ase_class.py --kaks FL=FL.KaKs_all.tsv --kaks TN=TN.KaKs_all.tsv \
           --classes gene_temporal_class.tsv --out-prefix fig3f
  KaKs tables from 06_allelic_kaks.py (gene_a = backbone gene id used in the ASE analysis);
  classes from 04_ase_test.py (analysis, gene_id, temporal_class).
"""
import argparse
import numpy as np, pandas as pd
from scipy import stats

CL = ["NoDiff", "HapDom", "Sub", "NoASE"]
rng = np.random.default_rng(20260924)


def P(*a):
    print(" ".join(str(x) for x in a), flush=True)


def hl_moses(x, y, conf):
    d = np.subtract.outer(x, y).ravel()
    n = d.size; n1, n2 = len(x), len(y)
    z = stats.norm.ppf(1 - (1 - conf) / 2)
    k = int(np.floor(n1 * n2 / 2 - z * np.sqrt(n1 * n2 * (n1 + n2 + 1) / 12)))
    k = max(k, 1)
    lo, hi = np.partition(d, [k - 1, n - k])[[k - 1, n - k]]
    return lo, hi, d

def cliff(x, y):
    u = stats.mannwhitneyu(x, y, alternative="two-sided").statistic
    return 2 * u / (len(x) * len(y)) - 1

def holm(p):
    p = np.asarray(p, float); o = np.argsort(p); m = len(p)
    adj = np.maximum.accumulate((m - np.arange(m)) * p[o]); out = np.empty(m); out[o] = np.minimum(adj, 1); return out


ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("--kaks", action="append", required=True, help="HYBRID=path")
ap.add_argument("--classes", required=True)
ap.add_argument("--out-prefix", required=True)
args = ap.parse_args()
cls = pd.read_csv(args.classes, sep="\t").rename(columns={"temporal_class": "overall_class"})
rows, bg_rows, kw_rows = [], [], []
for spec in args.kaks:
    an, path = spec.split("=", 1)
    raw = pd.read_csv(path, sep="\t", low_memory=False)
    g = pd.DataFrame({"gene_id": raw["gene_a"].astype(str).str.replace("evm.model.", "evm.TU.", regex=False),
                      "Ka": pd.to_numeric(raw.Ka, errors="coerce"), "Ks": pd.to_numeric(raw.Ks, errors="coerce"),
                      "KaKs": pd.to_numeric(raw["Ka/Ks"], errors="coerce")})
    g = g[g.Ka.notna() & (g.Ks > 0) & (g.Ks <= 0.10) & (g.KaKs >= 0) & (g.KaKs <= 3)]
    d = g.merge(cls[cls.analysis == an][["gene_id", "overall_class"]], on="gene_id").assign(analysis=an)
    for met in ["KaKs", "Ks"]:
        grp = [d.loc[d.overall_class == c, met] for c in CL]
        pw = [(CL[i], CL[j], stats.mannwhitneyu(grp[i], grp[j], alternative="two-sided").pvalue)
              for i in range(4) for j in range(i + 1, 4) if len(grp[i]) and len(grp[j])]
        ph = holm([p for *_, p in pw])
        kw_rows.append(dict(hybrid=an, metric=met, kruskal_P=stats.kruskal(*[x for x in grp if len(x)]).pvalue,
                            **{f"{a}_vs_{b}_P_holm": q for (a, b, _), q in zip(pw, ph)}))
    x = d[d.analysis == an]
    nocls = g[~g.gene_id.isin(x.gene_id)]
    for met in ["KaKs", "Ks"]:
        med_all = float(g[met].median()); med_out = float(nocls[met].median())
        q1, q3 = np.quantile(g[met], [.25, .75])
        margin = 0.10 * med_all
        bg_rows.append(dict(hybrid=an, metric=met, n_all_pairs=len(g), median_all_pairs=med_all, q1=q1, q3=q3,
                            n_pairs_without_ASE_class=len(nocls), median_pairs_without_ASE_class=med_out,
                            margin_10pct_of_genome_median=margin))
        P(f"{an} {met}: genome-wide 1:1 pairs n={len(g)} median={med_all:.5g} IQR={q1:.5g}-{q3:.5g}; "
          f"pairs without ASE class n={len(nocls)} median={med_out:.5g}; margin={margin:.5g}")
        grp = {c: x.loc[x.overall_class == c, met].to_numpy(float) for c in CL}
        for i in range(4):
            for j in range(i + 1, 4):
                a, b = grp[CL[i]], grp[CL[j]]
                if len(a) == 0 or len(b) == 0:
                    continue
                lo90, hi90, dd = hl_moses(a, b, 0.90)
                hl = float(np.median(dd))
                lo95, hi95, _ = hl_moses(a, b, 0.95)
                del dd
                bm = np.array([np.median(rng.choice(a, len(a))) - np.median(rng.choice(b, len(b))) for _ in range(2000)])
                p_lo = stats.mannwhitneyu(a + margin, b, alternative="greater").pvalue
                p_hi = stats.mannwhitneyu(a - margin, b, alternative="less").pvalue
                cd = cliff(a, b)
                cb = []
                for _ in range(500):
                    cb.append(cliff(rng.choice(a, len(a)), rng.choice(b, len(b))))
                c90 = np.percentile(cb, [5, 95])
                # TOST on Cliff's delta (+/-0.147) via bootstrap percentile CI (90%)
                rows.append(dict(hybrid=an, metric=met, class_1=CL[i], class_2=CL[j], n_1=len(a), n_2=len(b),
                                 median_1=np.median(a), median_2=np.median(b),
                                 median_difference=np.median(a) - np.median(b),
                                 median_difference_boot95_low=np.percentile(bm, 2.5),
                                 median_difference_boot95_high=np.percentile(bm, 97.5),
                                 HL_shift=hl, HL_CI95_low=lo95, HL_CI95_high=hi95, HL_CI90_low=lo90, HL_CI90_high=hi90,
                                 margin=margin, TOST_P=max(p_lo, p_hi),
                                 smallest_equivalence_margin=max(abs(lo90), abs(hi90)),
                                 MWU_P_two_sided=stats.mannwhitneyu(a, b, alternative="two-sided").pvalue,
                                 cliffs_delta=cd, cliffs_delta_CI90_low=c90[0], cliffs_delta_CI90_high=c90[1],
                                 cliffs_equivalent_0147=bool(c90[0] > -0.147 and c90[1] < 0.147)))
                r = rows[-1]
                P(f"  {an} {met} {CL[i]} vs {CL[j]}: HL={hl:.4g} 95%CI[{lo95:.4g},{hi95:.4g}] 90%CI[{lo90:.4g},{hi90:.4g}] "
                  f"TOST P={r['TOST_P']:.3g} (margin {margin:.4g}); dmed={r['median_difference']:.4g} "
                  f"boot95[{r['median_difference_boot95_low']:.4g},{r['median_difference_boot95_high']:.4g}] "
                  f"cliff={cd:.4f} 90%[{c90[0]:.4f},{c90[1]:.4f}]")
res = pd.DataFrame(rows)
pd.DataFrame(kw_rows).to_csv(args.out_prefix + ".kruskal_holm.tsv", sep="\t", index=False, float_format="%.6g")
res.to_csv(args.out_prefix + ".pairwise_equivalence.tsv", sep="\t", index=False, float_format="%.6g")
pd.DataFrame(bg_rows).to_csv(args.out_prefix + ".genome_background.tsv", sep="\t", index=False, float_format="%.6g")
for (an, met), s in res.groupby(["hybrid", "metric"]):
    P(f"SUMMARY {an} {met}: max|HL|={s.HL_shift.abs().max():.4g}; max TOST P={s.TOST_P.max():.3g}; "
      f"max smallest-margin={s.smallest_equivalence_margin.max():.4g}; max|cliff|={s.cliffs_delta.abs().max():.4f}; "
      f"all cliff-equivalent={s.cliffs_equivalent_0147.all()}")
    s2 = s[(s.class_1 != "NoASE") & (s.class_2 != "NoASE")]
    P(f"   excl. NoASE: max TOST P={s2.TOST_P.max():.3g}; max smallest-margin={s2.smallest_equivalence_margin.max():.4g}")
