#!/usr/bin/env python3
"""Fig. 3f robustness: pairwise effect sizes and equivalence (TOST) for Ka/Ks and Ks among the
four temporal ASE classes (NoDiff, HapDom, Sub, NoASE) in FL and TN.

Inputs (read-only):
  per-gene Ka/Ks + class : work/trace/Fig3A/out/recomputed_fig3f_kaks.tsv (s06_kaks.py; = Fig. 3f data)
  genome-wide 1:1 pairs  : raw KaKs_Calculator tables (all allelic 1:1 pairs, same filters)
Margins:
  (1) absolute margin = 10% of the genome-wide median of all allelic 1:1 pairs passing the Fig. 3f filters
      (0 < Ks <= 0.10, 0 <= Ka/Ks <= 3) in the same hybrid;
  (2) standardized margin: Cliff's delta within +/-0.147 ("negligible", Romano et al. 2006).
Tests:
  HL shift (Hodges-Lehmann) with distribution-free (Moses) 90% and 95% CIs; bootstrap 95% CI of the difference
  in medians (2,000 resamples); Wilcoxon-Mann-Whitney TOST (two one-sided shifted tests), P_TOST = max of the two;
  Cliff's delta with bootstrap 90% CI; smallest symmetric margin compatible with equivalence at alpha = 0.05
  (= max |bound| of the 90% Moses CI).
"""
import numpy as np, pandas as pd
from scipy import stats

A = "${ANALYSIS_DIR}/22_answer_reviews/00_ms/03_V3"
RAW = {"FL": (A + "/03_figure3/03_3c/final.Fig3c_American_hap1_Africa_hap2_KaKs_all.tsv", "africa_gene"),
       "TN": (A + "/02_figure/_runs/RUN-ASE-FIGURES-CURRENT-V3-001/work/tn_kaks/TN_Dura_Pisifera_KaKs_all_current.tsv", "gene_dura")}
IN = "${CLUSTER_WORK}/trace/Fig3A/out/recomputed_fig3f_kaks.tsv"
OUT = "${CLUSTER_WORK}/enh_3fgj/out"
CL = ["NoDiff", "HapDom", "Sub", "NoASE"]
rng = np.random.default_rng(20260924)
log = open(OUT + "/s01_log.txt", "w")
def P(*a):
    s = " ".join(str(x) for x in a); print(s, flush=True); log.write(s + "\n"); log.flush()

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

d = pd.read_csv(IN, sep="\t")
rows, bg_rows = [], []
for an in ["FL", "TN"]:
    path, col = RAW[an]
    raw = pd.read_csv(path, sep="\t", low_memory=False)
    g = pd.DataFrame({"gene_id": raw[col].astype(str).str.replace("evm.model.", "evm.TU.", regex=False),
                      "Ka": pd.to_numeric(raw.Ka, errors="coerce"), "Ks": pd.to_numeric(raw.Ks, errors="coerce"),
                      "KaKs": pd.to_numeric(raw["Ka/Ks"], errors="coerce")})
    g = g[g.Ka.notna() & (g.Ks > 0) & (g.Ks <= 0.10) & (g.KaKs >= 0) & (g.KaKs <= 3)]
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
res.to_csv(OUT + "/fig3f_pairwise_equivalence.tsv", sep="\t", index=False, float_format="%.6g")
pd.DataFrame(bg_rows).to_csv(OUT + "/fig3f_genome_background.tsv", sep="\t", index=False, float_format="%.6g")
for (an, met), s in res.groupby(["hybrid", "metric"]):
    P(f"SUMMARY {an} {met}: max|HL|={s.HL_shift.abs().max():.4g}; max TOST P={s.TOST_P.max():.3g}; "
      f"max smallest-margin={s.smallest_equivalence_margin.max():.4g}; max|cliff|={s.cliffs_delta.abs().max():.4f}; "
      f"all cliff-equivalent={s.cliffs_equivalent_0147.all()}")
    s2 = s[(s.class_1 != "NoASE") & (s.class_2 != "NoASE")]
    P(f"   excl. NoASE: max TOST P={s2.TOST_P.max():.3g}; max smallest-margin={s2.smallest_equivalence_margin.max():.4g}")
log.close()
