#!/usr/bin/env python3
"""Fig. 2d: FL vs TN statistics for the four composite molecular axes (n = 3 per material and stage).

Input  data/molecular_phenotype_scores_by_sample.tsv  (stage26 per-sample scores; identical to the trace
       recomputation, max |diff| <= 1.7e-15; per-stage means = Source Data Fig2d_axis_scores)
Tests
  1. per stage: two-sided Welch t test FL vs TN; BH within each axis over 19 stages (q_axis) and over all
     76 axis-stage tests (q_all)
  2. per stage: FL - TN contrast from the cell-means model score ~ material * stage fitted per axis
     (pooled residual variance, 76 df); BH within axis
  3. per axis: two-way ANOVA score ~ material * stage (sum-to-zero contrasts, type III): material main effect and
     material x stage interaction
"""
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats
import statsmodels.formula.api as smf
from statsmodels.stats.anova import anova_lm
from statsmodels.stats.multitest import multipletests

HERE = Path(__file__).resolve().parent
OUT = HERE / "out"; OUT.mkdir(exist_ok=True)
AX = {"P02": "storage lipids", "P01": "oleic balance", "P03": "hydrolytic deterioration",
      "P04": "oxidative deterioration"}
STAGES = ["0d", "15d", "35d", "50d", "65d", "80d", "95d", "110d", "125d", "140d", "155d", "170d", "185d",
          "12h", "24h", "36h", "48h", "60h", "72h"]
d = pd.read_csv(HERE / "data/molecular_phenotype_scores_by_sample.tsv", sep="\t")
d = d[d.axis_id.isin(AX)].rename(columns={"molecular_phenotype_score": "score", "genotype": "material"})
assert len(d) == 4 * 114

# cross-check against Source Data means
sd = pd.read_excel(HERE.parents[1] / "deliver/Source_Data/Source_Data.xlsx", sheet_name="Fig2d_axis_scores")
m = d.groupby(["axis_id", "material", "stage"]).score.mean().reset_index()
chk = m.merge(sd, left_on=["axis_id", "material", "stage"], right_on=["axis_id", "genotype", "stage"])
assert len(chk) == 152
print("max |mean - SourceData| =", (chk.score - chk["mean"]).abs().max())

rows, anova_rows = [], []
for ax, name in AX.items():
    a = d[d.axis_id == ax].copy()
    a["stage"] = pd.Categorical(a.stage, STAGES)
    fit = smf.ols("score ~ C(material, Sum) * C(stage, Sum)", data=a).fit()
    an = anova_lm(fit, typ=3)
    for term, lab in (("C(material, Sum)", "material"), ("C(stage, Sum)", "stage"),
                      ("C(material, Sum):C(stage, Sum)", "material x stage")):
        anova_rows.append(dict(axis_id=ax, axis=name, term=lab, df=int(an.loc[term, "df"]),
                               df_resid=int(fit.df_resid), F=an.loc[term, "F"], P=an.loc[term, "PR(>F)"]))
    # heteroscedasticity-robust (HC3) Wald F tests of the same terms (sensitivity)
    rob = smf.ols("score ~ C(material, Sum) * C(stage, Sum)", data=a).fit(cov_type="HC3")
    names = rob.model.exog_names
    for lab, pick in (("material", lambda n: n.startswith("C(material") and ":" not in n),
                      ("material x stage", lambda n: ":" in n)):
        idx = [i for i, n in enumerate(names) if pick(n)]
        R = np.zeros((len(idx), len(names)))
        for k, i in enumerate(idx):
            R[k, i] = 1
        wt = rob.wald_test(R, use_f=True, scalar=True)
        anova_rows.append(dict(axis_id=ax, axis=name, term=lab + " (HC3-robust)", df=len(idx),
                               df_resid=int(fit.df_resid), F=float(wt.statistic), P=float(wt.pvalue)))
    s2 = fit.mse_resid; dfr = fit.df_resid
    for s in STAGES:
        fl = a[(a.material == "FL") & (a.stage == s)].score.to_numpy()
        tn = a[(a.material == "TN") & (a.stage == s)].score.to_numpy()
        assert len(fl) == 3 and len(tn) == 3
        w = stats.ttest_ind(fl, tn, equal_var=False)
        diff = fl.mean() - tn.mean()
        tp = diff / np.sqrt(s2 * (1 / 3 + 1 / 3))
        pp = 2 * stats.t.sf(abs(tp), dfr)
        vf, vt = fl.var(ddof=1) / 3, tn.var(ddof=1) / 3
        wdf = (vf + vt) ** 2 / (vf ** 2 / 2 + vt ** 2 / 2)
        rows.append(dict(axis_id=ax, axis=name, stage=s, phase="post-harvest" if s.endswith("h") else "development",
                         FL_mean=fl.mean(), TN_mean=tn.mean(), FL_sd=fl.std(ddof=1), TN_sd=tn.std(ddof=1),
                         diff_FL_minus_TN=diff, higher="FL" if diff > 0 else "TN",
                         welch_t=w.statistic, welch_df=wdf, welch_P=w.pvalue,
                         model_t=tp, model_df=dfr, model_P=pp))
r = pd.DataFrame(rows)
for col, grp in (("welch_P", "welch"), ("model_P", "model")):
    r[f"{grp}_q_axis"] = r.groupby("axis_id")[col].transform(lambda p: multipletests(p, method="fdr_bh")[1])
    r[f"{grp}_q_all"] = multipletests(r[col], method="fdr_bh")[1]
an = pd.DataFrame(anova_rows)
# sign summary per axis
for ax in AX:
    x = r[r.axis_id == ax]
    print(ax, "FL higher at", int((x.diff_FL_minus_TN > 0).sum()), "of 19; TN higher at", int((x.diff_FL_minus_TN < 0).sum()),
          "; Welch q_axis<0.05:", x.loc[x.welch_q_axis < 0.05, "stage"].tolist())
r.to_csv(OUT / "fig2d_stage_tests.tsv", sep="\t", index=False, float_format="%.6g")
an.to_csv(OUT / "fig2d_twoway_anova.tsv", sep="\t", index=False, float_format="%.6g")

pd.set_option("display.width", 250)
print(an.to_string())
for ax in AX:
    x = r[r.axis_id == ax]
    print(f"\n== {ax} {AX[ax]}")
    print(x[["stage", "diff_FL_minus_TN", "welch_P", "welch_q_axis", "model_P", "model_q_axis"]].to_string(index=False,
          float_format=lambda v: f"{v:.3g}"))
