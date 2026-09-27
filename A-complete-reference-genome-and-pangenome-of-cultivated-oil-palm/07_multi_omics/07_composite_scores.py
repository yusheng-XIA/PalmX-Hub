#!/usr/bin/env python3
"""Composite molecular scores (Fig. 2c) from metabolite z-scores.

score(sample) = sum_i(d_i w_i z_i) / sum_i |d_i w_i|
  z_i : consensus z-score of member compound i across the 114 samples (missing -> 0, the across-sample mean)
  d_i : predefined direction on the axis (+1, -1 or 0)
  w_i : annotation-confidence weight (1 Level 2; 0.35 MS1 exact-mass lipid-class proxy supported across ionization
        modes or by adduct pairs; 0.15 proxy based on a default adduct only)
The denominator is summed over all member compounds. FL vs TN at each stage: two-sided Welch t test (n = 3 + 3),
BH across the 19 stages of each axis; per axis a two-way ANOVA of score ~ material * stage (type III SS,
sum-to-zero contrasts).

usage: 07_composite_scores.py zscores.tsv members.tsv samples.tsv OUT_PREFIX
  zscores.tsv : compound x sample consensus z-scores; members.tsv : axis, compound, direction, weight
  samples.tsv : sample, material (FL/TN), stage
"""
import sys

import numpy as np
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf
from scipy import stats
from statsmodels.stats.multitest import multipletests

z = pd.read_csv(sys.argv[1], sep="\t", index_col=0)
mem = pd.read_csv(sys.argv[2], sep="\t")
meta = pd.read_csv(sys.argv[3], sep="\t", index_col=0)
out = sys.argv[4]

scores = {}
for axis, m in mem.groupby("axis"):
    zz = z.reindex(m.compound).fillna(0.0).to_numpy()
    dw = (m.direction * m.weight).to_numpy()
    scores[axis] = pd.Series(dw @ zz / np.abs(dw).sum(), index=z.columns)
S = pd.DataFrame(scores).join(meta)
S.to_csv(f"{out}.scores.tsv", sep="\t")

rows, anova = [], []
for axis in scores:
    for st, d in S.groupby("stage"):
        a, b = d.loc[d.material == "FL", axis], d.loc[d.material == "TN", axis]
        rows.append(dict(axis=axis, stage=st, mean_FL=a.mean(), mean_TN=b.mean(),
                         P=stats.ttest_ind(a, b, equal_var=False).pvalue))
    fit = smf.ols(f"Q('{axis}') ~ C(material, Sum) * C(stage, Sum)", data=S).fit()
    t = sm.stats.anova_lm(fit, typ=3); t["axis"] = axis; anova.append(t)
R = pd.DataFrame(rows)
R["P_adj"] = R.groupby("axis").P.transform(lambda p: multipletests(p, method="fdr_bh")[1])
R.to_csv(f"{out}.stage_tests.tsv", sep="\t", index=False)
pd.concat(anova).to_csv(f"{out}.anova_typeIII.tsv", sep="\t")
