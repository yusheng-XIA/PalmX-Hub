#!/usr/bin/env python3
"""FAD2 qRT-PCR (EgFAD2.1, three primer pairs; ACTIN reference).

Per stage and material: wells without a Ct are excluded and the remaining wells averaged (four wells per primer
pair on the same plate as the ACTIN wells). log2(FAD2/ACTIN) = -dCt, dCt = mean FAD2 Ct - mean ACTIN Ct; -dCt is
averaged across the three primer pairs. FL is compared with TN, NS and TK with two-sided exact Wilcoxon
signed-rank tests paired by stage (13 developmental and 6 post-harvest pairs), BH across the three comparisons
within each phase.

usage: 08_qpcr_fad2.py ct_wells.tsv OUT_PREFIX
  ct_wells.tsv: material, stage, phase (development/postharvest), target (FAD2_P1/FAD2_P2/FAD2_P3/ACTIN), Ct
"""
import sys

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

ct = pd.read_csv(sys.argv[1], sep="\t")
out = sys.argv[2]
ct["Ct"] = pd.to_numeric(ct.Ct, errors="coerce")
m = ct.dropna(subset=["Ct"]).groupby(["material", "stage", "phase", "target"]).Ct.mean().unstack("target")
primers = [c for c in m.columns if c.startswith("FAD2")]
dct = pd.DataFrame({p: -(m[p] - m["ACTIN"]) for p in primers})
m["neg_dCt"] = dct.mean(axis=1)
m.reset_index().to_csv(f"{out}.neg_dCt.tsv", sep="\t", index=False)

v = m["neg_dCt"].unstack("material")
rows = []
for phase, d in v.groupby(level="phase"):
    for other in ["TN", "NS", "TK"]:
        x = d[["FL", other]].dropna()
        res = stats.wilcoxon(x.FL, x[other], alternative="two-sided", method="exact")
        rows.append(dict(phase=phase, comparison=f"FL_vs_{other}", n_pairs=len(x),
                         median_diff=float(np.median(x.FL - x[other])), P=res.pvalue))
R = pd.DataFrame(rows)
R["P_adj"] = R.groupby("phase").P.transform(lambda p: multipletests(p, method="fdr_bh")[1])
R.to_csv(f"{out}.wilcoxon.tsv", sep="\t", index=False)
print(R.to_string(index=False))
