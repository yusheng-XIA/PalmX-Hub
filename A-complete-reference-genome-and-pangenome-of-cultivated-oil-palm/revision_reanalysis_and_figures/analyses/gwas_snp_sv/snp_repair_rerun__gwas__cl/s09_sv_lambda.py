"""lambda_GC of the published SV scans over the 370,136 unique SVs (first occurrence of duplicated IDs), 60 traits."""
import numpy as np, pandas as pd
from scipy.stats import chi2
from config import *
tasks = pd.read_csv(RUN / "manifests/sv_tasks.tsv", sep="\t"); rows = []
for c, t in zip(tasks.category, tasks.trait):
    ps = pd.read_csv(RUN / f"sv/results/{c}/{t}/{t}.ps", sep="\t", header=None, usecols=[0, 3]).drop_duplicates(0)
    rows.append(dict(trait=t, n_sv=len(ps), lambda_gc_sv=float(chi2.isf(np.median(ps[3]), 1) / chi2.ppf(.5, 1)), n_sig_sv=int((ps[3] < BONF_SV).sum())))
pd.DataFrame(rows).to_csv(W / "out/sv_lambda_unique.tsv", sep="\t", index=False); print(pd.DataFrame(rows).describe())
