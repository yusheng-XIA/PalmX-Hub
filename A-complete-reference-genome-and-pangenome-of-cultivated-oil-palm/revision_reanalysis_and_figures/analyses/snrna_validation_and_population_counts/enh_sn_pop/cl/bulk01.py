#!/usr/bin/env python3
"""Score snRNA C6/C9/C10 marker programmes in the 114 bulk RNA-seq libraries (FL, TN; 19 stages x 3 replicates).
Normalisation as Fig. 2e (DESeq2 median-of-ratios on genes with >=10 counts in >=3 libraries), log2(x + 1).
Scores per library: (i) mean z-score of marker genes (z across the 114 libraries); (ii) rank score (mean within-library
percentile rank of marker genes among expressed genes; singscore-like).
Tests: FL vs TN at 95, 125 and 185 d (Welch t, n = 3 vs 3; BH over the 9 cluster x stage tests of each score/marker set);
genotype x stage interaction (185 d vs 95 d) by OLS on the 12 libraries."""
import numpy as np, pandas as pd
from scipy import stats
import statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests
B = "${ANALYSIS_DIR}"
RUN = B + "/22_answer_reviews/00_ms/03_V3/03_figure3/05_multiomics_integration/runs/RUN-MULTIOMICS-INTEGRATION-20260721-001/outputs"
W = "${CLUSTER_WORK}/enh_sn_pop"
O = W + "/bulk/"
STAGES = ["0d","15d","35d","50d","65d","80d","95d","110d","125d","140d","155d","170d","185d","12h","24h","36h","48h","60h","72h"]
cw = pd.read_csv(RUN + "/stage1_identity/sample_crosswalk_114.tsv", sep="\t", dtype=str)
counts = pd.read_csv(B + "/22_answer_reviews/00_ms/03_V3/03_figure3/00_minipan/03_rnaseq_mapping/count_matrices/pangraphrna_hisat2_graph.gene_counts.tsv", sep="\t", index_col=0)
S = cw.rna_sample.tolist(); Cn = counts[S].astype(float)
keep = (Cn >= 10).sum(axis=1) >= 3; Ck = Cn[keep]
lg = np.log(Ck.values); lgm = lg.mean(1); ok = np.isfinite(lgm)
sf = np.array([np.exp(np.median((lg[:, j] - lgm)[ok & (Ck.values[:, j] > 0)])) for j in range(lg.shape[1])])
N = np.log2(Ck.div(sf, axis=1) + 1)
Z = N.sub(N.mean(1), axis=0).div(N.std(1, ddof=1).replace(0, np.nan), axis=0)
R = N.rank(axis=0, pct=True)
info = cw.set_index("rna_sample")[["genotype", "stage"]]
print("genes kept", keep.sum(), "libraries", len(S))
# marker sets
up = pd.read_csv(W + "/sn/C6_C9_C10_top120_markers_with_annotation.tsv", sep="\t")
sets = {}
for k in ("C6", "C9", "C10"):
    g = up[up.cluster == k].sort_values("marker_rank_score", ascending=False).gene_id.tolist()
    sets[(k, "published_top50")] = g[:50]
    sets[(k, "published_top120")] = g[:120]
import os
wl = pd.read_csv(W + "/sn/out/marker_sets_top200.tsv", sep="\t") if os.path.exists(W + "/sn/out/marker_sets_top200.tsv") else pd.DataFrame(columns=["set"])
for s_, d in wl.groupby("set"):
    k = s_.split("_")[0]
    if "within" in s_: sets[(k, "within_library_top50")] = d.sort_values("rank").gene_id.tolist()[:50]
rows = []; per = []
for (k, nm), g in sets.items():
    gi = [x for x in g if x in N.index]
    z = Z.loc[gi].mean(0); r = R.loc[gi].mean(0)
    for j in S:
        per.append(dict(cluster=k, marker_set=nm, n_genes_used=len(gi), library=j, genotype=info.loc[j, "genotype"],
                        stage=info.loc[j, "stage"], mean_z=z[j], rank_score=r[j]))
per = pd.DataFrame(per)
per.to_csv(O + "bulk_scores_per_library.tsv", sep="\t", index=False)
for (k, nm), d in per.groupby(["cluster", "marker_set"]):
    for sc in ("mean_z", "rank_score"):
        for st in ("95d", "125d", "185d"):
            a = d[(d.genotype == "FL") & (d.stage == st)][sc].values; b = d[(d.genotype == "TN") & (d.stage == st)][sc].values
            assert len(a) == 3 and len(b) == 3
            t, p = stats.ttest_ind(a, b, equal_var=False)
            rows.append(dict(cluster=k, marker_set=nm, score=sc, stage=st, FL_mean=a.mean(), FL_sd=a.std(ddof=1), TN_mean=b.mean(),
                             TN_sd=b.std(ddof=1), diff_FL_minus_TN=a.mean() - b.mean(), t=t, P=p, n_genes=d.n_genes_used.iloc[0]))
res = pd.DataFrame(rows)
res["BH"] = np.nan
for (nm, sc), idx in res.groupby(["marker_set", "score"]).groups.items():
    res.loc[idx, "BH"] = multipletests(res.loc[idx, "P"], method="fdr_bh")[1]
# interaction: (FL-TN at 185 d) - (FL-TN at 95 d)
inter = []
for (k, nm), d in per.groupby(["cluster", "marker_set"]):
    for sc in ("mean_z", "rank_score"):
        dd = d[d.stage.isin(["95d", "185d"])].copy(); dd["y"] = dd[sc]
        f = smf.ols("y ~ C(genotype, Treatment('TN')) * C(stage, Treatment('95d'))", data=dd).fit()
        key = [x for x in f.params.index if ":" in x][0]
        a185 = dd[(dd.genotype == "FL") & (dd.stage == "185d")].y.values; a95 = dd[(dd.genotype == "FL") & (dd.stage == "95d")].y.values
        t2, p2 = stats.ttest_ind(a185, a95, equal_var=False)
        inter.append(dict(cluster=k, marker_set=nm, score=sc, interaction_estimate=f.params[key], interaction_P=f.pvalues[key],
                          FL_185_minus_95=a185.mean() - a95.mean(), FL_185_vs_95_P=p2))
inter = pd.DataFrame(inter)
res.to_csv(O + "bulk_FL_vs_TN_tests.tsv", sep="\t", index=False); inter.to_csv(O + "bulk_interaction_tests.tsv", sep="\t", index=False)
pd.set_option("display.width", 250)
print(res.to_string()); print(inter.to_string())
# 19-stage means for plotting
per.groupby(["cluster", "marker_set", "genotype", "stage"])[["mean_z", "rank_score"]].agg(["mean", "std"]).to_csv(O + "bulk_scores_19stage_summary.tsv", sep="\t")
