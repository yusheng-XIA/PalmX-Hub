#!/usr/bin/env python3
"""snRNA C6/C9/C10 checks: Scrublet doublets per library, clustering stability on the Seurat SNN graph
(Louvain, resolution 0.4-1.2, plus 80% subsampling), within-library marker sets for bulk scoring."""
import sys, gzip, json
from pathlib import Path
import numpy as np, pandas as pd, scipy.sparse as sp, scipy.io
W = Path("${CLUSTER_WORK}/enh_sn_pop/sn"); X = W / "export"; O = W / "out"; O.mkdir(exist_ok=True)
log = open(O / "sn02_log.txt", "w")
def P(*a):
    print(*a, flush=True); print(*a, file=log, flush=True)
meta = pd.read_csv(X / "meta.tsv", sep="\t", index_col=0)
cells = open(X / "cells.txt").read().split(); genes = open(X / "genes.txt").read().split()
assert cells == meta.index.tolist()
C = scipy.io.mmread(X / "counts.mtx").tocsc()        # genes x cells
P("counts", C.shape, "nnz", C.nnz)
assert (np.asarray(C.sum(0)).ravel() == meta.nCount_RNA.values).all()
cl = meta.seurat_clusters.astype(int).values
lib = meta["orig.ident"].values
LAB = {"E_1_95_2": "FL 95 d", "E_1_125_4": "FL 125 d", "E_1_185_9": "FL 185 d", "E_4_95_5": "TN 95 d", "E_4_125_4": "TN 125 d", "E_4_185_5": "TN 185 d"}
P("seurat_clusters == RNA_snn_res.1:", (meta.seurat_clusters == meta["RNA_snn_res.1"]).all(), "; == res.0.8:", (meta.seurat_clusters == meta["RNA_snn_res.0.8"]).mean())

# ---------------- 1. Scrublet per library (scores cached); calls = top expected-rate fraction of each library
import os
cache = O / "scrublet_per_nucleus.tsv.gz"
if cache.exists():
    ds = pd.read_csv(cache, sep="\t", index_col=0).loc[cells, "dbl_score"].values; auto = None
else:
    import scrublet as scr
    Ct = C.T.tocsr(); ds = np.zeros(len(cells))
    for L in sorted(set(lib)):
        idx = np.where(lib == L)[0]
        s = scr.Scrublet(Ct[idx], expected_doublet_rate=0.008 * len(idx) / 1000, random_state=10086)
        sc, pr = s.scrub_doublets(min_counts=2, min_cells=3, min_gene_variability_pctl=85, n_prin_comps=30, verbose=False)
        ds[idx] = sc
dp = np.zeros(len(cells), bool); rows = []
for L in sorted(set(lib)):
    idx = np.where(lib == L)[0]; n = len(idx); rate = 0.008 * n / 1000
    k = int(round(rate * n)); top = idx[np.argsort(-ds[idx])[:k]]; dp[top] = True
    rows.append(dict(library=LAB[L], nuclei=n, expected_rate=round(rate, 4), called_doublets=k,
                     score_cutoff=round(float(ds[top].min()), 4), median_score=round(float(np.median(ds[idx])), 4)))
pd.DataFrame(rows).to_csv(O / "scrublet_by_library.tsv", sep="\t", index=False)
meta["dbl_score"] = ds; meta["dbl_pred"] = dp
# score percentile within library
meta["dbl_pctl"] = meta.groupby("orig.ident").dbl_score.rank(pct=True)
per = meta.groupby("seurat_clusters").agg(nuclei=("dbl_pred", "size"), called_doublets=("dbl_pred", "sum"),
                                          median_score=("dbl_score", "median"), median_score_pctl=("dbl_pctl", "median"),
                                          median_UMI=("nCount_RNA", "median"), median_genes=("nFeature_RNA", "median"))
per["doublet_pct"] = per.called_doublets / per.nuclei * 100
per.to_csv(O / "scrublet_by_cluster.tsv", sep="\t")
cmp = []
for k, L in ((9, "E_1_185_9"), (9, "E_4_125_4"), (10, "E_4_185_5"), (6, "E_4_185_5")):
    a = meta[(meta.seurat_clusters == k) & (meta["orig.ident"] == L)]; b = meta[(meta.seurat_clusters != k) & (meta["orig.ident"] == L)]
    cmp.append(dict(cluster=f"C{k}", library=LAB[L], n_in=len(a), n_other=len(b), doublet_pct_in=a.dbl_pred.mean() * 100,
                    doublet_pct_other=b.dbl_pred.mean() * 100, median_score_in=a.dbl_score.median(), median_score_other=b.dbl_score.median(),
                    median_pctl_in=a.dbl_pctl.median(), median_UMI_in=a.nCount_RNA.median(), median_UMI_other=b.nCount_RNA.median(),
                    median_genes_in=a.nFeature_RNA.median(), median_genes_other=b.nFeature_RNA.median()))
pd.DataFrame(cmp).to_csv(O / "scrublet_same_library.tsv", sep="\t", index=False)
P(pd.DataFrame(rows).to_string()); P(per.to_string()); P(pd.DataFrame(cmp).to_string())
meta[["orig.ident", "seurat_clusters", "dbl_score", "dbl_pctl", "dbl_pred"]].to_csv(cache, sep="\t", compression="gzip")

# ---------------- 2. Clustering stability on the Seurat SNN graph
import igraph as ig, random
from sklearn.metrics import adjusted_rand_score as ARI
T = pd.read_csv(X / "graph_RNA_snn.tsv.gz", sep="\t")
gn = open(X / "graph_RNA_snn_names.txt").read().split(); assert gn == cells
T = T[T.i < T.j]
G = ig.Graph(n=len(cells), edges=list(zip(T.i.values - 1, T.j.values - 1)), directed=False)
G.es["weight"] = T.x.values
P("graph edges", G.ecount())
def run(res, sub=None, seed=10086):
    g = G if sub is None else G.induced_subgraph(sub)
    random.seed(seed); ig.set_random_number_generator(random)
    best_q, best_m = None, None
    for st in range(3):   # Seurat n.start analogue: keep the partition with the highest modularity
        part = g.community_multilevel(weights="weight", resolution=res)
        q = g.modularity(part.membership, weights="weight", resolution=res)
        if best_q is None or q > best_q: best_q, best_m = q, part.membership
    return np.array(best_m)
def jac(a, b):
    return len(a & b) / len(a | b)
def best(orig_set, lab, ids):
    """best-matching new cluster for an original cluster: Jaccard, and share of the original nuclei in it"""
    s = pd.Series(lab, index=ids)
    top = s[list(orig_set)].value_counts()
    out = []
    for c in top.index[:3]:
        new = set(s.index[s.values == c]); out.append((jac(orig_set, new), len(orig_set & new) / len(orig_set), int(c), len(new)))
    return max(out)
orig = {k: set(np.where(cl == k)[0]) for k in (6, 9, 10)}
rows = []
for res in (0.4, 0.5, 0.6, 0.8, 1.0, 1.2):
    m = run(res)
    r = dict(resolution=res, n_clusters=len(set(m)), ARI_vs_published=ARI(cl, m))
    for k in (6, 9, 10):
        j, f, c, n = best(orig[k], m, np.arange(len(cells)))
        r[f"C{k}_jaccard"] = j; r[f"C{k}_share_in_best"] = f; r[f"C{k}_best_size"] = n
    rows.append(r); P(r)
# stored Seurat solutions (same graph, Seurat's own Louvain)
for col in ("RNA_snn_res.0.5", "RNA_snn_res.0.8", "RNA_snn_res.1"):
    m = meta[col].values
    r = dict(resolution=f"Seurat stored {col.split('res.')[1]}", n_clusters=len(set(m)), ARI_vs_published=ARI(cl, m))
    for k in (6, 9, 10):
        j, f, c, n = best(orig[k], m, np.arange(len(cells)))
        r[f"C{k}_jaccard"] = j; r[f"C{k}_share_in_best"] = f; r[f"C{k}_best_size"] = n
    rows.append(r); P(r)
pd.DataFrame(rows).to_csv(O / "cluster_stability_resolution.tsv", sep="\t", index=False)
# subsampling at the published resolution (1.0)
rng = np.random.default_rng(10086); rows = []
for b in range(20):
    sub = np.sort(rng.choice(len(cells), int(0.8 * len(cells)), replace=False))
    m = run(1.0, sub=sub, seed=10086 + b)
    r = dict(rep=b + 1, n_clusters=len(set(m)), ARI_vs_published=ARI(cl[sub], m))
    subset = set(sub)
    for k in (6, 9, 10):
        o = orig[k] & subset
        j, f, c, n = best(o, m, sub)
        r[f"C{k}_jaccard"] = j; r[f"C{k}_share_in_best"] = f
    rows.append(r); P(r)
pd.DataFrame(rows).to_csv(O / "cluster_stability_subsample80.tsv", sep="\t", index=False)

# ---------------- 3. Within-library marker sets (for bulk scoring)
lib_size = np.asarray(C.sum(0)).ravel()
def lognorm(M, cols):
    S = M[:, cols].multiply(1e4 / lib_size[cols]).tocsc(); S.data = np.log1p(S.data); return S
res = []
for k, L in ((9, "E_1_185_9"), (10, "E_4_185_5"), (6, "E_4_185_5"), (9, None), (10, None), (6, None)):
    cols = np.where(lib == L)[0] if L else np.arange(len(cells))
    E = lognorm(C, cols); ink = cl[cols] == k
    mi = np.asarray(E[:, ink].mean(1)).ravel(); mo = np.asarray(E[:, ~ink].mean(1)).ravel()
    pi = np.asarray((E[:, ink] > 0).mean(1)).ravel(); po = np.asarray((E[:, ~ink] > 0).mean(1)).ravel()
    lf = np.log2((np.expm1(mi) + 1e-6) / (np.expm1(mo) + 1e-6))
    score = lf * np.maximum(pi - po, 0)
    d = pd.DataFrame(dict(set=f"C{k}_" + ("within_" + LAB[L].replace(" ", "") if L else "all_nuclei"), gene_id=genes, n_in=int(ink.sum()),
                          log2fc=lf, pct_in=pi, pct_out=po, score=score))
    d = d[np.isfinite(d.score) & (d.pct_in >= 0.1) & (d.log2fc > 0.25)].sort_values("score", ascending=False).head(200)
    d["rank"] = np.arange(1, len(d) + 1); res.append(d)
pd.concat(res).to_csv(O / "marker_sets_top200.tsv", sep="\t", index=False)
P("done")
