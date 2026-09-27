#!/usr/bin/env python3
"""Clustering stability of C6/C9/C10 on the published Seurat SNN graph (Harmony dims 1-15, k = 20):
Louvain (igraph multilevel with resolution; best of 3 starts) at resolution 0.4-1.2 and 20 x 80% subsamples at 1.0.
For each original cluster: best-matching new cluster (max Jaccard among the new clusters holding its nuclei), share of
the original nuclei in it, and the FL fraction of that cluster's 185-d nuclei."""
from pathlib import Path
import numpy as np, pandas as pd, igraph as ig, random
from sklearn.metrics import adjusted_rand_score as ARI
import os
W = Path(os.environ.get("SN_WORK", "snRNA_work")); X = W / "export"; O = W / "out"
meta = pd.read_csv(X / "meta.tsv", sep="\t", index_col=0); cells = meta.index.tolist()
cl = meta.seurat_clusters.astype(int).values
is185 = (meta.timepoint == "185d").values; isFL = (meta.variety == "FL").values
T = pd.read_csv(X / "graph_RNA_snn.tsv.gz", sep="\t"); T = T[T.i < T.j]
G = ig.Graph(n=len(cells), edges=list(zip(T.i.values - 1, T.j.values - 1)), directed=False); G.es["weight"] = T.x.values
def run(res, sub=None, seed=10086):
    g = G if sub is None else G.induced_subgraph(sub)
    random.seed(seed); ig.set_random_number_generator(random)
    bq, bm = None, None
    for _ in range(3):
        m = g.community_multilevel(weights="weight", resolution=res).membership
        q = g.modularity(m, weights="weight", resolution=res)
        if bq is None or q > bq: bq, bm = q, m
    return np.array(bm)
def match(k, lab, ids):
    o = set(ids[cl[ids] == k]); s = pd.Series(lab, index=ids)
    best = None
    for c in s[list(o)].value_counts().index[:3]:
        n = set(s.index[s.values == c]); j = len(o & n) / len(o | n)
        if best is None or j > best[0]:
            nn = np.array(sorted(n)); fl = isFL[nn][is185[nn]].mean() if is185[nn].any() else np.nan
            best = (j, len(o & n) / len(o), len(n), fl)
    return best
rows = []; memb = {}
for res in (0.4, 0.6, 0.8, 1.0, 1.2):
    m = run(res); memb[f"louvain_res{res}"] = m
    r = dict(analysis="resolution", resolution=res, rep=0, n_clusters=len(set(m)), ARI_vs_published=ARI(cl, m))
    for k in (6, 9, 10):
        j, f, n, fl = match(k, m, np.arange(len(cells)))
        r.update({f"C{k}_jaccard": j, f"C{k}_share_in_best": f, f"C{k}_best_size": n, f"C{k}_best_FLfrac185": fl})
    rows.append(r); print(r, flush=True)
for col in ("RNA_snn_res.0.3", "RNA_snn_res.0.5", "RNA_snn_res.0.8"):
    m = meta[col].values
    r = dict(analysis="Seurat stored", resolution=float(col.split("res.")[1]), rep=0, n_clusters=len(set(m)), ARI_vs_published=ARI(cl, m))
    for k in (6, 9, 10):
        j, f, n, fl = match(k, m, np.arange(len(cells)))
        r.update({f"C{k}_jaccard": j, f"C{k}_share_in_best": f, f"C{k}_best_size": n, f"C{k}_best_FLfrac185": fl})
    rows.append(r); print(r, flush=True)
rng = np.random.default_rng(10086)
for b in range(20):
    sub = np.sort(rng.choice(len(cells), int(0.8 * len(cells)), replace=False))
    m = run(1.0, sub=sub, seed=10086 + b)
    r = dict(analysis="subsample80", resolution=1.0, rep=b + 1, n_clusters=len(set(m)), ARI_vs_published=ARI(cl[sub], m))
    for k in (6, 9, 10):
        j, f, n, fl = match(k, m, sub)
        r.update({f"C{k}_jaccard": j, f"C{k}_share_in_best": f, f"C{k}_best_size": n, f"C{k}_best_FLfrac185": fl})
    rows.append(r); print(r, flush=True)
pd.DataFrame(rows).to_csv(O / "cluster_stability.tsv", sep="\t", index=False)
pd.DataFrame(memb, index=cells).to_csv(O / "louvain_memberships.tsv.gz", sep="\t", compression="gzip")
print("done")
