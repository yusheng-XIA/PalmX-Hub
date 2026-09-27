# One-vs-rest Wilcoxon markers for the 20 snRNA clusters, following Seurat v5 FindAllMarkers
# (LogNormalize, scale 1e4; avg_log2FC = log2((sum(expm1) + 1)/n) difference; Bonferroni p_val_adj).
import sys, numpy as np, pandas as pd, scipy.io, scipy.sparse as sp
from scipy.stats import norm
E = "${CLUSTER_WORK}/enh_sn_pop/sn/export"
O = "${CLUSTER_WORK}/go_enrich/out"
genes = [l.strip() for l in open(f"{E}/genes.txt")]
meta = pd.read_csv(f"{E}/meta.tsv", sep="\t")
for c in [c for c in meta.columns if c.startswith("RNA_snn") or c == "seurat_clusters"]:
    print(c, meta[c].nunique(), flush=True)
cl = meta["seurat_clusters"].values
print(pd.Series(cl).value_counts().sort_index().to_dict(), flush=True)
X = scipy.io.mmread(f"{E}/counts.mtx").tocsc().astype(np.float64)  # genes x cells
lib = np.asarray(X.sum(0)).ravel()
X = X.multiply(1e4 / lib).tocsc(); X.data = np.log1p(X.data)
X = X.T.tocsc()  # cells x genes
n = X.shape[0]; G = X.shape[1]
clusters = sorted(np.unique(cl))
masks = {k: (cl == k) for k in clusters}
nin = {k: masks[k].sum() for k in clusters}
rows = []
for g in range(G):
    s, e = X.indptr[g], X.indptr[g + 1]
    idx = X.indices[s:e]; val = X.data[s:e]
    nz = len(idx)
    if nz == 0: continue
    # ranks: zeros tied at average rank (n-nz+1)/2; nonzeros ranked above
    nzero = n - nz
    order = np.argsort(val, kind="mergesort"); sv = val[order]
    r = np.empty(nz);
    # average ranks for ties among nonzeros
    uniq, inv, cnt = np.unique(sv, return_inverse=True, return_counts=True)
    cum = np.cumsum(cnt); start = cum - cnt
    avg = nzero + (start + cum + 1) / 2.0
    r[order] = avg[inv]
    zr = (nzero + 1) / 2.0
    tie = (nzero**3 - nzero) + np.sum(cnt**3 - cnt)
    ex = np.expm1(val)
    for k in clusters:
        m = masks[k][idx]
        n1 = nin[k]; n2 = n - n1
        pct1 = m.sum() / n1; pct2 = (nz - m.sum()) / n2
        if max(pct1, pct2) < 0.01: continue
        # Seurat v5 mean.fxn: log2((rowSums(expm1(x)) + 1) / ncol(x))
        lfc = np.log2((ex[m].sum() + 1) / n1) - np.log2((ex[~m].sum() + 1) / n2)
        if lfc < 0.1: continue  # only.pos, logfc.threshold 0.1
        R1 = r[m].sum() + (n1 - m.sum()) * zr
        U = R1 - n1 * (n1 + 1) / 2.0
        mu = n1 * n2 / 2.0
        sd = np.sqrt(n1 * n2 / 12.0 * ((n + 1) - tie / (n * (n - 1))))
        z = (U - mu - 0.5 * np.sign(U - mu)) / sd if sd > 0 else 0
        p = 2 * norm.sf(abs(z))
        rows.append((k, genes[g], lfc, round(pct1, 3), round(pct2, 3), p))
    if g % 2000 == 0: print("gene", g, flush=True)
df = pd.DataFrame(rows, columns=["cluster", "gene", "avg_log2FC", "pct.1", "pct.2", "p_val"])
df["p_val_adj"] = np.minimum(df["p_val"] * G, 1.0)
df.sort_values(["cluster", "p_val", "avg_log2FC"], ascending=[True, True, False]).to_csv(f"{O}/sn_markers_all.tsv", sep="\t", index=False)
# detection per gene (for background)
det = np.diff(X.indptr) / n
pd.DataFrame({"gene": genes, "pct_all": det}).to_csv(f"{O}/sn_gene_detection.tsv", sep="\t", index=False)
print("done", len(df), flush=True)
