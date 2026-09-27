"""EMMAX scans of one chromosome for all trait x model combinations.
M0  kinship only (published model; reproduction)            all 60 traits
M1  kinship + intercept + SNP PC1-5 (PLINK PCA, Fig. 4b)   all 60 traits
M1v kinship + intercept + SV PC1-5 (as in the SV model)    focal traits
M2  M1 without K4-Pop4 accessions                           focal traits
M3  M1 without accessions with max K = 4 ancestry < 0.7    focal traits
M1i M1 on rank-based inverse-normal-transformed phenotype  focal traits
Output per model/trait/chrom: -log10 P (float32 .npy, EMMAX marker order) and rows with P < 1e-5."""
import sys, json, time, numpy as np, pandas as pd
from scipy import stats
from config import *
from emmaxpy import GLS, read_bed_rows, read_matrix, read_pheno
chrom = sys.argv[1]
only_models = sys.argv[2].split(",") if len(sys.argv) > 2 else None
ids = [l.split()[1] for l in open(FAM)]; N = len(ids)
K = read_matrix(KIN_SNP); assert K.shape == (N, N)
Xp = pd.read_csv(W / "in/X_snp5pc.tsv", sep="\t", header=None, index_col=0).loc[ids].to_numpy()
Xv = pd.read_csv(W / "in/X_sv5pc.tsv", sep="\t", header=None, index_col=0).loc[ids].to_numpy()
grp = pd.read_csv(W / "in/k4_groups.tsv", sep="\t", index_col=0).loc[ids]
tasks = pd.read_csv(RUN / "manifests/sv_tasks.tsv", sep="\t")
models = []
for cat, tr in zip(tasks.category, tasks.trait):
    y = read_pheno(RUN / f"phenotypes/{cat}/{tr}.txt", ids)
    models.append(("M0", tr, y, np.ones((N, 1))))
    models.append(("M1", tr, y, Xp))
    if tr in FOCAL:
        models.append(("M1v", tr, y, Xv))
        models.append(("M2", tr, np.where(grp["pop"].to_numpy() == 4, np.nan, y), Xp))
        models.append(("M3", tr, np.where(grp["maxQ"].to_numpy() < 0.7, np.nan, y), Xp))
        ok = np.isfinite(y); z = np.full(N, np.nan)          # rank-based inverse normal transform (Blom)
        z[ok] = stats.norm.ppf((stats.rankdata(y[ok]) - 0.375) / (ok.sum() + 0.25))
        models.append(("M1i", tr, z, Xp))
if only_models:
    models = [m for m in models if m[0] in only_models]
if len(sys.argv) > 3:                                     # part k/n of the model list (parallel jobs per chromosome)
    k_, n_ = map(int, sys.argv[3].split("/")); models = models[k_::n_]
gl = []
for mod, tr, y, X in models:
    keep = np.isfinite(y)
    g = GLS(y[keep], X[keep], K[np.ix_(keep, keep)])
    od = W / "out/scan" / mod / tr; od.mkdir(parents=True, exist_ok=True)
    if chrom == "chr01B":   # (each part writes its own models)
        json.dump(dict(n=int(keep.sum()), q=int(X.shape[1]), **g.r), open(od / "reml.json", "w"))
    gl.append((mod, tr, keep, g, od))
rng = pd.read_csv(W / "in/chr_ranges.tsv", sep="\t").set_index("chr").loc[chrom]
pos = np.load(W / f"in/pos_{chrom}.npy")
nlp = {k: np.empty(rng.n, dtype=np.float32) for k in range(len(gl))}
hits = {k: [] for k in range(len(gl))}
CH = 50000; t0 = time.time()
for s in range(0, int(rng.n), CH):
    e = min(int(rng.n), s + CH)
    G = read_bed_rows(BED, N, int(rng.start) + s, int(rng.start) + e)
    for k, (mod, tr, keep, g, od) in enumerate(gl):
        b, se, p = g.scan(G[:, keep])
        nlp[k][s:e] = -np.log10(np.clip(p, 1e-300, 1))
        h = np.nonzero(p < 1e-5)[0]
        if len(h):
            hits[k].append(pd.DataFrame({"idx": s + h, "pos": pos[s + h], "beta": b[h], "se": se[h], "p": p[h]}))
    print(chrom, e, "/", rng.n, f"{time.time() - t0:.0f}s", flush=True)
for k, (mod, tr, keep, g, od) in enumerate(gl):
    np.save(od / f"{chrom}.npy", nlp[k])
    (pd.concat(hits[k]) if hits[k] else pd.DataFrame(columns=["idx", "pos", "beta", "se", "p"])).to_csv(od / f"{chrom}.hits.tsv", sep="\t", index=False)
print("done", chrom, f"{time.time() - t0:.0f}s")
