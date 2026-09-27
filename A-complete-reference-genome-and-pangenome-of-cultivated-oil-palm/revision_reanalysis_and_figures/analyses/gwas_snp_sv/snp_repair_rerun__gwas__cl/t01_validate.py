"""Check emmaxpy against published EMMAX output: REML delta and P for the first 50,000 markers of a chromosome."""
import sys, numpy as np, pandas as pd
from config import *
from emmaxpy import GLS, read_bed_rows, read_matrix, read_pheno
ids = [l.split()[1] for l in open(FAM)]; N = len(ids); K = read_matrix(KIN_SNP)
rng = pd.read_csv(W / "in/chr_ranges.tsv", sep="\t").set_index("chr")
tasks = pd.read_csv(RUN / "manifests/sv_tasks.tsv", sep="\t").set_index("trait")
for tr, c in [("Nut_weight_g", "chr01B"), ("C12_0_Lauric_acid", "chr07B"), ("Flesh_thickness_mm", "chr16B")]:
    cat = tasks.loc[tr, "category"]
    y = read_pheno(RUN / f"phenotypes/{cat}/{tr}.txt", ids); keep = np.isfinite(y)
    g = GLS(y[keep], np.ones((keep.sum(), 1)), K[np.ix_(keep, keep)])
    reml = [float(x) for x in open(RUN / f"snp/results/{cat}/{tr}/emmax_{c}.reml").read().split()]
    s = int(rng.loc[c, "start"]); M = 50000
    if tr == "Nut_weight_g":   # include the lead SNP chr01B:3,153,030
        pos = np.load(W / f"in/pos_{c}.npy"); off = int(np.searchsorted(pos, 3153030)) - 25000
    else:
        off = 0
    G = read_bed_rows(BED, N, s + off, s + off + M)
    b, se, p = g.scan(G[:, keep])
    ps = pd.read_csv(RUN / f"snp/results/{cat}/{tr}/emmax_{c}.ps", sep="\t", header=None, skiprows=off, nrows=M)
    d = np.abs(np.log10(p) - np.log10(ps[3].to_numpy()))
    print(tr, c, "n", keep.sum(), "delta mine/EMMAX", round(g.r["delta"], 6), reml[2], "REML LL", round(g.r["ll"], 5), reml[0],
          "| max|dlog10P|", d.max(), "median", np.median(d), "| min P mine/EMMAX", p.min(), ps[3].min(),
          "| |beta| corr", np.corrcoef(np.abs(b), np.abs(ps[1]))[0, 1])
