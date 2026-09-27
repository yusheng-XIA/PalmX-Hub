"""SV-GWAS for the eight SV-significant traits under the published SV model (SV kinship, intercept, SV PC1-5) with the
raw phenotype (reproduction; compared row by row with the published .ps) and with a rank-based inverse-normal
transformed phenotype (Blom). Bonferroni 0.05/370,136; 250-kb anchor intervals as published."""
import numpy as np, pandas as pd
from scipy import stats
from scipy.stats import chi2
from config import *
from emmaxpy import GLS, read_matrix, read_pheno
ids = [l.split()[1] for l in open(FAM)]; N = len(ids)
TR = ["C12_0_Lauric_acid", "C14_0_Myristic_acid", "C18_3n3_Alpha_linolenic_acid", "Flesh_thickness_mm", "Nut_length_mm",
      "Shell_thickness_mm", "Shell_weight_g", "Stem_height_cm"]
tasks = pd.read_csv(RUN / "manifests/sv_tasks.tsv", sep="\t"); cat_of = dict(zip(tasks.trait, tasks.category))
svid, chrom, pos, G = [], [], [], []
code = {"1": 0, "2": 1, "0": -9, "N": -9}
with open(SV_TPED) as fh:
    for line in fh:
        f = line.split()
        a = np.array([code.get(x, -9) for x in f[4:]], dtype=np.int16).reshape(-1, 2)
        d = a.sum(1).astype(np.float32); d[(a < 0).any(1)] = np.nan
        svid.append(f[1]); chrom.append(f[0]); pos.append(int(f[3])); G.append(d)
G = np.vstack(G); chrom = np.array(chrom); pos = np.array(pos); svid = np.array(svid)
print("SV rows", len(svid), "unique", len(set(svid)), flush=True)
K = read_matrix(KIN_SV); X = pd.read_csv(W / "in/X_sv5pc.tsv", sep="\t", header=None, index_col=0).loc[ids].to_numpy()
def intervals(c, p):
    k = 0
    for cc in np.unique(c):
        anchor = None
        for x in np.sort(p[c == cc]):
            if anchor is None or x - anchor > 250_000: k += 1; anchor = x
    return k
rows = []
for tr in TR:
    y = read_pheno(RUN / f"phenotypes/{cat_of[tr]}/{tr}.txt", ids); keep = np.isfinite(y)
    z = np.full(N, np.nan); z[keep] = stats.norm.ppf((stats.rankdata(y[keep]) - 0.375) / (keep.sum() + 0.25))
    for lab, yy in [("raw", y), ("INT", z)]:
        g = GLS(yy[keep], X[keep], K[np.ix_(keep, keep)])
        p = np.concatenate([g.scan(G[s:s + 20000, keep])[2] for s in range(0, len(G), 20000)])
        _, first = np.unique(svid, return_index=True); u = np.zeros(len(p), bool); u[first] = True
        sig = (p < BONF_SV) & u
        lam = chi2.isf(np.median(p[u]), 1) / chi2.ppf(.5, 1)
        r = dict(trait=tr, phenotype=lab, n=int(keep.sum()), lambda_gc=lam, n_sig_sv=int(sig.sum()), n_intervals=intervals(chrom[sig], pos[sig]),
                 chroms=",".join(sorted(set(chrom[sig]))), min_p=float(p[u].min()))
        if lab == "raw":
            pub = pd.read_csv(RUN / f"sv/results/{cat_of[tr]}/{tr}/{tr}.ps", sep="\t", header=None, usecols=[0, 3])
            assert (pub[0].to_numpy() == svid).all()
            r["max_abs_dlog10P_vs_published"] = float(np.abs(np.log10(p) - np.log10(pub[3].to_numpy())).max())
        rows.append(r); print(r, flush=True)
pd.DataFrame(rows).to_csv(W / "out/sv_int_summary.tsv", sep="\t", index=False)
