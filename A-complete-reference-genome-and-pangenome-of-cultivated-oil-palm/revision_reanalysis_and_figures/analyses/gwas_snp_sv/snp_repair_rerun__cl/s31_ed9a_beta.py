"""beta/se/P for every SNP displayed in ED9a (nut weight, model M0), so the Source Data sheet keeps its P/beta/se columns."""
import sys, numpy as np, pandas as pd
sys.path.insert(0, "${CLUSTER_WORK}/snp_repair_rerun/gwas/cl")
from config import *
from emmaxpy import GLS, read_bed_rows, read_matrix, read_pheno
O = W / "ed9"; tr = "Nut_weight_g"; D = W / "out/scan/M0" / tr
rng = pd.read_csv(W / "in/chr_ranges.tsv", sep="\t").set_index("chr")
ids = [l.split()[1] for l in open(FAM)]; N = len(ids)
y = read_pheno(RUN / f"phenotypes/yield/{tr}.txt", ids); k = np.isfinite(y)
K = read_matrix(KIN_SNP); g = GLS(y[k], np.ones((k.sum(), 1)), K[np.ix_(k, k)])
out = []
for c in CHROMS:
    nlp = np.load(D / f"{c}.npy").astype(float); pos = np.load(W / f"in/pos_{c}.npy")
    idx = np.nonzero((np.arange(len(nlp)) % 200 == 0) | (nlp > 5))[0]
    st = int(rng.loc[c, "start"]); B = 200000
    for b0 in range(0, len(nlp), B):
        sel = idx[(idx >= b0) & (idx < b0 + B)]
        if not len(sel): continue
        G = read_bed_rows(BED, N, st + b0, st + min(b0 + B, len(nlp)))[sel - b0].astype(float); G[G < 0] = np.nan
        be, se, p = g.scan(G[:, k])
        out.append(pd.DataFrame({"variant_type": "SNP", "chrom": c, "pos": pos[sel], "P": p, "beta": be, "se": se, "nlp_scan": nlp[sel]}))
    print(c, flush=True)
d = pd.concat(out)
dev = np.abs(-np.log10(d.P) - d.nlp_scan).max(); print("rows", len(d), "max |dlog10P| vs scan", dev)
assert len(d) == 131072 and dev < 0.01
d.drop(columns="nlp_scan").to_csv(O / "ED9a_SNP_display_beta.tsv", sep="\t", index=False, float_format="%.6g")
print("done")
