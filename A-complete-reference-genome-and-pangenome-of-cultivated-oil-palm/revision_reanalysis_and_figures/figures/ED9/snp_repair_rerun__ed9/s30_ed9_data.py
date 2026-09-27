"""New-scan inputs for ED9a-c (nut weight, published SNP model M0, repaired call set)."""
import sys, numpy as np, pandas as pd
sys.path.insert(0, "${CLUSTER_WORK}/snp_repair_rerun/gwas/cl")
from config import *
from emmaxpy import GLS, read_bed_rows, read_matrix, read_pheno
O = W / "ed9"; O.mkdir(exist_ok=True)
tr = "Nut_weight_g"; D = W / "out/scan/M0" / tr
rng = pd.read_csv(W / "in/chr_ranges.tsv", sep="\t").set_index("chr")
rows = []; allp = []
for c in CHROMS:
    nlp = np.load(D / f"{c}.npy").astype(float); pos = np.load(W / f"in/pos_{c}.npy")
    allp.append(nlp)
    keep = (np.arange(len(nlp)) % 200 == 0) | (nlp > 5)
    rows.append(pd.DataFrame({"variant_type": "SNP", "chrom": c, "pos": pos[keep], "neglog10_p": nlp[keep]}))
man = pd.concat(rows); man.to_csv(O / "ED9a_SNP_display.tsv", sep="\t", index=False)
a = np.sort(np.concatenate(allp))[::-1]; n = len(a); r = np.arange(1, n + 1)
exp = -np.log10(r / (n + 1))
sel = (a > 5) | (r % 250 == 0)
pd.DataFrame({"variant_type": "SNP", "expected_neglog10P": exp[sel], "observed_neglog10P": a[sel],
              "rank_source": np.where(a[sel] > 5, "exact_rank_all_P_lt_1e-5", "every_250th_rank")}).to_csv(O / "ED9b_SNP_QQ.tsv", sep="\t", index=False)
ids = [l.split()[1] for l in open(FAM)]; N = len(ids)
y = read_pheno(RUN / f"phenotypes/yield/{tr}.txt", ids); k = np.isfinite(y)
K = read_matrix(KIN_SNP); g = GLS(y[k], np.ones((k.sum(), 1)), K[np.ix_(k, k)])
c = "chr01B"; pos = np.load(W / f"in/pos_{c}.npy"); s, e = np.searchsorted(pos, [3050000, 3410001])
G = read_bed_rows(BED, N, int(rng.loc[c, "start"]) + s, int(rng.loc[c, "start"]) + e).astype(float); G[G < 0] = np.nan
b, se, p = g.scan(G[:, k])
li = int(np.nonzero(pos[s:e] == 3153030)[0][0]); L = G[li, k]
def r2(x):
    ok = ~np.isnan(x) & ~np.isnan(L)
    return np.corrcoef(x[ok], L[ok])[0, 1] ** 2 if ok.sum() > 3 and np.nanstd(x[ok]) > 0 else np.nan
reg = pd.DataFrame({"marker": [f"{c}:{x}" for x in pos[s:e]], "variant_type": "SNP", "chrom": c, "pos": pos[s:e], "p": p,
                    "neglog10_p": -np.log10(p), "beta": b, "se": se, "r2_to_lead": [r2(G[i, k]) for i in range(e - s)]})
reg.to_csv(O / "ED9c_regional_SNPs.tsv", sep="\t", index=False)
thr = 0.05 / N_SNP_TESTS
print("n tests", N_SNP_TESTS, "threshold %.4g" % thr, "display SNPs", len(man), "regional SNPs", len(reg),
      "lead P %.3g" % p[li], "n sig nut weight", int((a > -np.log10(thr)).sum()))
