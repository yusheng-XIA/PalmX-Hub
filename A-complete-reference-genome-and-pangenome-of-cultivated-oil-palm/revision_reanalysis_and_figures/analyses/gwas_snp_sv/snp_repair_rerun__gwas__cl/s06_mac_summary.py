"""lambda_GC, Bonferroni hits and intervals after restricting to markers with in-sample MAF >= 0.05 (phenotyped
accessions of each trait), and the in-sample minor-allele count (MAC) of Bonferroni-significant SNPs; models M0, M1."""
import numpy as np, pandas as pd
from scipy.stats import chi2
from config import *
GAP = 250_000
def n_intervals(ch, pos):
    k = 0
    for c in np.unique(ch):
        p = np.sort(pos[ch == c]); anchor = None
        for x in p:
            if anchor is None or x - anchor > GAP: k += 1; anchor = x
    return k
rows = []
for tdir in sorted((W / "out/mac").iterdir()):
    tr = tdir.name
    mac = np.concatenate([np.load(tdir / f"{c}.mac.npy") for c in CHROMS]).astype(np.int32)
    n2 = np.concatenate([np.load(tdir / f"{c}.n2.npy") for c in CHROMS]).astype(np.int32)
    maf = np.where(n2 > 0, mac / np.maximum(n2, 1), 0)
    chrom = np.concatenate([[c] * len(np.load(tdir / f"{c}.mac.npy")) for c in CHROMS])
    pos = np.concatenate([np.load(W / f"in/pos_{c}.npy") for c in CHROMS])
    for mod in ["M0", "M1"]:
        nlp = np.concatenate([np.load(W / f"out/scan/{mod}/{tr}/{c}.npy") for c in CHROMS]).astype(np.float64)
        sig = nlp > -np.log10(BONF_SNP)
        common = maf >= 0.05
        lam_all = chi2.isf(10 ** -np.median(nlp), 1) / chi2.ppf(.5, 1)
        lam_c = chi2.isf(10 ** -np.median(nlp[common]), 1) / chi2.ppf(.5, 1)
        s_c = sig & common
        rows.append(dict(trait=tr, model=mod, n_markers=len(nlp), n_insample_maf05=int(common.sum()),
                         lambda_all=lam_all, lambda_maf05=lam_c, n_sig=int(sig.sum()),
                         sig_mac_median=float(np.median(mac[sig])) if sig.any() else np.nan,
                         sig_mac_le5=int((sig & (mac <= 5)).sum()), sig_maf_lt05=int((sig & ~common).sum()),
                         n_sig_maf05=int(s_c.sum()), n_int_maf05=n_intervals(chrom[s_c], pos[s_c]) if s_c.any() else 0,
                         n_chr_maf05=len(np.unique(chrom[s_c]))))
        print(rows[-1], flush=True)
pd.DataFrame(rows).to_csv(W / "out/mac_summary.tsv", sep="\t", index=False)
