"""Genomic inflation of the published SNP-GWAS (all SNP tests per trait; EMMAX .ps files of 02_run_emmax.sh) and REML variance components."""
import sys, numpy as np, pandas as pd
from multiprocessing import Pool
from scipy.stats import chi2
from config import *
tasks = pd.read_csv(RUN / "manifests/sv_tasks.tsv", sep="\t")
def one(t):
    cat, tr = t
    d = RUN / "snp/results" / cat / tr
    p = np.concatenate([pd.read_csv(d / f"emmax_{c}.ps", sep="\t", header=None, usecols=[3], dtype=np.float64, engine="c")[3].to_numpy() for c in CHROMS])
    p = np.clip(p, 1e-300, 1)
    lam = float(np.median(chi2.isf(p, 1)) / chi2.ppf(0.5, 1))
    reml = [float(x) for x in open(d / "emmax_chr01B.reml").read().split()]
    return dict(category=cat, trait=tr, n_tests=len(p), lambda_gc_all=lam, n_sig=int((p < BONF_SNP).sum()),
                min_p=float(p.min()), reml_delta=reml[2], reml_vg=reml[3], reml_ve=reml[4], reml_h2=reml[5])
if __name__ == "__main__":
    with Pool(int(sys.argv[1]) if len(sys.argv) > 1 else 8) as pool:
        res = pool.map(one, list(zip(tasks.category, tasks.trait)))
    pd.DataFrame(res).to_csv(W / "out/orig_snp_lambda.tsv", sep="\t", index=False)
    print(pd.DataFrame(res).sort_values("lambda_gc_all").to_string())
