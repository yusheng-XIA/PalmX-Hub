"""In-sample minor-allele count of every Bonferroni-significant SV (published SV-GWAS) among the phenotyped accessions."""
import numpy as np, pandas as pd
from config import *
from emmaxpy import read_pheno, tped_dosage
ids = [l.split()[1] for l in open(FAM)]
tasks = pd.read_csv(RUN / "manifests/sv_tasks.tsv", sep="\t")
sig = pd.concat([pd.read_csv(RUN / f"sv/results/{c}/{t}/significant_svs.tsv", sep="\t").assign(trait=t, category=c)
                 for c, t in zip(tasks.category, tasks.trait)
                 if (RUN / f"sv/results/{c}/{t}/significant_svs.tsv").exists()], ignore_index=True)
sig = sig.drop_duplicates(["trait", "SV"])
need = set(sig.SV); g = {}
with open(SV_TPED) as fh:
    for line in fh:
        f = line.split(maxsplit=4)
        if f[1] in need and f[1] not in g: g[f[1]] = tped_dosage(f[4].split())
rows = []
for r in sig.itertuples(index=False):
    k = np.isfinite(read_pheno(RUN / f"phenotypes/{r.category}/{r.trait}.txt", ids)); d = g[r.SV][k]; ok = ~np.isnan(d)
    ac = d[ok].sum(); mac = min(ac, 2 * ok.sum() - ac)
    rows.append(dict(trait=r.trait, SV=r.SV, P=r.P, n=int(k.sum()), n_called=int(ok.sum()), insample_mac=int(mac),
                     n_minor_hom_or_het=int(((d[ok] != (2 if ac > ok.sum() else 0))).sum())))
O = pd.DataFrame(rows); O.to_csv(W / "out/sv_sig_insample_mac.tsv", sep="\t", index=False)
print(O.groupby("trait").insample_mac.describe().to_string())
print(O.groupby("trait").apply(lambda x: (x.insample_mac <= 5).sum()).to_string())
