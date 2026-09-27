"""r2 between the two ED9c focal SV records and the lead SNP chr01B:3,153,030 (repaired call set), nut-weight cohort."""
import numpy as np, pandas as pd
from config import *
from emmaxpy import read_bed_rows, read_pheno, tped_dosage
ids = [l.split()[1] for l in open(FAM)]; N = len(ids)
y = read_pheno(RUN / "phenotypes/yield/Nut_weight_g.txt", ids); k = np.isfinite(y)
rng = pd.read_csv(W / "in/chr_ranges.tsv", sep="\t").set_index("chr"); pos = np.load(W / "in/pos_chr01B.npy")
i = int(np.nonzero(pos == 3153030)[0][0]); L = read_bed_rows(BED, N, int(rng.loc["chr01B", "start"]) + i, int(rng.loc["chr01B", "start"]) + i + 1)[0].astype(float)
L[L < 0] = np.nan
ids_sv = [l.split()[1] for l in open(SV_TFAM)]; assert ids_sv == ids
want = {}
with open(SV_TPED) as fh:
    for line in fh:
        f = line.split(maxsplit=4)
        if f[0] == "chr01B" and 3050000 <= int(f[3]) <= 3410000: want[f[1]] = (int(f[3]), tped_dosage(f[4].split()))
for sid, (p, g) in want.items():
    ok = k & ~np.isnan(g) & ~np.isnan(L)
    print(sid, p, ("%.6f" % (np.corrcoef(g[ok], L[ok])[0, 1] ** 2)) if ok.sum() > 3 and np.std(g[ok]) > 0 else "nan", int(ok.sum()), sep="\t")
