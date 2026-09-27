"""Chromosome row ranges of the PLINK bed (checked against the EMMAX marker order), covariates and K = 4 groups."""
import numpy as np, pandas as pd
from config import *
(W / "in").mkdir(exist_ok=True)
bim = pd.read_csv(BIM, sep="\t", header=None, usecols=[0, 1, 3], names=["chr", "id", "pos"], dtype={"chr": str, "id": str, "pos": np.int64})
print("bim rows", len(bim))
rows = []
for c, g in bim.groupby("chr", sort=False):
    s, e = g.index.min(), g.index.max() + 1
    assert e - s == len(g), f"{c} not contiguous"
    np.save(W / f"in/pos_{c}.npy", g.pos.to_numpy())
    man = pd.read_csv(RUN / f"manifests/snp_markers_{c}.tsv", sep="\t", header=None, names=["id", "pos"])
    same = len(man) == len(g) and (man.id.to_numpy() == g.id.to_numpy()).all()
    rows.append((c, s, e, len(g), same))
    print(c, s, e, len(g), "order==EMMAX manifest:", same)
pd.DataFrame(rows, columns=["chr", "start", "stop", "n", "order_ok"]).to_csv(W / "in/chr_ranges.tsv", sep="\t", index=False)
ids = [l.split()[1] for l in open(FAM)]
assert ids == [l.split()[1] for l in open(SV_TFAM)], "SNP/SV sample order differ"
pca = pd.read_csv(SNP_PCA, sep=r"\s+", header=None).set_index(1)
pd.DataFrame(np.c_[np.ones(len(ids)), pca.loc[ids, [2, 3, 4, 5, 6]].to_numpy()], index=ids).to_csv(W / "in/X_snp5pc.tsv", sep="\t", header=False)
sv = pd.read_csv(SV_COV, sep=r"\s+", header=None).set_index(1)
assert list(sv.index) == ids
pd.DataFrame(sv.iloc[:, 1:].to_numpy(), index=ids).to_csv(W / "in/X_sv5pc.tsv", sep="\t", header=False)
qf = [l.split()[1] for l in open(QFAM)]
q = pd.read_csv(Q4, sep=r"\s+", header=None); q.index = qf
g = pd.DataFrame({"pop": q.to_numpy().argmax(1) + 1, "maxQ": q.to_numpy().max(1)}, index=qf).loc[ids]
g.to_csv(W / "in/k4_groups.tsv", sep="\t")
print("pop sizes", g["pop"].value_counts().sort_index().to_dict(), "maxQ<0.7", int((g.maxQ < 0.7).sum()))
r = np.corrcoef(np.c_[pca.loc[ids, [2, 3, 4, 5, 6]].to_numpy(), sv.iloc[:, 2:].to_numpy()].T)[:5, 5:]
print("|corr| SNP PCs (rows) vs SV PCs (cols)\n", np.round(np.abs(r), 3))
