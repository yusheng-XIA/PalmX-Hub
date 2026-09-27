"""Compare emmaxpy M1 (kinship + SNP PC1-5) with the EMMAX binary run with the same covariates (lauric acid, chr07B)."""
import numpy as np, pandas as pd, json
from config import *
mine = np.load(W / "out/scan/M1/C12_0_Lauric_acid/chr07B.npy").astype(float)
b = pd.read_csv(W / "tmp/bin/lauric_pc_chr07B.ps", sep="\t", header=None, usecols=[3])[3].to_numpy()
d = np.abs(mine + np.log10(np.clip(b, 1e-300, 1)))
reml = [float(x) for x in open(W / "tmp/bin/lauric_pc_chr07B.reml").read().split()]
print("markers", len(b), "max|dlog10P|", d.max(), "median", np.median(d), "delta binary", reml[2], "mine", json.load(open(W / "out/scan/M1/C12_0_Lauric_acid/reml.json"))["delta"],
      "| n P<Bonf binary", int((b < BONF_SNP).sum()), "mine", int((mine > -np.log10(BONF_SNP)).sum()), "| min P binary", b.min(), "mine", 10 ** -mine.max())
