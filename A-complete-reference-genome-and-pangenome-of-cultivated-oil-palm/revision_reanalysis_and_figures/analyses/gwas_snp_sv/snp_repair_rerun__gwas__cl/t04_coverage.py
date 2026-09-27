"""SNP-GWAS marker coverage per chromosome: max position, and gaps > 2 Mb between consecutive markers; SV marker span."""
import numpy as np, pandas as pd
from config import *
sv = pd.read_csv(SVK / "sv_assoc.map", sep=r"\s+", header=None, names=["chr", "id", "cm", "pos"])
for c in CHROMS:
    p = np.load(W / f"in/pos_{c}.npy"); d = np.diff(p); big = np.nonzero(d > 2_000_000)[0]
    s = sv[sv.chr == c].pos
    print(c, "SNPs", len(p), "span", p.min(), p.max(), "| SVs", len(s), "span", s.min(), s.max(),
          "| gaps>2Mb:", [(int(p[i]), int(p[i + 1])) for i in big][:6])
