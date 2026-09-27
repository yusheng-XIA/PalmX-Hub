#!/usr/bin/env python3
"""W sweep (P = 15, coverage mask T = 0.5, final Nigerian calls) with favourable loci re-derived from the reported
SV-GWAS; identical DP/masking to fix/fig5hi_mask (core_mask.py, run_mask.py)."""
import sys, os
H = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(H), "fig5hi_mask"))
import numpy as np, pandas as pd
import core_mask as K
from core_mask import chroms, donors, di
E5, _ = K.eligibility(0.5)
SD, _ = K.single_donor_loads(K.L1, E5)
BSL = min(v["Imputed_Total"] for v in SD.values()); best = min(donors, key=lambda d: (SD[d]["Imputed_Total"], d))
C5 = K.masked_cost(K.L1, E5)
old_F, old_n = K.F0, len(K.loci)

new_loci = pd.read_csv(f"{H}/data/fav_loci_reported.tsv", sep="\t")
al = pd.read_csv(f"{H}/data/donor_alleles.tsv", sep="\t")
K.alle = {}
for r in al.itertuples():
    K.alle.setdefault(r.sv_id, {})[r.donor] = r.carries_ALT == 1
K.loci = new_loci[["SV", "Chrom", "Pos", "Target"]]
new_F = K.fav_matrix(False)

rows, paths = [], {}
for lab, F, n in (("old_284", old_F, old_n), ("reported_gwas", new_F, len(new_loci))):
    for W in (0, 1, 2, 4, 8):
        p, o = K.run_weighted(W, K.P0, F, C5)
        ld, dsv, bp, mx, used = K.stats(p, K.L1)
        c1 = int(F["chr01B"][np.arange(len(p["chr01B"])), p["chr01B"]].sum())
        rows.append(dict(Set=lab, W=W, Load=ld, dSV=dsv, dSNP=ld - dsv, Breakpoints=bp, Segments=bp + len(chroms),
                         Donors=len(used), Reduction_pct=round(100 * (BSL - ld) / BSL, 2),
                         Captured=K.capture(p, F), Total=n, Captured_chr01B=c1,
                         MaxAttainable=int(sum((F[c] * E5[c]).max(1).sum() for c in chroms)),
                         Ineligible=K.n_ineligible(p, E5)))
        paths[(lab, W)] = p
S = pd.DataFrame(rows)
for lab, g in S.groupby("Set"):
    b = g[g.W == 0].iloc[0]
    S.loc[g.index, "dLoad_vs_W0"] = g.Load - b.Load; S.loc[g.index, "dCapture_vs_W0"] = g.Captured - b.Captured
print("best single:", best, BSL)
print(S.to_string())
S.to_csv(f"{H}/weight_sweep_compare.tsv", sep="\t", index=False)
# window-level path difference W=4
a, b = paths[("old_284", 4)], paths[("reported_gwas", 4)]
diff = sum(int((a[c] != b[c]).sum()) for c in chroms)
print("windows with a different donor (W=4, old vs new):", diff, "of", sum(len(a[c]) for c in chroms),
      {c: int((a[c] != b[c]).sum()) for c in chroms if (a[c] != b[c]).any()})
# old-path capture of new loci and new-path capture of old loci
print("old W4 path captures new loci:", K.capture(a, new_F), "/", len(new_loci),
      "| new W4 path captures old loci:", K.capture(b, old_F), "/", old_n)
# loci without donor call on new W4 path
nc = 0
for r in K.loci.itertuples():
    w = r.Pos // K.WIN; d = donors[b[r.Chrom][w]]
    if K.alle.get(r.SV, {}).get(d) is None: nc += 1
print("new W4: loci whose selected donor has no call:", nc)
np.save(f"{H}/path_W4_reported.npy", np.concatenate([b[c] for c in chroms]))
