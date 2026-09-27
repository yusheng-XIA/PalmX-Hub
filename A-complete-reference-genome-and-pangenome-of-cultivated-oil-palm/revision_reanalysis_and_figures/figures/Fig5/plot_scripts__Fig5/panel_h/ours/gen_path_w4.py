#!/usr/bin/env python3
"""Generate the African35 donor path at W = 4, P = 15 with the enh_C objective and write the
window-, segment-, chromosome- and locus-level tables used by Figure 5h,i, Source Data and the text.

Objective (enh_C.py / 71_build_loadonly_iph.py::exact_dp, per chromosome, 500-kb windows):
    sum_w [ L(w, d_w) - W * F(w, d_w) ] + P * (#switches),   W = 4, P = 15
F(w, d) = favourable loci in window w for which donor d carries the target allele (step9 donor calls;
donor without a call -> 0). Capture on the resulting path uses the same definition (compute_capture.py).
The functions are executed verbatim from ../enh_C/enh_C.py (sections 0-1), so tie-breaking is identical.
"""
import csv, os, sys
from collections import Counter, defaultdict
import numpy as np, pandas as pd

H = os.path.dirname(os.path.abspath(__file__))
ENH = os.path.join(os.path.dirname(H), "enh_C", "enh_C.py")
W4, P = 4, 15
src = open(ENH).read()
ns = {"__file__": ENH, "__name__": "enh_C_exec"}
exec(src[:src.index("# ---------- 2a.")], ns)          # data, F0/Fp, exact_dp, run_weighted, stats, capture
L, DSV, F0, Fp, chroms, donors, di = (ns[k] for k in ("L", "DSV", "F0", "Fp", "chroms", "donors", "di"))
BSL, best_single, loci, alle, PROXY = ns["BSL"], ns["best_single"], ns["loci"], ns["alle"], ns["PROXY"]
m = ns["m"]
assert ns["repro"]["Identical_windows"] == 3486 and ns["repro"]["Load"] == 23520   # W = 0 path reproduced

p4, obj4 = ns["run_weighted"](W4, P, F0)
p0, obj0 = ns["run_weighted"](0, P, F0)
load, dsv, bp, mx, used = ns["stats"](p4)
cap4, capP4 = ns["capture"](p4, F0), ns["capture"](p4, Fp)
print(f"W=4: load {load} dsv {dsv} dsnp {load-dsv} bp {bp} seg {bp+len(chroms)} donors {len(used)} obj {obj4} "
      f"capture {cap4} proxy {capP4} red {100*(BSL-load)/BSL:.4f}% (best single {best_single} {BSL})")

# ---- window-level path (same columns as results-8.9 ideal_loadonly_path.tsv + reward columns)
DSNP = {}
for c in chroms:
    s = m[m.Chrom == c]
    DSNP[c] = s.pivot(index="Window_Index", columns="Sample_ID", values="DSNP_Count")[donors].to_numpy(np.int64)
ends = {c: m[m.Chrom == c].groupby("Window_Index").Window_End_0based.first().to_numpy() for c in chroms}
rows = []; gcum = 0; bych = []; segs = []
for c in chroms:
    path = p4[c]; ccum = 0; n = len(path)
    for w in range(n):
        d = path[w]; sw = w > 0 and path[w] != path[w - 1]
        f = int(F0[c][w, d]); l = int(L[c][w, d])
        step = l - W4 * f + (P if sw else 0); ccum += step; gcum += step
        rows.append(dict(Panel="African35", Chrom=c, Window_Index=w, Window_Start_0based=w * 500_000,
                         Window_End_0based=int(ends[c][w]), Donor_ID=donors[d], DSV_Count=int(DSV[c][w, d]),
                         DSNP_Count=int(DSNP[c][w, d]), Total_Load=l, Fav_Loci_Carried=f,
                         Donor_Switch_From_Previous="Yes" if sw else "No", Cumulative_Chrom_Objective=ccum,
                         Cumulative_Genome_Objective=gcum, Breakpoint_Penalty_P=P, GWAS_Weight=W4))
    b = int((np.diff(path) != 0).sum())
    ld = int(L[c][np.arange(n), path].sum()); dv = int(DSV[c][np.arange(n), path].sum())
    bych.append(dict(Panel="African35", Chrom=c, Window_Count=n, Residual_DSV=dv, Residual_DSNP=ld - dv,
                     Residual_Total_Load=ld, Breakpoint_Count=b, Segment_Count=b + 1,
                     Donor_Count=len(set(path.tolist())), Objective=ccum, Breakpoint_Penalty_P=P, GWAS_Weight=W4))
    st = 0; k = 0
    for w in range(1, n + 1):
        if w == n or path[w] != path[st]:
            k += 1; d = path[st]
            segs.append(dict(Panel="African35", Segment_ID=f"{c}_SEG{k:04d}", Chrom=c, Start_Window_Index=st,
                             End_Window_Index_Exclusive=w, Segment_Start_0based=st * 500_000,
                             Segment_End_0based=int(ends[c][w - 1]), Donor_ID=donors[d], Window_Count=w - st,
                             Segment_Length_bp=int(ends[c][w - 1]) - st * 500_000,
                             DSV_Count=int(DSV[c][st:w, d].sum()), DSNP_Count=int(DSNP[c][st:w, d].sum()),
                             Total_Load=int(L[c][st:w, d].sum()), Breakpoint_Penalty_P=P, GWAS_Weight=W4))
            st = w
PW = pd.DataFrame(rows); BC = pd.DataFrame(bych); SG = pd.DataFrame(segs)
assert gcum == obj4 and BC.Objective.sum() == obj4 and len(SG) == bp + len(chroms)
PW.to_csv(f"{H}/path_W4.tsv", sep="\t", index=False)
SG.to_csv(f"{H}/segments_W4.tsv", sep="\t", index=False)
tot = dict(Panel="African35", Chrom="TOTAL", Window_Count=int(BC.Window_Count.sum()), Residual_DSV=int(BC.Residual_DSV.sum()),
           Residual_DSNP=int(BC.Residual_DSNP.sum()), Residual_Total_Load=int(BC.Residual_Total_Load.sum()),
           Breakpoint_Count=int(BC.Breakpoint_Count.sum()), Segment_Count=int(BC.Segment_Count.sum()),
           Donor_Count=len(used), Objective=int(BC.Objective.sum()), Breakpoint_Penalty_P=P, GWAS_Weight=W4)
pd.concat([BC, pd.DataFrame([tot])]).to_csv(f"{H}/by_chrom_W4.tsv", sep="\t", index=False)

# ---- donor contribution (for the 5h shading order: windows per donor, descending)
cnt = PW.Donor_ID.value_counts()
pd.DataFrame(dict(Donor_ID=cnt.index, Window_Count=cnt.values)).to_csv(f"{H}/donor_windows_W4.tsv", sep="\t", index=False)

# ---- locus-level capture (compute_capture.py definition)
WIN = 500_000
pathd = {c: {w: donors[d] for w, d in enumerate(p4[c])} for c in chroms}
path0 = {c: {w: donors[d] for w, d in enumerate(p0[c])} for c in chroms}
def state(sv, tgt, donor):
    ca = alle.get(sv, {}).get(donor)
    if ca is None: return None
    return int((tgt == "ALT" and ca) or (tgt == "REF" and not ca))
out = []
for r in loci.itertuples():
    w = r.Pos // WIN; don = pathd[r.Chrom][w]; s = state(r.SV, r.Target, don)
    s0 = state(r.SV, r.Target, path0[r.Chrom][w])
    pdn = PROXY.get(don); sp = state(r.SV, r.Target, pdn) if pdn else s
    out.append(dict(SV=r.SV, Chrom=r.Chrom, Pos=r.Pos, Window_Index=w, Target=r.Target, W4_Selected_Donor=don,
                    Selected_Donor_Carries_ALT="NA" if alle.get(r.SV, {}).get(don) is None else int(alle[r.SV][don]),
                    Captured_W4=1 if s == 1 else 0,
                    Capture_Status_W4="Captured" if s == 1 else ("Missed" if s == 0 else "Unknown_no_donor_call"),
                    Captured_W4_with_collapsed_proxy="NA" if sp is None else sp,
                    W0_Selected_Donor=path0[r.Chrom][w], Captured_W0=1 if s0 == 1 else 0))
CP = pd.DataFrame(out); CP.to_csv(f"{H}/capture_W4_path.tsv", sep="\t", index=False)
assert CP.Captured_W4.sum() == cap4 and CP.Captured_W0.sum() == 168
per = []
for c in chroms:
    g = CP[CP.Chrom == c]
    per.append(dict(Chrom=c, Fav_Total=len(g), Captured_W4=int(g.Captured_W4.sum()),
                    Missed_W4=int((g.Capture_Status_W4 == "Missed").sum()),
                    Unknown_W4_no_donor_call=int((g.Capture_Status_W4 == "Unknown_no_donor_call").sum()),
                    Capture_Rate_W4_pct=round(100 * g.Captured_W4.sum() / len(g), 2) if len(g) else None,
                    Captured_W4_collapsed_proxy=int((g.Captured_W4_with_collapsed_proxy == 1).sum()),
                    Captured_W0=int(g.Captured_W0.sum())))
PC = pd.DataFrame(per)
PC.loc[len(PC)] = dict(Chrom="TOTAL", Fav_Total=len(CP), Captured_W4=int(CP.Captured_W4.sum()),
                       Missed_W4=int((CP.Capture_Status_W4 == "Missed").sum()),
                       Unknown_W4_no_donor_call=int((CP.Capture_Status_W4 == "Unknown_no_donor_call").sum()),
                       Capture_Rate_W4_pct=round(100 * CP.Captured_W4.sum() / len(CP), 2),
                       Captured_W4_collapsed_proxy=int((CP.Captured_W4_with_collapsed_proxy == 1).sum()),
                       Captured_W0=int(CP.Captured_W0.sum()))
PC.to_csv(f"{H}/per_chrom_capture_W4.tsv", sep="\t", index=False)
mk = {}; mk0 = {}
for o in out:
    k = (o["Chrom"], o["Pos"]); mk[k] = max(mk.get(k, 0), o["Captured_W4"]); mk0[k] = max(mk0.get(k, 0), o["Captured_W0"])
with open(f"{H}/position_marks_W4.tsv", "w") as fh:
    fh.write("Chrom\tPos\tMark_W4\tMark_W0\tChanged\n")
    for k in sorted(mk, key=lambda k: (int(k[0][3:5]), k[1])):
        fh.write(f"{k[0]}\t{k[1]}\t{mk[k]}\t{mk0[k]}\t{int(mk[k] != mk0[k])}\n")

# ---- W0 vs W4 window differences
diff = [(c, w, donors[p0[c][w]], donors[p4[c][w]]) for c in chroms for w in range(len(p4[c])) if p0[c][w] != p4[c][w]]
with open(f"{H}/path_diff_W0_vs_W4.tsv", "w") as fh:
    fh.write("Chrom\tWindow_Index\tDonor_W0\tDonor_W4\n")
    for d in diff: fh.write("\t".join(map(str, d)) + "\n")
print("windows changed vs W0:", len(diff), Counter(d[0] for d in diff))
print(PC.to_string())
print(BC[["Chrom", "Residual_DSV", "Residual_DSNP", "Breakpoint_Count", "Segment_Count", "Donor_Count", "Objective"]].to_string())
print("positions", len(mk), "gold", sum(mk.values()), "changed vs W0", sum(mk[k] != mk0[k] for k in mk))
print("unknown W4 by donor", Counter(o["W4_Selected_Donor"] for o in out if o["Capture_Status_W4"].startswith("Unknown")))
