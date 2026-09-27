#!/usr/bin/env python3
"""Window-, segment-, chromosome- and locus-level tables of the coverage-masked W = 4, P = 15 path (Fig. 5h,i).

Same layout as ../fig5hi_w4/gen_path_w4.py, on the masked matrix of core_mask.py (T = 0.5; Nigerian dSNP from
the final assemblies). The current (unmasked) W = 4 path is recomputed from the original matrix for comparison.
Outputs (this directory): path_mask.tsv, segments_mask.tsv, by_chrom_mask.tsv, donor_windows_mask.tsv,
capture_mask_path.tsv, per_chrom_capture_mask.tsv, position_marks_mask.tsv, path_diff_current_vs_mask.tsv,
zoom_matrix_chr01B.tsv (5i heatmap input: load, coverage, eligibility).
"""
from collections import Counter
import numpy as np, pandas as pd
import core_mask as K
from core_mask import chroms, donors, di

H = K.H; W4, P, T = 4, 15, 0.5
E, ALL = K.eligibility(T)
C = K.masked_cost(K.L1, E)
p4, obj4 = K.run_weighted(W4, P, K.F0, C)          # masked path (Fig. 5h,i)
pc, objc = K.run_weighted(W4, P, K.F_OLD, C)       # current figure path (masked, old 284 loci)
p0, obj0 = K.run_weighted(0, P, K.F0, C)           # masked W = 0 (capture comparison)
assert K.stats(pc, K.L1)[0] == 24723 and K.stats(pc, K.L1)[2] == 374 and K.capture(pc, K.F_OLD) == 202
L = K.L1
load, dsv, bp, mx, used = K.stats(p4, L)
assert K.n_ineligible(p4, E) == 0
sdl, _ = K.single_donor_loads(L, E)
bsd = min(donors, key=lambda d: (sdl[d]["Imputed_Total"], d)); BSL = sdl[bsd]["Imputed_Total"]
cap4 = K.capture(p4)
print(f"fav233 W4: load {load} dsv {dsv} dsnp {load-dsv} bp {bp} seg {bp+16} donors {len(used)} obj {obj4} capture {cap4} "
      f"best single {bsd} {BSL} red {100*(BSL-load)/BSL:.4f}%")

rows = []; gcum = 0; bych = []; segs = []
for c in chroms:
    path = p4[c]; ccum = 0; n = len(path)
    for w in range(n):
        d = path[w]; sw = w > 0 and path[w] != path[w - 1]
        f = int(K.F0[c][w, d]); l = int(L[c][w, d])
        step = l - W4 * f + (P if sw else 0); ccum += step; gcum += step
        rows.append(dict(Panel="African35", Chrom=c, Window_Index=w, Window_Start_0based=w * K.WIN,
                         Window_End_0based=int(K.ends[c][w]), Donor_ID=donors[d], DSV_Count=int(K.DSV[c][w, d]),
                         DSNP_Count=int(K.DSNP1[c][w, d]), Total_Load=l, Fav_Loci_Carried=f,
                         Donor_Switch_From_Previous="Yes" if sw else "No", Cumulative_Chrom_Objective=ccum,
                         Cumulative_Genome_Objective=gcum, Breakpoint_Penalty_P=P, GWAS_Weight=W4,
                         Selected_Donor_Coverage=round(float(K.COV[c][w, d]), 4),
                         Eligible_Donors=int(E[c][w].sum()), No_Donor_Above_Threshold="Yes" if ALL[c][w] else "No"))
    b = int((np.diff(path) != 0).sum())
    ld = int(L[c][np.arange(n), path].sum()); dv = int(K.DSV[c][np.arange(n), path].sum())
    bych.append(dict(Panel="African35", Chrom=c, Window_Count=n, Residual_DSV=dv, Residual_DSNP=ld - dv,
                     Residual_Total_Load=ld, Breakpoint_Count=b, Segment_Count=b + 1,
                     Donor_Count=len(set(path.tolist())), Objective=ccum, Breakpoint_Penalty_P=P, GWAS_Weight=W4,
                     Windows_No_Donor_Above_Threshold=int(ALL[c].sum())))
    st = 0; k = 0
    for w in range(1, n + 1):
        if w == n or path[w] != path[st]:
            k += 1; d = path[st]
            segs.append(dict(Panel="African35", Segment_ID=f"{c}_SEG{k:04d}", Chrom=c, Start_Window_Index=st,
                             End_Window_Index_Exclusive=w, Segment_Start_0based=st * K.WIN,
                             Segment_End_0based=int(K.ends[c][w - 1]), Donor_ID=donors[d], Window_Count=w - st,
                             Segment_Length_bp=int(K.ends[c][w - 1]) - st * K.WIN,
                             DSV_Count=int(K.DSV[c][st:w, d].sum()), DSNP_Count=int(K.DSNP1[c][st:w, d].sum()),
                             Total_Load=int(L[c][st:w, d].sum()), Breakpoint_Penalty_P=P, GWAS_Weight=W4))
            st = w
PW = pd.DataFrame(rows); BC = pd.DataFrame(bych); SG = pd.DataFrame(segs)
assert gcum == obj4 and len(SG) == bp + 16
PW.to_csv(f"{H}/path_mask.tsv", sep="\t", index=False)
SG.to_csv(f"{H}/segments_mask.tsv", sep="\t", index=False)
tot = dict(Panel="African35", Chrom="TOTAL", Window_Count=int(BC.Window_Count.sum()), Residual_DSV=int(BC.Residual_DSV.sum()),
           Residual_DSNP=int(BC.Residual_DSNP.sum()), Residual_Total_Load=int(BC.Residual_Total_Load.sum()),
           Breakpoint_Count=int(BC.Breakpoint_Count.sum()), Segment_Count=int(BC.Segment_Count.sum()),
           Donor_Count=len(used), Objective=int(BC.Objective.sum()), Breakpoint_Penalty_P=P, GWAS_Weight=W4,
           Windows_No_Donor_Above_Threshold=int(BC.Windows_No_Donor_Above_Threshold.sum()))
pd.concat([BC, pd.DataFrame([tot])]).to_csv(f"{H}/by_chrom_mask.tsv", sep="\t", index=False)
cnt = PW.Donor_ID.value_counts()
pd.DataFrame(dict(Donor_ID=cnt.index, Window_Count=cnt.values)).to_csv(f"{H}/donor_windows_mask.tsv", sep="\t", index=False)

# ---- locus-level capture
pathd = {c: {w: donors[d] for w, d in enumerate(p4[c])} for c in chroms}
pathc = {c: {w: donors[d] for w, d in enumerate(pc[c])} for c in chroms}
def state(sv, tgt, donor):
    ca = K.alle.get(sv, {}).get(donor)
    if ca is None: return None
    return int((tgt == "ALT" and ca) or (tgt == "REF" and not ca))
out = []
for r in K.loci.itertuples():
    w = r.Pos // K.WIN; don = pathd[r.Chrom][w]; s = state(r.SV, r.Target, don)
    sc = state(r.SV, r.Target, pathc[r.Chrom][w])
    pdn = K.PROXY.get(don); sp = state(r.SV, r.Target, pdn) if pdn else s
    out.append(dict(SV=r.SV, Chrom=r.Chrom, Pos=r.Pos, Window_Index=w, Target=r.Target, Selected_Donor=don,
                    Selected_Donor_Carries_ALT="NA" if K.alle.get(r.SV, {}).get(don) is None else int(K.alle[r.SV][don]),
                    Captured=1 if s == 1 else 0,
                    Capture_Status="Captured" if s == 1 else ("Missed" if s == 0 else "Unknown_no_donor_call"),
                    Captured_with_collapsed_proxy="NA" if sp is None else sp,
                    Current_Selected_Donor=pathc[r.Chrom][w], Captured_Current=1 if sc == 1 else 0))
CP = pd.DataFrame(out); CP.to_csv(f"{H}/capture_mask_path.tsv", sep="\t", index=False)
assert CP.Captured.sum() == cap4
print('old-284 path captures', int(CP.Captured_Current.sum()), 'of the 233 loci')
per = []
for c in chroms:
    g = CP[CP.Chrom == c]
    per.append(dict(Chrom=c, Fav_Total=len(g), Captured=int(g.Captured.sum()),
                    Missed=int((g.Capture_Status == "Missed").sum()),
                    Unknown_no_donor_call=int((g.Capture_Status == "Unknown_no_donor_call").sum()),
                    Capture_Rate_pct=round(100 * g.Captured.sum() / len(g), 2) if len(g) else None,
                    Captured_collapsed_proxy=int((g.Captured_with_collapsed_proxy == 1).sum()),
                    Captured_Current=int(g.Captured_Current.sum())))
PC = pd.DataFrame(per)
PC.loc[len(PC)] = dict(Chrom="TOTAL", Fav_Total=len(CP), Captured=int(CP.Captured.sum()),
                       Missed=int((CP.Capture_Status == "Missed").sum()),
                       Unknown_no_donor_call=int((CP.Capture_Status == "Unknown_no_donor_call").sum()),
                       Capture_Rate_pct=round(100 * CP.Captured.sum() / len(CP), 2),
                       Captured_collapsed_proxy=int((CP.Captured_with_collapsed_proxy == 1).sum()),
                       Captured_Current=int(CP.Captured_Current.sum()))
PC.to_csv(f"{H}/per_chrom_capture_mask.tsv", sep="\t", index=False)
mk = {}; mkc = {}
for o in out:
    k = (o["Chrom"], o["Pos"]); mk[k] = max(mk.get(k, 0), o["Captured"]); mkc[k] = max(mkc.get(k, 0), o["Captured_Current"])
with open(f"{H}/position_marks_mask.tsv", "w") as fh:
    fh.write("Chrom\tPos\tMark\tMark_Current\tChanged\n")
    for k in sorted(mk, key=lambda k: (int(k[0][3:5]), k[1])):
        fh.write(f"{k[0]}\t{k[1]}\t{mk[k]}\t{mkc[k]}\t{int(mk[k] != mkc[k])}\n")

diff = [(c, w, donors[pc[c][w]], donors[p4[c][w]], round(float(K.COV[c][w, pc[c][w]]), 4), int(L[c][w, pc[c][w]]),
         int(L[c][w, p4[c][w]])) for c in chroms for w in range(len(p4[c])) if pc[c][w] != p4[c][w]]
with open(f"{H}/path_diff_current_vs_mask.tsv", "w") as fh:
    fh.write("Chrom\tWindow_Index\tDonor_Current\tDonor_Masked\tCoverage_of_Current_Donor\tLoad_Current\tLoad_Masked\n")
    for d in diff: fh.write("\t".join(map(str, d)) + "\n")

# ---- 5i matrix (chr01B): load (final Nigerian calls), coverage, eligibility, selection
z = []
for w in range(K.nwin["chr01B"]):
    for d in donors:
        j = di[d]
        z.append(dict(Donor_ID=d, Donor=K.disp(d), Window_Index=w, Window_Start_0based=w * K.WIN,
                      Window_End_0based=int(K.ends["chr01B"][w]), DSV_Count=int(K.DSV["chr01B"][w, j]),
                      DSNP_Count=int(K.DSNP1["chr01B"][w, j]), Total_Load=int(L["chr01B"][w, j]),
                      Coverage=round(float(K.COV["chr01B"][w, j]), 4), Eligible=int(E["chr01B"][w, j]),
                      Selected_on_path=int(p4["chr01B"][w] == j)))
pd.DataFrame(z).to_csv(f"{H}/zoom_matrix_chr01B.tsv", sep="\t", index=False)

cur_zero = sum(int((L[c][np.arange(len(pc[c])), pc[c]] == 0).sum()) for c in chroms)
new_zero = sum(int((L[c][np.arange(len(p4[c])), p4[c]] == 0).sum()) for c in chroms)
cur_lowcov = sum(int((K.COV[c][np.arange(len(pc[c])), pc[c]] < T).sum()) for c in chroms)
print("windows changed vs current:", len(diff), "| current path windows with selected-donor coverage < 50%:", cur_lowcov,
      "| zero-load windows on path: current", cur_zero, "masked", new_zero)
print("Nigerian windows current", Counter(donors[i] for c in chroms for i in pc[c] if donors[i].startswith('nrly')),
      "masked", Counter(donors[i] for c in chroms for i in p4[c] if donors[i].startswith('nrly')))
print(PC.to_string())
print(BC[["Chrom", "Residual_DSV", "Residual_DSNP", "Breakpoint_Count", "Segment_Count", "Donor_Count", "Objective"]].to_string())
