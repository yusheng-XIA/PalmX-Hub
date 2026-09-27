#!/usr/bin/env python3
"""Coverage-masked African35 donor designs: mask summary, W sweep (P = 15), constrained designs, sensitivity.

Outputs (results/):
  mask_summary_by_donor.tsv      windows below 25/50/75% coverage per donor, and zero-load windows among them
  all_donor_lowcov_windows.tsv   windows in which no donor reaches the threshold (left eligible for all)
  single_donor_masked.tsv        single-donor totals (definition in core_mask.single_donor_loads)
  weight_sweep_mask.tsv          W = 0,1,2,4,8 at P = 15, T = 0.5 (plus the collapsed-proxy reward set)
  constrained_paths_mask.tsv     k <= 1,2,3,5 breakpoints/chrom; m <= 2..5 donors (<= 2 bp/chrom); W = 0, T = 0.5
  sensitivity_mask.tsv           W = 0 and W = 4 paths for T = 0 (no mask), 0.25, 0.5, 0.75, with old/new Nigerian calls
"""
import itertools, sys, time
import numpy as np, pandas as pd
import core_mask as K
from core_mask import chroms, donors, di

OUT = f"{K.H}/results"
import os; os.makedirs(OUT, exist_ok=True)
T0 = 0.5
t0 = time.time()

# ---------- sanity: unmasked, old calls reproduce the current Fig. 5h path (W = 4) and W = 0
E_all = {c: np.ones_like(K.L0[c], bool) for c in chroms}
p, o = K.run_weighted(4, 15, K.F_OLD, K.L0); s = K.stats(p, K.L0)
assert (s[0], s[2], K.capture(p, K.F_OLD), o) == (23530, 349, 202, 27957), (s, o)
p, o = K.run_weighted(0, 15, K.F_OLD, K.L0); s = K.stats(p, K.L0)
assert (s[0], s[2], K.capture(p, K.F_OLD), o) == (23520, 346, 168, 28710), (s, o)
print("unmasked old-284 paths reproduced (W4 23530/349/202; W0 23520/346/168); now using the 233 reported-GWAS loci")

# ---------- mask summary
rows = []
for d in donors:
    j = di[d]; r = dict(Donor_ID=d, Donor=K.disp(d),
                        Coverage_Source="final assembly (re-aligned)" if d.startswith("nrly") else
                        ("original PAF (re-aligned, same command)" if d.startswith(("dura", "pisifera")) else "original PAF"),
                        Mean_Window_Coverage=round(float(np.mean(np.concatenate([K.COV[c][:, j] for c in chroms]))), 4))
    for T in (0.25, 0.5, 0.75):
        E, ALL = K.eligibility(T)
        lo = np.concatenate([~E[c][:, j] for c in chroms])
        z = np.concatenate([K.L1[c][:, j] == 0 for c in chroms])
        r[f"Masked_Windows_cov_lt_{int(T*100)}"] = int(lo.sum())
        r[f"Masked_ZeroLoad_Windows_cov_lt_{int(T*100)}"] = int((lo & z).sum())
        r[f"Masked_Load_cov_lt_{int(T*100)}"] = int(np.concatenate([K.L1[c][:, j] for c in chroms])[lo].sum())
    rows.append(r)
MS = pd.DataFrame(rows); MS.to_csv(f"{OUT}/mask_summary_by_donor.tsv", sep="\t", index=False)
print(MS.to_string())
al = []
for T in (0.25, 0.5, 0.75):
    E, ALL = K.eligibility(T)
    for c in chroms:
        for w in np.where(ALL[c])[0]:
            al.append(dict(Threshold=T, Chrom=c, Window_Index=int(w), Window_Start_0based=int(w) * K.WIN,
                           Max_Donor_Coverage=round(float(K.COV[c][w].max()), 4)))
AL = pd.DataFrame(al); AL.to_csv(f"{OUT}/all_donor_lowcov_windows.tsv", sep="\t", index=False)
print("windows with no eligible donor (left open):", AL.groupby("Threshold").size().to_dict())

# ---------- single donor
E5, _ = K.eligibility(T0)
SD, CI5 = K.single_donor_loads(K.L1, E5)
SDt = pd.DataFrame([dict(Donor_ID=d, Donor=K.disp(d), **v) for d, v in SD.items()]).sort_values("Imputed_Total")
SDt.to_csv(f"{OUT}/single_donor_masked.tsv", sep="\t", index=False)
SDt = SDt.sort_values(["Imputed_Total", "Donor_ID"])
best_single = SDt.iloc[0].Donor_ID; BSL = int(SDt.iloc[0].Imputed_Total)
BSL_i = BSL
print("best single (masked, imputed):", best_single, BSL, "| raw best:",
      min(donors, key=lambda d: (SD[d]["Raw_Total"], d)), min(v["Raw_Total"] for v in SD.values()))
print(SDt.head(8).to_string())

C5 = K.masked_cost(K.L1, E5)

# ---------- W sweep
def pathrow(path, obj, L, E, **kw):
    ld, dsv, bp, mx, used = K.stats(path, L)
    return dict(**kw, Residual_Total_Load=ld, Residual_DSV=dsv, Residual_DSNP=ld - dsv, Breakpoints=bp,
                Segments=bp + len(chroms), Donor_Count=len(used), Objective=obj,
                Fav_Captured=K.capture(path, K.F0), Fav_Captured_collapsed_proxy=K.capture(path, K.Fp),
                Fav_Total=len(K.loci), Ineligible_Windows_Selected=K.n_ineligible(path, E))

sw = []
paths = {}
for F, lab in ((K.F0, "step9_calls_only"), (K.Fp, "collapsed_proxy_for_6_phased_donors")):
    for W in (0, 1, 2, 4, 8):
        path, obj = K.run_weighted(W, K.P0, F, C5)
        r = pathrow(path, obj, K.L1, E5, Panel="African35", Coverage_Threshold=T0, Reward_Call_Set=lab,
                    GWAS_Weight_W=W, Breakpoint_Penalty_P=K.P0)
        r["Reduction_vs_Best_Single_pct"] = round(100 * (BSL - r["Residual_Total_Load"]) / BSL, 2)
        r["Max_Attainable_Capture_this_call_set"] = int(sum((F[c] * E5[c]).max(1).sum() for c in chroms))
        sw.append(r)
        if lab == "step9_calls_only": paths[W] = path
S = pd.DataFrame(sw)
for lab, g in S.groupby("Reward_Call_Set"):
    base = g[g.GWAS_Weight_W == 0].iloc[0]
    S.loc[g.index, "Delta_Load_vs_W0"] = g.Residual_Total_Load - base.Residual_Total_Load
    S.loc[g.index, "Delta_Capture_vs_W0"] = g.Fav_Captured - base.Fav_Captured
S.to_csv(f"{OUT}/weight_sweep_mask.tsv", sep="\t", index=False)
print(S.drop(columns=["Panel"]).to_string())
assert (S.Ineligible_Windows_Selected == 0).all()
np.save(f"{OUT}/path_W4_T50.npy", np.concatenate([paths[4][c] for c in chroms]))
np.save(f"{OUT}/path_W0_T50.npy", np.concatenate([paths[0][c] for c in chroms]))

# ---------- sensitivity (threshold x Nigerian calls)
sens = []
for calls, L in (("original_7-25_assembly", K.L0), ("final_assembly", K.L1)):
    for T in (0.0, 0.25, 0.5, 0.75):
        E, ALL = K.eligibility(T) if T > 0 else (E_all, {c: np.zeros(K.L0[c].shape[0], bool) for c in chroms})
        C = K.masked_cost(L, E)
        sdl, _ = K.single_donor_loads(L, E); bsd = min(donors, key=lambda d: (sdl[d]["Imputed_Total"], d))
        for W in (0, 4):
            path, obj = K.run_weighted(W, K.P0, K.F0, C)
            r = pathrow(path, obj, L, E, Nigerian_dSNP_Calls=calls, Coverage_Threshold=T, GWAS_Weight_W=W,
                        Breakpoint_Penalty_P=K.P0, Masking="excluded" if T > 0 else "none")
            r.update(Best_Single_Donor=bsd, Best_Single_Donor_Load=sdl[bsd]["Imputed_Total"],
                     Reduction_vs_Best_Single_pct=round(100 * (sdl[bsd]["Imputed_Total"] - r["Residual_Total_Load"]) / sdl[bsd]["Imputed_Total"], 2),
                     Windows_No_Eligible_Donor=int(sum(ALL[c].sum() for c in chroms)),
                     Nigerian_Windows=int(sum(np.isin(path[c], [di["nrly_hap1"], di["nrly_hap2"]]).sum() for c in chroms)),
                     Nigerian_Hap2_Windows=int(sum((path[c] == di["nrly_hap2"]).sum() for c in chroms)),
                     ZeroLoad_Windows_On_Path=int(sum((L[c][np.arange(len(path[c])), path[c]] == 0).sum() for c in chroms)))
            sens.append(r)
E, ALL = K.eligibility(T0); sdl, CI = K.single_donor_loads(K.L1, E); bsd = min(donors, key=lambda d: (sdl[d]["Imputed_Total"], d))
for W in (0, 4):
    path, obj = K.run_weighted(W, K.P0, K.F0, CI)
    r = pathrow(path, obj, CI, E, Nigerian_dSNP_Calls="final_assembly", Coverage_Threshold=T0, GWAS_Weight_W=W,
                Breakpoint_Penalty_P=K.P0)
    r.update(Masking="median-imputed instead of excluded", Best_Single_Donor=bsd, Best_Single_Donor_Load=sdl[bsd]["Imputed_Total"],
             Reduction_vs_Best_Single_pct=round(100 * (sdl[bsd]["Imputed_Total"] - r["Residual_Total_Load"]) / sdl[bsd]["Imputed_Total"], 2))
    sens.append(r)
SE = pd.DataFrame(sens); SE.to_csv(f"{OUT}/sensitivity_mask.tsv", sep="\t", index=False)
print(SE.to_string())

# ---------- constrained designs (W = 0, T = 0.5)
P15, _ = K.run_weighted(0, K.P0, K.F0, C5); U = K.stats(P15, K.L1)[0]
crow = []
def add(label, cls, path, note="", Lc=None):
    ld, dsv, bp, mx, used = K.stats(path, CI5 if Lc is None else Lc)
    inel = K.n_ineligible(path, E5)
    cnt = pd.Series([donors[i] for c in chroms for i in path[c]]).value_counts()
    crow.append(dict(Panel="African35", Design=label, Best_Single_Donor_Load=BSL_i, Constraint_Class=cls,
                     Residual_Total_Load=ld, Residual_DSV=dsv, Residual_DSNP=ld - dsv,
                     Reduction_vs_Best_Single=round(BSL - ld, 1),
                     Reduction_vs_Best_Single_pct=round(100 * (BSL - ld) / BSL, 2),
                     Share_of_Unconstrained_P15_Reduction_pct=round(100 * (BSL - ld) / (BSL - U), 1),
                     Total_Breakpoints=bp, Max_Breakpoints_per_Chrom=mx, Donor_Count=len(used),
                     Donors_by_windows=";".join(f"{d}:{n}" for d, n in cnt.items()),
                     Fav_Captured=K.capture(path, K.F0), Fav_Captured_collapsed_proxy=K.capture(path, K.Fp),
                     Fav_Total=len(K.loci), Ineligible_Windows_Selected=inel, Coverage_Threshold=T0, Note=note))

bj = di[best_single]
crow.append(dict(Panel="African35", Design=f"Best single donor ({best_single})", Best_Single_Donor_Load=BSL_i,
                 Constraint_Class="reference", Residual_Total_Load=BSL_i,
                 Residual_DSV=int(sum(K.DSV[c][E5[c][:, bj], bj].sum() for c in chroms)),
                 Residual_DSNP=BSL_i - int(sum(K.DSV[c][E5[c][:, bj], bj].sum() for c in chroms)),
                 Reduction_vs_Best_Single=0, Reduction_vs_Best_Single_pct=0.0,
                 Share_of_Unconstrained_P15_Reduction_pct=0.0, Total_Breakpoints=0, Max_Breakpoints_per_Chrom=0,
                 Donor_Count=1, Donors_by_windows=f"{best_single}:{sum(K.nwin.values())}",
                 Fav_Captured=K.capture({c: np.full(K.nwin[c], bj) for c in chroms}), Fav_Captured_collapsed_proxy=None,
                 Fav_Total=len(K.loci), Ineligible_Windows_Selected=SD[best_single]["Imputed_Windows"],
                 Coverage_Threshold=T0,
                 Note=f"own load in {sum(K.nwin.values()) - SD[best_single]['Imputed_Windows']} eligible windows "
                      f"({SD[best_single]['Own_Eligible_Load']}) + median eligible-donor load in "
                      f"{SD[best_single]['Imputed_Windows']} windows below the coverage threshold"))
for Kb in (1, 2, 3, 5):
    add(f"<= {Kb} breakpoint(s)/chrom, donors unrestricted", "a_per_chrom_breakpoints",
        {c: K.kbreak_dp(CI5[c], Kb) for c in chroms}, "" if Kb in (1, 2) else "extra point for curve")
    print("k", Kb, crow[-1]["Residual_Total_Load"], crow[-1]["Ineligible_Windows_Selected"], round(time.time() - t0)); sys.stdout.flush()

nd = len(donors)
trip = np.array(list(itertools.combinations(range(nd), 3)))
pair = np.array(list(itertools.combinations(range(nd), 2)))
ftrip = np.stack([K.subset_dp_values(CI5[c], trip) for c in chroms])
fpair = np.stack([K.subset_dp_values(CI5[c], pair) for c in chroms])
tidx = {tuple(t): i for i, t in enumerate(trip)}
best = {}
tot2 = fpair.sum(0); i = int(np.argmin(tot2)); best[2] = (int(tot2[i]), tuple(pair[i]))
tot3 = ftrip.sum(0); i = int(np.argmin(tot3)); best[3] = (int(tot3[i]), tuple(trip[i]))
for msz in (4, 5):
    bestv, bests = None, None
    sub3 = list(itertools.combinations(range(msz), 3))
    chunk = []
    def flush(chunk):
        global bestv, bests
        A = np.array(chunk)
        idx = np.array([[tidx[tuple(A[j, list(t)])] for t in sub3] for j in range(len(A))])
        vals = ftrip[:, idx].min(2).sum(0)
        j = int(np.argmin(vals))
        if bestv is None or vals[j] < bestv:
            bestv, bests = int(vals[j]), tuple(A[j])
    for sset in itertools.combinations(range(nd), msz):
        chunk.append(sset)
        if len(chunk) == 20000:
            flush(chunk); chunk = []
    if chunk: flush(chunk)
    best[msz] = (bestv, bests)
    print("m", msz, bestv, [donors[i] for i in bests], round(time.time() - t0)); sys.stdout.flush()
for msz in (2, 3, 4, 5):
    v, sub = best[msz]
    path = {c: K.subset_dp_path(CI5[c], list(sub)) for c in chroms}
    assert K.stats(path, CI5)[0] == v, (msz, K.stats(path, CI5)[0], v)
    add(f"<= {msz} donors genome-wide, <= 2 breakpoints/chrom", "b_donor_budget", path,
        "Optimal donor set: " + ",".join(donors[i] for i in sub) + ("" if msz in (2, 3, 5) else "; extra point for curve"))
add("Unconstrained IPH (W = 0, P = 15)", "reference", P15, "W = 0 comparison path", Lc=K.L1)
p4, o4 = K.run_weighted(4, K.P0, K.F0, C5)
add("Unconstrained IPH (W = 4, P = 15)", "reference", p4, "Fig. 5h path", Lc=K.L1)
CP = pd.DataFrame(crow); CP.to_csv(f"{OUT}/constrained_paths_mask.tsv", sep="\t", index=False)
print(CP.drop(columns=["Donors_by_windows"]).to_string())
print("done", round(time.time() - t0))
