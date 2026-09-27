#!/usr/bin/env python3
"""Enhancement C: feasibility-constrained IPH designs and GWAS-weight sweep on the African35 panel.
Inputs (read-only copies): results-8.9/09_ideal_parent_haplotypes/african35/load_matrix_500kb.tsv,
ideal_loadonly_path.tsv, penalty_sweep.tsv; step9_donor_alleles.tsv; ideal_step9_ALL_loci.tsv (284 favourable loci).
Objective (original 71_build_loadonly_iph.py::exact_dp): sum(dSV+dSNP) - W*F + P*switches, per chromosome, 500-kb windows.
F[w,d] = number of favourable loci in window w for which donor d carries the target allele
(compute_capture.py definition; donor without a step9 call -> not captured).
"""
import itertools, os, sys
import numpy as np, pandas as pd

H = os.path.dirname(os.path.abspath(__file__))
FIX = os.path.dirname(H)
D = os.path.join(H, "data")
WIN = 500_000
P0 = 15
PROXY = {"dura_hap1": "EG_dura", "dura_hap2": "EG_dura", "pisifera_hap1": "EG_pisifera",
         "pisifera_hap2": "EG_pisifera", "nrly_hap1": "EG_niriliya", "nrly_hap2": "EG_niriliya"}

m = pd.read_csv(f"{D}/load_matrix_500kb.tsv", sep="\t")
donors = sorted(m.Sample_ID.unique())                      # index order = name order (tie-break as original)
di = {d: i for i, d in enumerate(donors)}
chroms = sorted(m.Chrom.unique(), key=lambda c: (int(c[3:5]), c))
L = {}; DSV = {}
for c in chroms:
    s = m[m.Chrom == c]
    L[c] = s.pivot(index="Window_Index", columns="Sample_ID", values="Total_Load")[donors].to_numpy(np.int64)
    DSV[c] = s.pivot(index="Window_Index", columns="Sample_ID", values="DSV_Count")[donors].to_numpy(np.int64)
    assert (np.diff(s.Window_Index.drop_duplicates().sort_values().to_numpy()) == 1).all()
nW = sum(L[c].shape[0] for c in chroms)
single = {d: sum(int(L[c][:, di[d]].sum()) for c in chroms) for d in donors}
best_single = min(donors, key=lambda d: (single[d], d)); BSL = single[best_single]
print("windows", nW, "best single", best_single, BSL)

# ---------- favourable loci ----------
alle = {}
for r in pd.read_csv(f"{FIX}/fig5hi/src/step9_donor_alleles.tsv", sep="\t").itertuples():
    alle.setdefault(r.sv_id, {})[r.donor] = r.carries_ALT == 1
loci = pd.read_csv(f"{FIX}/misc9/data/ideal_step9_ALL_loci.tsv", sep="\t")

def fav_matrix(proxy=False):
    F = {c: np.zeros_like(L[c]) for c in chroms}
    for r in loci.itertuples():
        w = r.Pos // WIN
        if r.Chrom not in F or w >= F[r.Chrom].shape[0]:
            continue
        for d in donors:
            ca = alle.get(r.SV, {}).get(PROXY.get(d, d) if proxy else d)
            if ca is None:
                continue
            if (r.Target == "ALT" and ca) or (r.Target == "REF" and not ca):
                F[r.Chrom][w, di[d]] += 1
    return F
F0 = fav_matrix(False); Fp = fav_matrix(True)

def capture(path, F):
    return int(sum(F[c][np.arange(len(path[c])), path[c]].sum() for c in chroms))

def stats(path):
    load = sum(int(L[c][np.arange(len(path[c])), path[c]].sum()) for c in chroms)
    dsv = sum(int(DSV[c][np.arange(len(path[c])), path[c]].sum()) for c in chroms)
    bp_by = {c: int((np.diff(path[c]) != 0).sum()) for c in chroms}
    used = sorted({donors[i] for c in chroms for i in np.unique(path[c])})
    return load, dsv, sum(bp_by.values()), max(bp_by.values()), used

# ---------- 1. unconstrained penalised DP (exact replica of exact_dp tie-breaking) ----------
def exact_dp(cost, pen):
    n, k = cost.shape
    prev = cost[0].astype(np.int64).copy()
    back = np.zeros((n, k), np.int64)
    for w in range(1, n):
        order = np.argsort(prev, kind="stable")           # lowest cost, then lowest name
        b1, b2 = order[0], order[1]
        sw_src = np.where(np.arange(k) == b1, b2, b1)     # best other donor
        sw_cost = prev[sw_src] + pen
        stay = prev <= sw_cost                            # stay wins ties
        back[w] = np.where(stay, np.arange(k), sw_src)
        prev = np.where(stay, prev, sw_cost) + cost[w]
    end = int(np.lexsort((np.arange(k), prev))[0])
    path = np.zeros(n, np.int64); path[-1] = end
    for w in range(n - 1, 0, -1):
        path[w - 1] = back[w][path[w]]
    return path, int(prev[end])

def run_weighted(W, P, F):
    path = {}; obj = 0
    for c in chroms:
        path[c], o = exact_dp(L[c] - W * F[c], P); obj += o
    return path, obj

p15, obj15 = run_weighted(0, P0, F0)
ld, dsv, bp, mx, used = stats(p15)
ref = pd.read_csv(f"{D}/ideal_loadonly_path.tsv", sep="\t" if True else None)
ident = sum(int(donors[p15[r.Chrom][r.Window_Index]] == r.Donor_ID) for r in ref.itertuples())
repro = dict(Load=ld, DSV=dsv, DSNP=ld - dsv, Breakpoints=bp, Segments=bp + len(chroms), Donors=len(used),
             Objective=obj15, Identical_windows=ident, Total_windows=len(ref), Capture=capture(p15, F0),
             Capture_proxy=capture(p15, Fp), Reduction_pct=100 * (BSL - ld) / BSL)
print("REPRO", repro)
ps = pd.read_csv(f"{D}/penalty_sweep.tsv", sep="\t")
psc = []
for P in ps.Breakpoint_Penalty_P:
    pp, oo = run_weighted(0, int(P), F0); l_, _, b_, _, u_ = stats(pp)
    r = ps[ps.Breakpoint_Penalty_P == P].iloc[0]
    psc.append((int(P), l_, b_, oo, int(r.Residual_Total_Load), int(r.Breakpoint_Count), int(r.Objective)))
print("penalty sweep recheck (P, load, bp, obj | ref load, bp, obj):"); [print(" ", x) for x in psc]

# ---------- 2a. <=k breakpoints per chromosome, donors unrestricted (min load) ----------
def kbreak_dp(cost, K):
    """exact min-load path with at most K switches; returns path."""
    n, k = cost.shape
    INF = np.iinfo(np.int64).max // 4
    V = np.full((K + 1, k), INF, np.int64); V[0] = cost[0]
    back = np.zeros((n, K + 1, k, 2), np.int64)           # (prev donor, prev b)
    for w in range(1, n):
        NV = np.full_like(V, INF)
        for b in range(K + 1):
            stay = V[b]
            if b > 0:
                o = np.argsort(V[b - 1], kind="stable"); b1, b2 = o[0], o[1]
                src = np.where(np.arange(k) == b1, b2, b1); sw = V[b - 1][src]
            else:
                src = np.zeros(k, np.int64); sw = np.full(k, INF)
            use_stay = stay <= sw
            NV[b] = np.where(use_stay, stay, sw) + cost[w]
            back[w, b, :, 0] = np.where(use_stay, np.arange(k), src)
            back[w, b, :, 1] = np.where(use_stay, b, b - 1)
        V = NV
    b, d = np.unravel_index(np.argmin(V), V.shape)
    path = np.zeros(n, np.int64); path[-1] = d
    for w in range(n - 1, 0, -1):
        d, b = back[w, b, d]; path[w - 1] = d
    return path

rows = []
def add(label, cls, path, note=""):
    ld, dsv, bp, mx, used = stats(path)
    cnt = pd.Series([donors[i] for c in chroms for i in path[c]]).value_counts()
    rows.append(dict(Design=label, Constraint_Class=cls, Residual_Total_Load=ld, Residual_DSV=dsv,
                     Residual_DSNP=ld - dsv, Reduction_vs_Best_Single=BSL - ld,
                     Reduction_vs_Best_Single_pct=round(100 * (BSL - ld) / BSL, 2),
                     Share_of_Unconstrained_P15_Reduction_pct=round(100 * (BSL - ld) / (BSL - repro["Load"]), 1),
                     Total_Breakpoints=bp, Max_Breakpoints_per_Chrom=mx, Donor_Count=len(used),
                     Donors_by_windows=";".join(f"{d}:{n}" for d, n in cnt.items()),
                     Fav_Captured=capture(path, F0), Fav_Captured_collapsed_proxy=capture(path, Fp),
                     Fav_Total=len(loci), Note=note))

add("Best single donor (EG_107)", "reference", {c: np.full(L[c].shape[0], di[best_single]) for c in chroms})
for K in (1, 2, 3, 5):
    add(f"<= {K} breakpoint(s)/chrom, donors unrestricted", "a_per_chrom_breakpoints",
        {c: kbreak_dp(L[c], K) for c in chroms}, "" if K in (1, 2) else "extra point for curve")

# ---------- 2b. <= m donors genome-wide, <= 2 breakpoints/chrom ----------
def subset_dp_values(cost, subsets, K=2):
    """vectorised exact DP: min load with <=K switches using only donors in each subset (rows of `subsets`)."""
    S, s = subsets.shape
    INF = np.iinfo(np.int64).max // 4
    C = cost[:, subsets]                                   # n x S x s
    V = np.full((S, K + 1, s), INF, np.int64); V[:, 0, :] = C[0]
    for w in range(1, C.shape[0]):
        NV = np.empty_like(V)
        NV[:, 0, :] = V[:, 0, :]
        for b in range(1, K + 1):
            prevb = V[:, b - 1, :]
            # best other donor within subset
            tot_min = prevb.min(1, keepdims=True)
            # min excluding self: compute via sorted two smallest
            srt = np.sort(prevb, 1); m1, m2 = srt[:, :1], srt[:, 1:2]
            excl = np.where(prevb == m1, np.where((prevb == m1).sum(1, keepdims=True) > 1, m1, m2), m1)
            NV[:, b, :] = np.minimum(V[:, b, :], excl)
        V = NV + C[w][:, None, :]
    return V.reshape(S, -1).min(1)

def subset_dp_path(cost, sub, K=2):
    c2 = np.full_like(cost, 10**9); c2[:, sub] = cost[:, sub]
    return kbreak_dp(c2, K)

nd = len(donors)
trip = np.array(list(itertools.combinations(range(nd), 3)))
pair = np.array(list(itertools.combinations(range(nd), 2)))
ftrip = np.stack([subset_dp_values(L[c], trip) for c in chroms])   # chrom x triples
fpair = np.stack([subset_dp_values(L[c], pair) for c in chroms])
tidx = {tuple(t): i for i, t in enumerate(trip)}
best = {}
# m=1
best[1] = min(((single[d], (di[d],)) for d in donors))
# m=2
tot2 = fpair.sum(0); i = int(np.argmin(tot2)); best[2] = (int(tot2[i]), tuple(pair[i]))
# m=3
tot3 = ftrip.sum(0); i = int(np.argmin(tot3)); best[3] = (int(tot3[i]), tuple(trip[i]))
# m=4, m=5 : cost(S) = sum_c min over 3-subsets T of S of f_c(T)  (a <=2-breakpoint path uses <=3 donors)
for msz in (4, 5):
    bestv, bests = None, None
    subs = itertools.combinations(range(nd), msz)
    sub3 = list(itertools.combinations(range(msz), 3))
    chunk = []
    def flush(chunk):
        global bestv, bests
        A = np.array(chunk)                                   # N x msz
        idx = np.array([[tidx[tuple(A[j, list(t)])] for t in sub3] for j in range(len(A))])
        vals = ftrip[:, idx].min(2).sum(0)                    # chrom x N x 10 -> N
        j = int(np.argmin(vals))
        if bestv is None or vals[j] < bestv:
            bestv, bests = int(vals[j]), tuple(A[j])
    for sset in subs:
        chunk.append(sset)
        if len(chunk) == 20000:
            flush(chunk); chunk = []
    if chunk: flush(chunk)
    best[msz] = (bestv, bests)
    print("m", msz, bestv, [donors[i] for i in bests]); sys.stdout.flush()

for msz in (2, 3, 4, 5):
    v, sub = best[msz]
    path = {c: subset_dp_path(L[c], list(sub)) for c in chroms}
    assert stats(path)[0] == v, (msz, stats(path)[0], v)
    add(f"<= {msz} donors genome-wide, <= 2 breakpoints/chrom", "b_donor_budget",
        path, "Optimal donor set: " + ",".join(donors[i] for i in sub) + ("" if msz in (2, 3, 5) else "; extra point for curve"))

add("Unconstrained IPH (W = 0, P = 15)", "reference", p15, "current Fig. 5 path")
C = pd.DataFrame(rows)
C.insert(0, "Panel", "African35"); C.insert(2, "Best_Single_Donor_Load", BSL)
C.to_csv(f"{H}/constrained_paths.tsv", sep="\t", index=False)
print(C.drop(columns=["Donors_by_windows"]).to_string())


# ---------- 3. W sweep ----------
sw = []
for F, lab in ((F0, "step9_calls_only"), (Fp, "collapsed_proxy_for_6_phased_donors")):
    for W in (0, 1, 2, 4, 8):
        path, obj = run_weighted(W, P0, F)
        ld, dsv, bp, mx, used = stats(path)
        sw.append(dict(Panel="African35", Reward_Call_Set=lab, GWAS_Weight_W=W, Breakpoint_Penalty_P=P0,
                       Residual_Total_Load=ld, Residual_DSV=dsv, Residual_DSNP=ld - dsv, Breakpoints=bp,
                       Segments=bp + len(chroms), Donor_Count=len(used), Objective=obj,
                       Reduction_vs_Best_Single_pct=round(100 * (BSL - ld) / BSL, 2),
                       Fav_Captured=capture(path, F0), Fav_Captured_collapsed_proxy=capture(path, Fp),
                       Fav_Total=len(loci),
                       Max_Attainable_Capture_this_call_set=int(sum(F[c].max(1).sum() for c in chroms))))
S = pd.DataFrame(sw)
for lab, g in S.groupby("Reward_Call_Set"):
    base = g[g.GWAS_Weight_W == 0].iloc[0]
    S.loc[g.index, "Delta_Load_vs_W0"] = g.Residual_Total_Load - base.Residual_Total_Load
    S.loc[g.index, "Delta_Capture_vs_W0"] = g.Fav_Captured - base.Fav_Captured
    S.loc[g.index, "Delta_Capture_proxy_vs_W0"] = g.Fav_Captured_collapsed_proxy - base.Fav_Captured_collapsed_proxy
S.to_csv(f"{H}/weight_sweep.tsv", sep="\t", index=False)
print(S.to_string())
print("REPRO", repro)
print("upper bound capture (per-window best donor): calls", int(sum(F0[c].max(1).sum() for c in chroms)), "proxy", int(sum(Fp[c].max(1).sum() for c in chroms)))
