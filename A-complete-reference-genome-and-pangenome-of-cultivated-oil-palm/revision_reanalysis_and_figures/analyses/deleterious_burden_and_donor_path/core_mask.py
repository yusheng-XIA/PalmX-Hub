#!/usr/bin/env python3
"""Shared data and dynamic-programming core for the coverage-masked African35 donor path (Fig. 5h,i).

Load matrix
    results-8.9 09_ideal_parent_haplotypes/african35/load_matrix_500kb.tsv (copy: ../enh_C/data/), with the
    DSNP_Count of nrly_hap1 / nrly_hap2 replaced by the re-call on the final Nigerian assemblies
    (data/dsnp/nrly_final_hap{1,2}.windows.tsv; same minimap2/paftools commands and the same stage 05/09 rules,
    scripts/10_call_dsnp_asm.sh + scripts/21_dsnp_windows.py). DSV_Count is unchanged (formal African35 dSV set).
Coverage mask
    data/coverage/<donor>.tsv: fraction of each 500-kb reference window covered by the donor's alignments
    in the PAF that produced its dSNP calls (paftools call filters: aligned length >= 1000, MAPQ >= 5,
    primary only). A donor is not eligible in a window with coverage < T (primary T = 0.5; 0.25 and 0.75
    as sensitivity). Windows in which no donor reaches T are left eligible for every donor (recorded).
Objective (71_build_loadonly_iph.py::exact_dp; enh_C.py): sum_w [L(w,d_w) - W*F(w,d_w)] + P*switches, with
    ineligible (window, donor) cells given cost BIG (never chosen while an eligible donor exists).
Functions exact_dp, kbreak_dp, subset_dp_values, subset_dp_path, fav_matrix are copied verbatim from
../enh_C/enh_C.py (tie-breaking identical).
"""
import os
import numpy as np
import pandas as pd

H = os.path.dirname(os.path.abspath(__file__))
FIX = os.path.dirname(H)
WIN = 500_000
P0 = 15
BIG = 10 ** 9
PROXY = {"dura_hap1": "EG_dura", "dura_hap2": "EG_dura", "pisifera_hap1": "EG_pisifera",
         "pisifera_hap2": "EG_pisifera", "nrly_hap1": "EG_niriliya", "nrly_hap2": "EG_niriliya"}

m0 = pd.read_csv(f"{FIX}/enh_C/data/load_matrix_500kb.tsv", sep="\t")
donors = sorted(m0.Sample_ID.unique())
di = {d: i for i, d in enumerate(donors)}
chroms = sorted(m0.Chrom.unique(), key=lambda c: (int(c[3:5]), c))
nwin = {c: int(m0[m0.Chrom == c].Window_Index.max()) + 1 for c in chroms}
ends = {c: m0[(m0.Chrom == c)].groupby("Window_Index").Window_End_0based.first().to_numpy() for c in chroms}


def pivot(m, col):
    return {c: m[m.Chrom == c].pivot(index="Window_Index", columns="Sample_ID", values=col)[donors].to_numpy(np.int64)
            for c in chroms}


# ---- original (unmasked, old Nigerian calls) matrix
L0, DSV, DSNP0 = pivot(m0, "Total_Load"), pivot(m0, "DSV_Count"), pivot(m0, "DSNP_Count")

# ---- Nigerian dSNP re-called on the final assemblies
m1 = m0.copy()
NEWCALL = {}
for d in ("nrly_hap1", "nrly_hap2"):
    n = pd.read_csv(f"{H}/data/dsnp/nrly_final_{d[-4:]}.windows.tsv", sep="\t")
    NEWCALL[d] = n
    key = m1.Sample_ID == d
    mm = m1[key].merge(n, on=["Chrom", "Window_Index"], suffixes=("", "_new"), how="left")
    assert len(mm) == key.sum() and mm.DSNP_Count_new.notna().all()
    m1.loc[key, "DSNP_Count"] = mm.DSNP_Count_new.astype(int).to_numpy()
m1["Total_Load"] = m1.DSV_Count + m1.DSNP_Count
L1, DSNP1 = pivot(m1, "Total_Load"), pivot(m1, "DSNP_Count")

# ---- coverage (Nigerian: final assemblies)
COVFILE = {d: f"{H}/data/coverage/{'nrly_final_' + d[-4:] if d.startswith('nrly') else d}.tsv" for d in donors}
cv = []
for d in donors:
    c = pd.read_csv(COVFILE[d], sep="\t"); c["Donor"] = d; cv.append(c)
cv = pd.concat(cv)
COV = {c: cv[cv.Chrom == c].pivot(index="Window_Index", columns="Donor", values="Coverage")[donors].to_numpy(float)
       for c in chroms}
for c in chroms:
    assert COV[c].shape == L0[c].shape


def eligibility(T):
    """E[c] boolean windows x donors; ALL[c] boolean windows where no donor reaches T (then all eligible)."""
    E, ALL = {}, {}
    for c in chroms:
        e = COV[c] >= T
        none = ~e.any(1)
        e[none] = True
        E[c], ALL[c] = e, none
    return E, ALL


def masked_cost(L, E):
    return {c: np.where(E[c], L[c], BIG) for c in chroms}


# ---------- favourable loci (enh_C.py)
alle = {}
for r in pd.read_csv(f"{FIX}/fig5hi/src/step9_donor_alleles.tsv", sep="\t").itertuples():
    alle.setdefault(r.sv_id, {})[r.donor] = r.carries_ALT == 1
loci = pd.read_csv(f"{FIX}/misc9/data/ideal_step9_ALL_loci.tsv", sep="\t")


def fav_matrix(proxy=False):
    F = {c: np.zeros_like(L0[c]) for c in chroms}
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


F0 = fav_matrix(False)
Fp = fav_matrix(True)


def capture(path, F=F0):
    return int(sum(F[c][np.arange(len(path[c])), path[c]].sum() for c in chroms))


def stats(path, L):
    load = sum(int(L[c][np.arange(len(path[c])), path[c]].sum()) for c in chroms)
    dsv = sum(int(DSV[c][np.arange(len(path[c])), path[c]].sum()) for c in chroms)
    bp_by = {c: int((np.diff(path[c]) != 0).sum()) for c in chroms}
    used = sorted({donors[i] for c in chroms for i in np.unique(path[c])})
    return load, dsv, sum(bp_by.values()), max(bp_by.values()), used


def n_ineligible(path, E):
    return sum(int((~E[c][np.arange(len(path[c])), path[c]]).sum()) for c in chroms)


# ---------- DP (verbatim from enh_C.py)
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


def run_weighted(W, P, F, C):
    path = {}; obj = 0
    for c in chroms:
        path[c], o = exact_dp(C[c] - W * F[c], P); obj += o
    return path, obj


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
            srt = np.sort(prevb, 1); m1_, m2_ = srt[:, :1], srt[:, 1:2]
            excl = np.where(prevb == m1_, np.where((prevb == m1_).sum(1, keepdims=True) > 1, m1_, m2_), m1_)
            NV[:, b, :] = np.minimum(V[:, b, :], excl)
        V = NV + C[w][:, None, :]
    return V.reshape(S, -1).min(1)


def subset_dp_path(cost, sub, K=2):
    c2 = np.full_like(cost, 10**12); c2[:, sub] = cost[:, sub]
    return kbreak_dp(c2, K)


def disp(d):
    import re
    m = re.fullmatch(r"(dura|pisifera|nrly|bk)_hap([12])", d)
    return f"{dict(dura='TK', pisifera='NS', nrly='Nigerian', bk='TN')[m.group(1)]}-Hap{m.group(2)}" if m else d


def lower_median(v):
    v = np.sort(v); return int(v[(len(v) - 1) // 2])


def window_medians(L, E):
    """per window: lower median of the loads of the donors eligible in that window (an observed load value)."""
    return {c: np.array([lower_median(L[c][w][E[c][w]]) for w in range(L[c].shape[0])], np.int64) for c in chroms}


def imputed_cost(L, E):
    """cost for designs with hard breakpoint/donor limits (and the single-donor reference): a donor's own load where
    it is eligible, the lower-median eligible-donor load where it is not."""
    med = window_medians(L, E)
    return {c: np.where(E[c], L[c], med[c][:, None]) for c in chroms}


def single_donor_loads(L, E):
    """Single-donor totals under the mask: own load in eligible windows + lower-median eligible-donor load in each
    window where the donor is below the coverage threshold (neutral imputation)."""
    out = {}
    C = imputed_cost(L, E)
    for d in donors:
        j = di[d]
        own = int(sum(L[c][E[c][:, j], j].sum() for c in chroms))
        out[d] = dict(Imputed_Total=int(sum(C[c][:, j].sum() for c in chroms)), Own_Eligible_Load=own,
                      Imputed_Windows=int(sum((~E[c][:, j]).sum() for c in chroms)),
                      Raw_Total=int(sum(L[c][:, j].sum() for c in chroms)))
    return out, C


# ======================= fav233 (2026-09-26): favourable loci re-derived from the reported SV-GWAS =======================
# ../fav_rederive/: 60 retained traits, tail-trimmed phenotypes, the author's recommendation rules
# (TRAIT_DIRECTION, MIN_RECOMMEND_N = 5, genotype-consistent set, homozygous target, REF/ALT of unequal length) -> 233 loci.
OLD_loci, OLD_alle = loci, alle
F_OLD, Fp_OLD = F0, Fp
alle = {}
for r in pd.read_csv(f"{FIX}/fav_rederive/data/donor_alleles.tsv", sep="\t").itertuples():
    alle.setdefault(r.sv_id, {})[r.donor] = r.carries_ALT == 1
loci = pd.read_csv(f"{FIX}/fav_rederive/data/fav_loci_reported.tsv", sep="\t")[["SV", "Chrom", "Pos", "Target"]]
F0 = fav_matrix(False)
Fp = fav_matrix(True)
NFAV = len(loci)
assert NFAV == 233


def capture(path, F=None):
    F = F0 if F is None else F
    return int(sum(F[c][np.arange(len(path[c])), path[c]].sum() for c in chroms))
