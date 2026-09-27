#!/usr/bin/env python3
"""Fig. 3g-j robustness (TN parent-hybrid expression modes and cis/trans proxy classes).

Rebuilds the published classification from world-readable upstream inputs (same rules as
RUN-TN-EXPR-HETEROSIS-V2-001 / RUN-ASE-HETEROSIS-DOWNSTREAM-V2-001, validated in trace/Fig3B):
  3g  12 modes + Conserved from legacy size-factor-normalised TK, NS, TN (max >= 10, relative tolerance 0.15)
  3h  7 classes: A = log2((TK+1)/(NS+1)), B = TN pooled allelic log2((A+0.5)/(B+0.5)) over qualifying replicates;
      cut-offs |A| >= 0.75, |B| >= 0.5, |A-B| >= 0.75; ASE-eligible gene-stages with max expression >= 10
  3i  cis contribution |B|/(|B|+|A-B|), medians by window x |A| bin
  3j  PDO/DO/ODO x 7 classes, Pearson residuals
Analyses:
  (1) read-depth downsampling (binomial thinning) of parents (TK, NS legacy libraries), hybrid (TN legacy library
      and TN per-replicate allele fragments) or both, to 50% and 75%, 10 repeats each; size factors re-estimated
      on all 76 legacy libraries as in the original pipeline; eligibility re-applied.
  (2) temporal consistency: adjacent sampled stages as pseudo-replicates (independent libraries), agreement vs
      within-pair label permutation (1,000 permutations).
  (3) TN biological replicates: per-replicate allelic ratios (>= 10 informative fragments), replicate correlation,
      direction agreement, and per-replicate reclassification of 3h/3i/3j.
"""
import numpy as np, pandas as pd, sys, time
from pathlib import Path

A3 = Path("${ANALYSIS_DIR}/22_answer_reviews/00_ms/03_V3/03_figure3")
EXPR = A3 / "01_ASE/00_shared/legacy_expression_normalized_long.tsv"
SF = A3 / "01_ASE/00_shared/legacy_expression_size_factors.tsv"
BK = A3 / "01_ASE/01_bk/ase/config/validated_one_to_one_gene_pairs.tsv"
SL = A3 / "01_ASE/02_seedless/ase/config/validated_one_to_one_gene_pairs.tsv"
UB = A3 / "05_multiomics_integration/runs/RUN-MULTIOMICS-ALLELE-CNS-V4-20260723-001/outputs/stage3_existing_ASE_unification_attempt001"
ASSIGN = Path("${CLUSTER_WORK}/trace/Fig3A/out/s11_TN_sample_stage_assignment.tsv")
REF3G = Path("${CLUSTER_WORK}/trace/Fig3B/out/recalc_3g_modes.tsv")
REF3H = Path("${CLUSTER_WORK}/trace/Fig3B/out/recalc_3h_classes_v2.tsv")
OUT = Path("${CLUSTER_WORK}/enh_3fgj/out")
NREP = int(sys.argv[1]) if len(sys.argv) > 1 else 10
STAGES = ["0d", "15d", "35d", "50d", "65d", "80d", "95d", "110d", "125d", "140d", "155d", "170d", "185d",
          "12h", "24h", "36h", "48h", "60h", "72h"]
SIDX = {s: i for i, s in enumerate(STAGES)}
PH = {**{s: "Days 0–65" for s in STAGES[:5]}, **{s: "Days 80–140" for s in STAGES[5:10]},
      **{s: "Days 155–185" for s in STAGES[10:13]}, **{s: "Hours 12–72" for s in STAGES[13:]}}
PHASES = ["Days 0–65", "Days 80–140", "Days 155–185", "Hours 12–72"]
BINS = ["0–1", "1–2", "2–3", "3–4", "4+"]
MODES = [f"M{i}" for i in range(1, 13)] + ["Conserved"]
REG = ["I.Cis_only", "II.Trans_only", "III.Cis_trans_enhancing", "IV.Cis_trans_compensating", "V.Compensatory",
       "VI.Conserved", "VII.Ambiguous"]
INH = ["PDO", "DO", "ODO"]
rng = np.random.default_rng(20260924)
log = open(OUT / "s02_log.txt", "w")
T0 = time.time()
def P(*a):
    s = " ".join(str(x) for x in a); print(f"[{time.time()-T0:7.1f}s] {s}", flush=True); log.write(s + "\n"); log.flush()

# ------------------------------------------------------------------ vectorised classifiers
def near(x, y, tol=0.15):
    return np.abs(x - y) <= tol * np.maximum(np.maximum(np.abs(x), np.abs(y)), 1.0)

def modes(a, b, f):
    eab, eaf, ebf = near(a, b), near(a, f), near(b, f)
    # strict order codes for the non-tied case
    order = np.where(a < b, np.where(b < f, "M7", np.where(a < f, "M8", "M11")),
                     np.where(a < f, "M9", np.where(b < f, "M10", "M12")))
    m1 = np.select([eab & eaf, eab, eaf, ebf],
                   ["Conserved", np.where(f > a, "M1", "M2"), np.where(b > a, "M3", "M4"), np.where(a > b, "M5", "M6")],
                   default=order)
    hi, lo = np.maximum(a, b), np.minimum(a, b)
    m2 = np.select([(f > hi) & ~near(f, hi), (f < lo) & ~near(f, lo), near(f, hi), near(f, lo)],
                   ["H2P", "L2P", "CHP", "CLP"], default="B2P")
    m3 = np.select([m1 == "Conserved", np.isin(m2, ["H2P", "L2P"]), np.isin(m2, ["CHP", "CLP"])],
                   ["Conserved", "ODO", "DO"], default="PDO")
    return m1, m3

def regclass(a, b):
    asig, bsig, absig = np.abs(a) >= .75, np.abs(b) >= .5, np.abs(a - b) >= .75
    return np.select([asig & bsig & ~absig, asig & ~bsig & absig, asig & bsig & absig & (b * (a - b) > 0),
                      asig & bsig & absig, ~asig & bsig & absig, ~asig & ~bsig & ~absig],
                     REG[:6], default="VII.Ambiguous")

# ------------------------------------------------------------------ inputs
expr = pd.read_csv(EXPR, sep="\t", usecols=["gene_id_africa_hap2", "stage", "TK_parent_A", "NS_parent_B",
                                             "TN_boke", "FL_seedless"])
expr = expr.drop_duplicates(["gene_id_africa_hap2", "stage"])
sf = pd.read_csv(SF, sep="\t", index_col=0).iloc[:, 0]
genes = np.array(sorted(expr.gene_id_africa_hap2.unique()))
gi = {g: i for i, g in enumerate(genes)}
COLS = [f"{g}{i:02d}" for g in ("FL", "NS", "TK", "TN") for i in range(1, 20)]
cidx = {c: k for k, c in enumerate(COLS)}
counts = np.zeros((len(genes), 76))
src = {"FL": "FL_seedless", "NS": "NS_parent_B", "TK": "TK_parent_A", "TN": "TN_boke"}
r = expr.gene_id_africa_hap2.map(gi).to_numpy()
s = expr.stage.map(SIDX).to_numpy()
for g, c in src.items():
    for st in range(19):
        m = s == st
        counts[r[m], cidx[f"{g}{st+1:02d}"]] = expr[c].to_numpy()[m] * sf[f"{g}{st+1:02d}"]
resid = np.abs(counts - np.round(counts)).max()
P("genes", len(genes), "reconstructed legacy counts: max |x - round(x)| =", f"{resid:.3g}")
counts_exact = counts.copy()
counts = np.round(counts).astype(np.int64)

def size_factors(cm):
    arr = cm.astype(float); pos = arr > 0
    pn = pos.sum(1); ls = np.log(np.where(pos, arr, 1.0)).sum(1)
    gm = np.zeros(arr.shape[0]); gm[pn > 0] = np.exp(ls[pn > 0] / pn[pn > 0])
    u = gm > 0
    ratios = np.where(pos[u], arr[u] / gm[u, None], np.nan)
    f = np.nanmedian(ratios, axis=0)
    return f / np.exp(np.mean(np.log(f)))

f0 = size_factors(counts)
P("size factors re-estimated from reconstructed counts vs published: max rel diff",
  f"{np.max(np.abs(f0 / sf[COLS].to_numpy() - 1)):.3g}")

bk = pd.read_csv(BK, sep="\t").rename(columns={"gene_a": "gene_dura"})
sl = pd.read_csv(SL, sep="\t").rename(columns={"gene_a": "gene_africa"})
bridge = bk[["orthogroup", "gene_dura"]].merge(sl[["orthogroup", "gene_africa"]], on="orthogroup")
asg = pd.read_csv(ASSIGN, sep="\t")
asg = asg.sort_values(["stage_index", "sample"])
asg["rep"] = asg.groupby("stage_index").cumcount() + 1
sc = pd.read_csv(UB / "sample_allele_counts_unified.tsv.gz", sep="\t",
                 usecols=["analysis", "sample", "backbone_gene_id", "allele_A_fragments", "allele_B_fragments"])
sc = sc[sc.analysis == "TN"].drop(columns="analysis").rename(columns={"backbone_gene_id": "gene_dura"})
sc = sc.merge(asg[["sample", "stage_index", "rep"]], on="sample")
sc["stage"] = sc.stage_index.map(lambda i: STAGES[i - 1])
sc = sc[sc.gene_dura.isin(bridge.gene_dura)].reset_index(drop=True)
P("TN per-replicate allele rows (bridged genes):", len(sc), "replicates per stage:", asg.groupby("stage_index").size().unique())
REFA = sc.allele_A_fragments.to_numpy(np.int64); REFB = sc.allele_B_fragments.to_numpy(np.int64)

def pooled_B(ra, rb):
    d = sc[["gene_dura", "stage"]].copy()
    q = (ra + rb) >= 10
    d["qa"] = np.where(q, ra, 0); d["qb"] = np.where(q, rb, 0); d["q"] = q.astype(int)
    g = d.groupby(["gene_dura", "stage"], sort=False).agg(qrep=("q", "sum"), ref=("qa", "sum"), alt=("qb", "sum")).reset_index()
    g = g[(g.qrep >= 2) & (g.ref + g.alt >= 30)]
    g["B"] = np.log2((g.ref + .5) / (g.alt + .5))
    return g[["gene_dura", "stage", "B"]]

SF_PUB = sf[COLS].to_numpy()
def expr_long(cm, exact=False):
    # published size factors, rescaled by the change in median-ratio factors after thinning
    f = SF_PUB if exact else SF_PUB * size_factors(cm) / f0
    nm = cm / f
    blocks = []
    for st in range(19):
        blocks.append(pd.DataFrame({"gene_africa": genes, "stage": STAGES[st],
                                    "TK": nm[:, cidx[f"TK{st+1:02d}"]], "NS": nm[:, cidx[f"NS{st+1:02d}"]],
                                    "TN": nm[:, cidx[f"TN{st+1:02d}"]]}))
    e = pd.concat(blocks, ignore_index=True)
    e = e[e[["TK", "NS", "TN"]].max(axis=1) >= 10].copy()
    e["mode"], e["inh"] = modes(e.TK.to_numpy(), e.NS.to_numpy(), e.TN.to_numpy())
    return e

def reg_table(e, Bt):
    m = Bt.merge(bridge, on="gene_dura").merge(e, on=["gene_africa", "stage"])
    m["A"] = np.log2((m.TK + 1) / (m.NS + 1))
    m["cls"] = regclass(m.A.to_numpy(), m.B.to_numpy())
    return m

def cis_medians(m):
    t = m.B.abs() + (m.A - m.B).abs()
    c = np.where(t > 0, m.B.abs() / t.replace(0, np.nan), np.nan)
    b = pd.cut(m.A.abs(), [-np.inf, 1, 2, 3, 4, np.inf], labels=BINS)
    d = pd.DataFrame({"phase": m.stage.map(PH), "bin": b.astype(str), "cis": c}).dropna()
    return d.groupby(["phase", "bin"]).cis.median()

def residuals(m):
    x = m[m.inh.isin(INH)]
    ct = pd.crosstab(x.inh, x.cls).reindex(index=INH, columns=REG, fill_value=0).to_numpy(float)
    ex = np.outer(ct.sum(1), ct.sum(0)) / ct.sum()
    return (ct - ex) / np.sqrt(ex)

def props(labels, levels):
    v = pd.Series(labels).value_counts(normalize=True)
    return np.array([100 * v.get(l, 0.0) for l in levels])

def kappa(a, b, levels):
    po = np.mean(a == b)
    pa = props(a, levels) / 100; pb = props(b, levels) / 100
    pe = float((pa * pb).sum())
    return po, (po - pe) / (1 - pe)

# ------------------------------------------------------------------ full data (validation)
E0 = expr_long(counts_exact, exact=True)
B0 = pooled_B(REFA, REFB)
M0 = reg_table(E0, B0)
P("3g gene-stage observations:", len(E0), " 3h gene-stage observations:", len(M0))
ref3g = pd.read_csv(REF3G, sep="\t"); ref3g = ref3g[ref3g.version == "dedup"]
my3g = E0.groupby(["stage", "mode"]).size().rename("genes").reset_index().merge(ref3g, on=["stage", "mode"], how="outer", suffixes=("", "_ref"))
P("3g validation vs trace recompute (247 cells): cells", len(my3g), "max |count diff|", int((my3g.genes - my3g.genes_ref).abs().max()))
ref3h = pd.read_csv(REF3H, sep="\t"); ref3h = ref3h[ref3h.variant == "eligible_max10"]
my3h = M0.groupby(["stage", "cls"]).size().rename("genes").reset_index().rename(columns={"cls": "regulatory_class"}).merge(
    ref3h, on=["stage", "regulatory_class"], how="outer", suffixes=("", "_ref"))
P("3h validation vs trace recompute (133 cells): cells", len(my3h), "max |count diff|", int((my3h.genes - my3h.genes_ref).abs().max()))
M0["inh"] = M0.inh  # already from E0
R0 = residuals(M0); C0 = cis_medians(M0)
full_prop = {"3g": props(E0["mode"], MODES), "3h": props(M0.cls, REG)}
P("full 3g %:", dict(zip(MODES, full_prop["3g"].round(2))))
P("full 3h %:", dict(zip(REG, full_prop["3h"].round(2))))
P("full 3i medians:", C0.round(3).to_dict())

# ------------------------------------------------------------------ (1) downsampling
ds_rows, ds_prop = [], []
tk_ns = [cidx[f"{g}{i:02d}"] for g in ("TK", "NS") for i in range(1, 20)]
tn = [cidx[f"TN{i:02d}"] for i in range(1, 20)]
key3g = E0.set_index(["gene_africa", "stage"])
key3h = M0.set_index(["gene_dura", "stage"])
for scen in ["parents", "hybrid", "both"]:
    for frac in [0.5, 0.75]:
        for rep in range(NREP):
            cm = counts.copy()
            cols = (tk_ns if scen in ("parents", "both") else []) + (tn if scen in ("hybrid", "both") else [])
            cm[:, cols] = rng.binomial(cm[:, cols], frac)
            if scen in ("hybrid", "both"):
                Bt = pooled_B(rng.binomial(REFA, frac), rng.binomial(REFB, frac))
            else:
                Bt = B0
            E = expr_long(cm); M = reg_table(E, Bt)
            j3g = key3g[["mode", "inh"]].join(E.set_index(["gene_africa", "stage"])[["mode", "inh"]], how="inner", rsuffix="_ds")
            j3h = key3h[["cls"]].join(M.set_index(["gene_dura", "stage"])[["cls"]], how="inner", rsuffix="_ds")
            po_g, k_g = kappa(j3g["mode"].to_numpy(), j3g["mode_ds"].to_numpy(), MODES)
            po_i, k_i = kappa(j3g["inh"].to_numpy(), j3g["inh_ds"].to_numpy(), INH + ["Conserved"])
            po_h, k_h = kappa(j3h["cls"].to_numpy(), j3h["cls_ds"].to_numpy(), REG)
            R = residuals(M); C = cis_medians(M)
            cc = pd.concat([C0.rename("full"), C.rename("ds")], axis=1).dropna()
            pg, ph = props(E["mode"], MODES), props(M.cls, REG)
            ds_rows.append(dict(scenario=scen, fraction=frac, repeat=rep + 1,
                                n_3g=len(E), retained_3g=len(j3g) / len(E0), agree_3g=po_g, kappa_3g=k_g,
                                agree_inheritance=po_i, kappa_inheritance=k_i,
                                n_3h=len(M), retained_3h=len(j3h) / len(M0), agree_3h=po_h, kappa_3h=k_h,
                                max_abs_pp_change_3g=np.abs(pg - full_prop["3g"]).max(),
                                max_abs_pp_change_3h=np.abs(ph - full_prop["3h"]).max(),
                                r_residuals_3j=np.corrcoef(R.ravel(), R0.ravel())[0, 1],
                                sign_agree_3j_cells_absres_ge2=float(np.mean(np.sign(R[np.abs(R0) >= 2]) == np.sign(R0[np.abs(R0) >= 2]))),
                                max_abs_diff_3i_median=float((cc.full - cc.ds).abs().max()),
                                cis_median_decreasing_all_windows=bool(all(np.all(np.diff(C.loc[p].reindex(BINS).to_numpy()) < 0) for p in PHASES))))
            for lv, v in zip(MODES, pg):
                ds_prop.append(dict(panel="3g", scenario=scen, fraction=frac, repeat=rep + 1, category=lv, percent=v))
            for lv, v in zip(REG, ph):
                ds_prop.append(dict(panel="3h", scenario=scen, fraction=frac, repeat=rep + 1, category=lv, percent=v))
            P(scen, frac, rep + 1, {k: (round(v, 4) if isinstance(v, float) else v) for k, v in ds_rows[-1].items() if k not in ("scenario", "fraction", "repeat")})
ds = pd.DataFrame(ds_rows); ds.to_csv(OUT / "fig3gj_downsampling_runs.tsv", sep="\t", index=False, float_format="%.6g")
dp = pd.DataFrame(ds_prop)
for lv, v in zip(MODES, full_prop["3g"]):
    ds_prop.append(dict(panel="3g", scenario="full", fraction=1.0, repeat=0, category=lv, percent=v))
for lv, v in zip(REG, full_prop["3h"]):
    ds_prop.append(dict(panel="3h", scenario="full", fraction=1.0, repeat=0, category=lv, percent=v))
pd.DataFrame(ds_prop).to_csv(OUT / "fig3gj_downsampling_proportions.tsv", sep="\t", index=False, float_format="%.6g")
summ = ds.groupby(["scenario", "fraction"]).agg(["mean", "min", "max"])
summ.to_csv(OUT / "fig3gj_downsampling_summary.tsv", sep="\t", float_format="%.6g")
P("DOWNSAMPLING SUMMARY\n" + ds.groupby(["scenario", "fraction"])[["agree_3g", "kappa_3g", "agree_inheritance", "agree_3h", "kappa_3h",
    "retained_3g", "retained_3h", "max_abs_pp_change_3g", "max_abs_pp_change_3h", "r_residuals_3j", "max_abs_diff_3i_median",
    "cis_median_decreasing_all_windows"]].mean().round(4).to_string())

# ------------------------------------------------------------------ (2) temporal consistency
DEV, POST = STAGES[:13], STAGES[13:]
pairs = [(DEV[i], DEV[i + 1]) for i in range(12)] + [(POST[i], POST[i + 1]) for i in range(5)]
def temporal(df, gcol, lab, levels, tag, nperm=1000):
    rows = []
    w = df.pivot(index=gcol, columns="stage", values=lab)
    for s1, s2 in pairs + [(DEV[i], DEV[i + k]) for k in range(2, 13) for i in range(13 - k)]:
        x = w[[s1, s2]].dropna()
        a, b = x[s1].to_numpy(), x[s2].to_numpy()
        obs = float(np.mean(a == b))
        adj = (s1, s2) in pairs
        if adj:
            null = np.array([np.mean(a == rng.permutation(b)) for _ in range(nperm)])
            pv = (1 + np.sum(null >= obs)) / (nperm + 1)
            nm, nl, nh = null.mean(), *np.percentile(null, [2.5, 97.5])
        else:
            pa = props(a, levels) / 100; pb = props(b, levels) / 100
            nm, nl, nh, pv = float((pa * pb).sum()), np.nan, np.nan, np.nan
        lag = abs(SIDX[s2] - SIDX[s1])
        rows.append(dict(panel=tag, stage_1=s1, stage_2=s2, adjacent=adj, lag=lag, n_genes=len(x), agreement=obs,
                         null_mean=nm, null_low95=nl, null_high95=nh, perm_P=pv, fold_over_null=obs / nm))
    return rows
trows = temporal(E0, "gene_africa", "mode", MODES, "3g_modes") + \
        temporal(E0[E0.inh.isin(INH)], "gene_africa", "inh", INH, "3j_inheritance") + \
        temporal(M0, "gene_dura", "cls", REG, "3h_regulatory")
tt = pd.DataFrame(trows); tt.to_csv(OUT / "fig3gj_temporal_consistency.tsv", sep="\t", index=False, float_format="%.6g")
ta = tt[tt.adjacent]
P("TEMPORAL (adjacent stages)\n" + ta.groupby("panel")[["agreement", "null_mean", "fold_over_null", "perm_P", "n_genes"]].agg(["mean", "min", "max"]).round(4).to_string())
P("TEMPORAL lag profile (developmental series)\n" + tt[tt.stage_1.isin(DEV) & tt.stage_2.isin(DEV)].groupby(["panel", "lag"]).agreement.mean().unstack(0).round(4).to_string())

# ------------------------------------------------------------------ (3) TN replicates
rp = sc.copy(); rp["n"] = rp.allele_A_fragments + rp.allele_B_fragments
rp = rp[rp.n >= 10].copy()
rp["Br"] = np.log2((rp.allele_A_fragments + .5) / (rp.allele_B_fragments + .5))
w = rp.pivot_table(index=["gene_dura", "stage"], columns="rep", values="Br")
cor_rows = []
for st in STAGES:
    x = w.xs(st, level="stage")
    for i, j in [(1, 2), (1, 3), (2, 3)]:
        y = x[[i, j]].dropna()
        nz = y[(y[i].abs() >= 0.5) & (y[j].abs() >= 0.5)]
        cor_rows.append(dict(stage=st, rep_i=i, rep_j=j, n=len(y), pearson=y[i].corr(y[j]),
                             spearman=y[i].corr(y[j], method="spearman"),
                             n_both_abs_ge_0_5=len(nz), direction_agreement=float(np.mean(np.sign(nz[i]) == np.sign(nz[j]))),
                             direction_agreement_all_nonzero=float(np.mean(np.sign(y[i][(y[i] != 0) & (y[j] != 0)]) == np.sign(y[j][(y[i] != 0) & (y[j] != 0)])))))
cr = pd.DataFrame(cor_rows); cr.to_csv(OUT / "fig3gj_replicate_correlation.tsv", sep="\t", index=False, float_format="%.6g")
P("REPLICATES: pearson median/min/max", cr.pearson.median().round(4), cr.pearson.min().round(4), cr.pearson.max().round(4),
  "spearman median", cr.spearman.median().round(4), "direction agreement (|B|>=0.5 both) median/min",
  cr.direction_agreement.median().round(4), cr.direction_agreement.min().round(4))
# eligible gene-stages (published set) with per-replicate B
rep_rows, rep_prop, rep_cis = [], [], []
base = M0.set_index(["gene_dura", "stage"])
lab = {}
for k in (1, 2, 3):
    b = w[k].dropna().rename("B").reset_index()
    b = b.merge(M0[["gene_dura", "stage"]], on=["gene_dura", "stage"])
    M = reg_table(E0, b)
    lab[k] = M.set_index(["gene_dura", "stage"]).cls
    j = base[["cls"]].join(lab[k].rename("cls_r"), how="inner")
    po, kp = kappa(j.cls.to_numpy(), j.cls_r.to_numpy(), REG)
    R = residuals(M); C = cis_medians(M)
    ph = props(M.cls, REG)
    rep_rows.append(dict(comparison=f"replicate {k} vs pooled", n=len(j), agreement=po, kappa=kp,
                         max_abs_pp_change_3h=np.abs(ph - full_prop["3h"]).max(),
                         r_residuals_3j=np.corrcoef(R.ravel(), R0.ravel())[0, 1],
                         cis_median_decreasing_all_windows=bool(all(np.all(np.diff(C.loc[p].reindex(BINS).to_numpy()) < 0) for p in PHASES))))
    for lv, v in zip(REG, ph):
        rep_prop.append(dict(replicate=k, category=lv, percent=v))
    for (p_, b_), v in C.items():
        rep_cis.append(dict(replicate=k, window=p_, abs_A_bin=b_, median_cis=v))
for i, j2 in [(1, 2), (1, 3), (2, 3)]:
    jj = pd.concat([lab[i].rename("a"), lab[j2].rename("b")], axis=1, join="inner")
    po, kp = kappa(jj.a.to_numpy(), jj.b.to_numpy(), REG)
    rep_rows.append(dict(comparison=f"replicate {i} vs replicate {j2}", n=len(jj), agreement=po, kappa=kp))
for (p_, b_), v in C0.items():
    rep_cis.append(dict(replicate="pooled", window=p_, abs_A_bin=b_, median_cis=v))
for lv, v in zip(REG, full_prop["3h"]):
    rep_prop.append(dict(replicate="pooled", category=lv, percent=v))
rr = pd.DataFrame(rep_rows); rr.to_csv(OUT / "fig3gj_replicate_reclassification.tsv", sep="\t", index=False, float_format="%.6g")
pd.DataFrame(rep_prop).to_csv(OUT / "fig3gj_replicate_proportions.tsv", sep="\t", index=False, float_format="%.6g")
pd.DataFrame(rep_cis).to_csv(OUT / "fig3gj_replicate_cis_medians.tsv", sep="\t", index=False, float_format="%.6g")
P("REPLICATE RECLASSIFICATION\n" + rr.round(4).to_string())
# replicate scatter subsample for plotting (stage 95d)
x = w.xs("95d", level="stage")[[1, 2]].dropna()
x.sample(min(4000, len(x)), random_state=1).reset_index().to_csv(OUT / "fig3gj_replicate_scatter_95d.tsv", sep="\t", index=False, float_format="%.5g")
P("DONE")
log.close()
