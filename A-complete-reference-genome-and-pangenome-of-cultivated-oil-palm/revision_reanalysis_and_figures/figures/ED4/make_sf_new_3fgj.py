#!/usr/bin/env python3
"""Supplementary_Fig_new_3fgj: robustness of Fig. 3f (Ka/Ks equivalence among ASE classes) and Fig. 3g-j
(TN parent-hybrid expression modes and cis/trans proxy classes).

Inputs : out/ (cluster work/enh_3fgj/out; s01_kaks_equiv.py, s02_modes_robust.py)
Outputs: Supplementary_Fig_new_3fgj.pdf (vector, Arial, 180 mm wide), .png (600 dpi),
         source_data/Supplementary_Fig_new_3fgj_<panel>.tsv
Colours: palette scheme A (fix/beautify/common/palA.py).
"""
import sys
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

H = Path(__file__).resolve().parent; O = H / "out"; SD = H / "source_data"; SD.mkdir(exist_ok=True)
sys.path.insert(0, str(H.parent / "beautify/common"))
from palA import FL, TN, ASE, EXPAND, CONTRACT, tint  # noqa: E402

plt.rcParams.update({"font.family": "Arial", "font.size": 6, "axes.linewidth": 0.5, "xtick.major.width": 0.5,
                     "ytick.major.width": 0.5, "xtick.major.size": 2, "ytick.major.size": 2, "axes.labelsize": 6,
                     "xtick.labelsize": 5.5, "ytick.labelsize": 5.5, "pdf.fonttype": 42, "ps.fonttype": 42,
                     "axes.spines.top": False, "axes.spines.right": False, "legend.frameon": False,
                     "legend.fontsize": 5.5, "axes.unicode_minus": True, "mathtext.fontset": "custom",
                     "mathtext.rm": "Arial", "mathtext.it": "Arial:italic", "mathtext.bf": "Arial:bold"})
MAT = {"FL": FL, "TN": TN}
C3H, C3G, C3J = CONTRACT, EXPAND, ASE["Sub"]          # cis/trans classes, 12 modes, PDO/DO/ODO
CL = ["NoDiff", "HapDom", "Sub", "NoASE"]
REG = ["I.Cis_only", "II.Trans_only", "III.Cis_trans_enhancing", "IV.Cis_trans_compensating", "V.Compensatory",
       "VI.Conserved", "VII.Ambiguous"]
REGL = ["I", "II", "III", "IV", "V", "VI", "VII"]
PHASES = ["Days 0–65", "Days 80–140", "Days 155–185", "Hours 12–72"]
PHL = ["0–65 d", "80–140 d", "155–185 d", "12–72 h"]
BINS = ["0–1", "1–2", "2–3", "3–4", "4+"]
KAKS, KS = r"$K_{\rm a}/K_{\rm s}$", r"$K_{\rm s}$"


def sci(p):
    if p >= 0.001:
        return f"{p:.3f}" if p < 0.01 else f"{p:.2f}"
    m, e = f"{p:.0e}".split("e")
    return rf"{m} $\times$ 10$^{{{int(e)}}}$"


LX = {"a": 0.01, "b": 0.50, "c": 0.01, "d": 0.355, "e": 0.665, "f": 0.01, "g": 0.30, "h": 0.60}


def letter(ax, s, x=None, y=None):
    top = ax.get_position().y1
    fig.text(LX[s], top + (0.035 if s in "ab" else 0.02), s, fontsize=8, fontweight="bold", va="bottom", ha="left")


fig = plt.figure(figsize=(180 / 25.4, 170 / 25.4))   # 2026-09-24: Extended Data limit 180 x 170 mm (was 172)

# ================= a, b: pairwise Hodges-Lehmann shifts with 90%/95% CIs and equivalence margins
eq = pd.read_csv(O / "fig3f_pairwise_equivalence.tsv", sep="\t")
bg = pd.read_csv(O / "fig3f_genome_background.tsv", sep="\t")
sd_ab = []
for k, (met, lab, lt, scale, unit) in enumerate([("KaKs", KAKS, "a", 1, ""), ("Ks", KS, "b", 1e3, r" ($\times$10$^{-3}$)")]):
    ax = fig.add_axes([0.175 + k * 0.49, 0.705, 0.235, 0.235])
    y = 0; yt, yl = [], []
    for an in ["FL", "TN"]:
        s = eq[(eq.hybrid == an) & (eq.metric == met)]
        m = bg[(bg.hybrid == an) & (bg.metric == met)].margin_10pct_of_genome_median.item() * scale
        y0 = y
        for r in s.itertuples():
            c = MAT[an]
            ax.plot([r.HL_CI95_low * scale, r.HL_CI95_high * scale], [y, y], color=c, lw=0.6, solid_capstyle="butt")
            ax.plot([r.HL_CI90_low * scale, r.HL_CI90_high * scale], [y, y], color=c, lw=2.0, solid_capstyle="butt")
            ax.plot(r.HL_shift * scale, y, "o", ms=2.6, mfc="white", mec=c, mew=0.6, zorder=3)
            ax.text(1.02, y, sci(r.TOST_P), transform=ax.get_yaxis_transform(), fontsize=5, va="center", ha="left",
                    color="k" if r.TOST_P < 0.05 else "0.55")
            yt.append(y); yl.append(f"{r.class_1} – {r.class_2}")
            sd_ab.append(dict(panel=lt, hybrid=an, metric=met, class_1=r.class_1, class_2=r.class_2, n_1=r.n_1, n_2=r.n_2,
                              HL_shift=r.HL_shift, HL_CI90_low=r.HL_CI90_low, HL_CI90_high=r.HL_CI90_high,
                              HL_CI95_low=r.HL_CI95_low, HL_CI95_high=r.HL_CI95_high,
                              median_difference=r.median_difference,
                              median_difference_boot95_low=r.median_difference_boot95_low,
                              median_difference_boot95_high=r.median_difference_boot95_high,
                              equivalence_margin=r.margin, TOST_P=r.TOST_P,
                              smallest_equivalence_margin=r.smallest_equivalence_margin,
                              MWU_P_two_sided=r.MWU_P_two_sided, cliffs_delta=r.cliffs_delta,
                              cliffs_delta_CI90_low=r.cliffs_delta_CI90_low, cliffs_delta_CI90_high=r.cliffs_delta_CI90_high))
            y -= 1
        ax.fill_between([-m, m], y0 + 0.5, y + 0.5, color=tint(MAT[an], 0.18), lw=0, zorder=0)
        ax.text(-0.02, (y0 + y + 1) / 2, an, transform=ax.get_yaxis_transform(), fontsize=6, fontweight="bold",
                color=MAT[an], ha="right", va="center")
        y -= 0.6
    ax.axvline(0, color="0.4", lw=0.4, ls=(0, (2, 2)), zorder=1)
    ax.set_yticks(yt); ax.set_yticklabels(yl, fontsize=5)
    ax.tick_params(axis="y", length=0, pad=18)
    ax.spines["left"].set_visible(False)
    ax.set_ylim(y + 0.3, 0.7)
    ax.set_xlabel(f"Hodges–Lehmann difference in {lab}{unit}")
    ax.text(1.02, 1.0, r"TOST $P$", transform=ax.transAxes, fontsize=5, ha="left", va="bottom")
    letter(ax, lt, x=-0.52)
    if k == 0:
        h = [Line2D([], [], color="0.3", lw=2.0), Line2D([], [], color="0.3", lw=0.6),
             plt.Rectangle((0, 0), 1, 1, color="0.88", lw=0)]
        ax.legend(h, ["90% CI", "95% CI", "Margin (±10% genome median)"], loc="lower left", bbox_to_anchor=(-0.62, 1.0),
                  ncol=3, fontsize=5, handlelength=1.6, columnspacing=1.0, borderaxespad=0.2)
pd.DataFrame(sd_ab).to_csv(SD / "Supplementary_Fig_new_3fgj_ab.tsv", sep="\t", index=False, float_format="%.6g")
pd.DataFrame(bg).to_csv(SD / "Supplementary_Fig_new_3fgj_ab_genome_background.tsv", sep="\t", index=False, float_format="%.6g")

# ================= c: gene-level agreement after downsampling
ds = pd.read_csv(O / "fig3gj_downsampling_runs.tsv", sep="\t")
ax = fig.add_axes([0.085, 0.385, 0.25, 0.225]); letter(ax, "c", x=-0.30)
SC = [("parents", "Parents"), ("hybrid", "Hybrid"), ("both", "Both")]
x = 0; xt, xl = [], []
sd_c = []
for sc, sl in SC:
    for fr in [0.75, 0.5]:
        s = ds[(ds.scenario == sc) & (ds.fraction == fr)]
        for j, (col, c, mk) in enumerate([("agree_3g", C3G, "s"), ("agree_inheritance", C3J, "^"), ("agree_3h", C3H, "o")]):
            v = 100 * s[col]
            ax.plot([x + (j - 1) * 0.22] * 2, [v.min(), v.max()], color=c, lw=0.6)
            ax.plot(x + (j - 1) * 0.22, v.mean(), mk, ms=3, mfc=c, mec="none")
        for r in s.itertuples():
            sd_c.append(dict(panel="c", scenario=sc, fraction=fr, repeat=r.repeat,
                             agreement_12_modes=r.agree_3g, kappa_12_modes=r.kappa_3g,
                             agreement_PDO_DO_ODO=r.agree_inheritance, kappa_PDO_DO_ODO=r.kappa_inheritance,
                             agreement_cis_trans=r.agree_3h, kappa_cis_trans=r.kappa_3h,
                             retained_12_modes=r.retained_3g, retained_cis_trans=r.retained_3h,
                             max_abs_pp_change_12_modes=r.max_abs_pp_change_3g,
                             max_abs_pp_change_cis_trans=r.max_abs_pp_change_3h,
                             r_residuals_3j=r.r_residuals_3j, max_abs_diff_cis_contribution_median=r.max_abs_diff_3i_median,
                             cis_contribution_decreasing_with_absA_all_windows=r.cis_median_decreasing_all_windows))
        xt.append(x); xl.append(f"{int(fr * 100)}%"); x += 1
    x += 0.4
for i, (sc, sl) in enumerate(SC):
    ax.text((xt[2 * i] + xt[2 * i + 1]) / 2, -0.2, sl, transform=ax.get_xaxis_transform(), ha="center", va="top", fontsize=5.5)
ax.set_xticks(xt); ax.set_xticklabels(xl)
ax.set_ylim(75, 100); ax.set_ylabel("Gene–stage assignments unchanged (%)")
ax.set_xlim(-0.6, x - 0.8)
ax.legend([Line2D([], [], marker=m, ls="", mfc=c, mec="none", ms=3) for m, c in [("o", C3H), ("^", C3J), ("s", C3G)]],
          ["Cis/trans classes (Fig. 3h)", "PDO/DO/ODO (Fig. 3j)", "12 modes (Fig. 3g)"], loc="lower left",
          fontsize=5, handletextpad=0.2, borderaxespad=0.1)
pd.DataFrame(sd_c).to_csv(SD / "Supplementary_Fig_new_3fgj_c.tsv", sep="\t", index=False, float_format="%.6g")

# ================= d: cis/trans class proportions, full data vs 50% downsampling and single replicates
dp = pd.read_csv(O / "fig3gj_downsampling_proportions.tsv", sep="\t")
rp = pd.read_csv(O / "fig3gj_replicate_proportions.tsv", sep="\t")
ax = fig.add_axes([0.42, 0.385, 0.25, 0.225]); letter(ax, "d", x=-0.22)
full = dp[(dp.panel == "3h") & (dp.scenario == "full")].set_index("category").percent.reindex(REG)
ax.bar(range(7), full.values, 0.72, color=tint(C3H, 0.35), lw=0, label="All data")
sd_d = [dict(panel="d", source="all data", category=c, percent=v) for c, v in full.items()]
offs = {"parents": -0.24, "hybrid": -0.08, "both": 0.08}
mk = {"parents": "o", "hybrid": "s", "both": "^"}
for sc in ["parents", "hybrid", "both"]:
    s = dp[(dp.panel == "3h") & (dp.scenario == sc) & (dp.fraction == 0.5)].groupby("category").percent.agg(["mean", "min", "max"]).reindex(REG)
    ax.plot(np.arange(7) + offs[sc], s["mean"], mk[sc], ms=2.4, mfc="white", mec="0.25", mew=0.5, ls="")
    for c, r in s.iterrows():
        sd_d.append(dict(panel="d", source=f"{sc} reads thinned to 50% (mean of 10)", category=c, percent=r["mean"],
                         percent_min=r["min"], percent_max=r["max"]))
rr = rp[rp.replicate.astype(str) != "pooled"]
for k in ["1", "2", "3"]:
    s = rr[rr.replicate.astype(str) == k].set_index("category").percent.reindex(REG)
    ax.plot(np.arange(7) + 0.24, s.values, "D", ms=2.0, mfc=C3H, mec="none", alpha=0.8, ls="")
    for c, v in s.items():
        sd_d.append(dict(panel="d", source=f"TN replicate {k} alone", category=c, percent=v))
ax.set_xticks(range(7)); ax.set_xticklabels(REGL)
ax.set_xlabel("Cis/trans class"); ax.set_ylabel("Gene–stage observations (%)")
ax.set_ylim(0, 45)
ax.legend([plt.Rectangle((0, 0), 1, 1, color=tint(C3H, 0.35), lw=0)] +
          [Line2D([], [], marker=mk[s], ls="", mfc="white", mec="0.25", mew=0.5, ms=2.6) for s in offs] +
          [Line2D([], [], marker="D", ls="", mfc=C3H, mec="none", ms=2.2)],
          ["All data", "Parents 50%", "Hybrid 50%", "Both 50%", "Single TN replicate"], loc="upper left", fontsize=5,
          handletextpad=0.2, borderaxespad=0.1, labelspacing=0.3)
pd.DataFrame(sd_d).to_csv(SD / "Supplementary_Fig_new_3fgj_d.tsv", sep="\t", index=False, float_format="%.6g")

# ================= e: adjacent-stage agreement vs permutation null
tt = pd.read_csv(O / "fig3gj_temporal_consistency.tsv", sep="\t")
ta = tt[tt.adjacent].copy()
ax = fig.add_axes([0.745, 0.385, 0.235, 0.225]); letter(ax, "e", x=-0.25)
pairs = ta[ta.panel == "3h_regulatory"][["stage_1", "stage_2"]].apply(lambda r: f"{r.stage_1}–{r.stage_2}", axis=1).tolist()
xx = np.arange(len(pairs)); xx = np.where(xx >= 12, xx + 0.6, xx)
for pan, c, mkr, lab in [("3h_regulatory", C3H, "o", "Cis/trans classes"), ("3j_inheritance", C3J, "^", "PDO/DO/ODO"),
                         ("3g_modes", C3G, "s", "12 modes")]:
    s = ta[ta.panel == pan]
    ax.fill_between(xx, 100 * s.null_low95, 100 * s.null_high95, color=tint(c, 0.3), lw=0, step=None)
    ax.plot(xx, 100 * s.null_mean, color=c, lw=0.5, ls=(0, (2, 1.5)))
    ax.plot(xx, 100 * s.agreement, marker=mkr, ms=2.4, color=c, lw=0.6, label=lab)
ax.set_xticks(xx); ax.set_xticklabels(pairs, rotation=90, fontsize=5)
ax.set_ylabel("Same assignment in adjacent stages (%)"); ax.set_ylim(0, 75)
ax.axvline(11.8, color="0.7", lw=0.4)
hh, ll = ax.get_legend_handles_labels()
hh += [Line2D([], [], color="0.4", lw=0.5, ls=(0, (2, 1.5))), plt.Rectangle((0, 0), 1, 1, color="0.85", lw=0)]
ll += ["Permutation mean", "Permutation 95% range"]
ax.legend(hh, ll, loc="upper left", fontsize=5, handletextpad=0.3, borderaxespad=0.1, labelspacing=0.25, ncol=2,
          columnspacing=0.6, handlelength=1.4)
ax.set_ylim(0, 85)
ta[["panel", "stage_1", "stage_2", "n_genes", "agreement", "null_mean", "null_low95", "null_high95", "perm_P", "fold_over_null"]] \
    .assign(panel=lambda d: "e:" + d.panel).to_csv(SD / "Supplementary_Fig_new_3fgj_e.tsv", sep="\t", index=False, float_format="%.6g")
tt.to_csv(SD / "Supplementary_Fig_new_3fgj_e_all_stage_pairs.tsv", sep="\t", index=False, float_format="%.6g")

# ================= f: replicate-level allelic ratios
cr = pd.read_csv(O / "fig3gj_replicate_correlation.tsv", sep="\t")
sc95 = pd.read_csv(O / "fig3gj_replicate_scatter_95d.tsv", sep="\t")
ax = fig.add_axes([0.085, 0.07, 0.20, 0.215]); letter(ax, "f", x=-0.36)
ax.scatter(sc95["1"], sc95["2"], s=1.2, color=TN, lw=0, alpha=0.45, rasterized=True)
lim = 8.5
ax.plot([-lim, lim], [-lim, lim], color="0.4", lw=0.4, ls=(0, (2, 2)))
ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim); ax.set_aspect("equal")
ax.set_xticks([-8, -4, 0, 4, 8]); ax.set_yticks([-8, -4, 0, 4, 8])
ax.set_xlabel(r"TN replicate 1, log$_2$(A/B)"); ax.set_ylabel(r"TN replicate 2, log$_2$(A/B)")
r95 = cr[(cr.stage == "95d") & (cr.rep_i == 1) & (cr.rep_j == 2)].iloc[0]
ax.text(0.03, 0.97, f"95 d, n = {int(r95.n):,}\n" + rf"$r$ = {r95.pearson:.3f}" +
        f"\nAll stages and pairs:\n" + rf"$r$ = {cr.pearson.min():.2f}–{cr.pearson.max():.2f}" +
        f"\nSame sign (|log$_2$(A/B)| ≥ 0.5\nin both): {100 * cr.direction_agreement.min():.1f}–{100 * cr.direction_agreement.max():.1f}%",
        transform=ax.transAxes, va="top", ha="left", fontsize=5, linespacing=1.15)
sc95.rename(columns={"1": "replicate_1_log2AB", "2": "replicate_2_log2AB"}).assign(panel="f", stage="95d") \
    .to_csv(SD / "Supplementary_Fig_new_3fgj_f_scatter_95d.tsv", sep="\t", index=False, float_format="%.5g")
cr.assign(panel="f").to_csv(SD / "Supplementary_Fig_new_3fgj_f_replicate_correlation.tsv", sep="\t", index=False, float_format="%.6g")

# ================= g: replicate reclassification (agreement, kappa) and direction agreement
rc = pd.read_csv(O / "fig3gj_replicate_reclassification.tsv", sep="\t")
ax = fig.add_axes([0.37, 0.07, 0.17, 0.215]); letter(ax, "g", x=-0.42)
labs = ["Rep. 1 vs pooled", "Rep. 2 vs pooled", "Rep. 3 vs pooled", "Rep. 1 vs 2", "Rep. 1 vs 3", "Rep. 2 vs 3"]
yy = np.arange(6)[::-1]
ax.barh(yy, 100 * rc.agreement, 0.62, color=[C3H] * 3 + [tint(C3H, 0.5)] * 3, lw=0)
for yv, a_, k_ in zip(yy, rc.agreement, rc.kappa):
    ax.text(100 * a_ - 1.5, yv, f"{100 * a_:.1f}% (κ = {k_:.2f})", ha="right", va="center", fontsize=5, color="white")
ax.set_yticks(yy); ax.set_yticklabels(labs, fontsize=5)
ax.set_xlim(0, 100); ax.set_xlabel("Cis/trans assignments unchanged (%)")
ax.tick_params(axis="y", length=0)
rc.assign(panel="g").to_csv(SD / "Supplementary_Fig_new_3fgj_g.tsv", sep="\t", index=False, float_format="%.6g")

# ================= h: cis contribution vs |A| by window: pooled, single replicates, 50% downsampling range
cm = pd.read_csv(O / "fig3gj_replicate_cis_medians.tsv", sep="\t")
sd_h = []
for w, (ph, pl) in enumerate(zip(PHASES, PHL)):
    ax = fig.add_axes([0.655 + w * 0.083, 0.07, 0.07, 0.215])
    if w == 0:
        letter(ax, "h", x=-0.95)
        ax.set_ylabel("Median cis contribution")
    else:
        ax.set_yticklabels([])
    for k, col, lw in [("1", tint(TN, 0.55), 0.6), ("2", tint(TN, 0.55), 0.6), ("3", tint(TN, 0.55), 0.6), ("pooled", "k", 0.9)]:
        s = cm[(cm.replicate.astype(str) == k) & (cm.window == ph)].set_index("abs_A_bin").median_cis.reindex(BINS)
        ax.plot(range(5), s.values, color=col, lw=lw, marker="o" if k == "pooled" else None, ms=1.8,
                zorder=3 if k == "pooled" else 2)
        for b, v in s.items():
            sd_h.append(dict(panel="h", window=pl, abs_A_bin=b, source="pooled replicates" if k == "pooled" else f"TN replicate {k}",
                             median_cis_contribution=v))
    ax.set_ylim(0.28, 0.52); ax.set_xticks(range(5)); ax.set_xticklabels(BINS, rotation=90, fontsize=5)
    ax.set_title(pl, fontsize=5.5, pad=2)
    if w == 1:
        ax.set_xlabel("|A| bin (parental log$_2$(TK/NS))", x=1.1)
fig.axes[-4].legend([Line2D([], [], color="k", marker="o", ms=1.8, lw=0.9), Line2D([], [], color=tint(TN, 0.55), lw=0.6)],
                    ["Pooled", "Replicate"], loc="lower left", fontsize=5, handlelength=1.2, borderaxespad=0.1)
pd.DataFrame(sd_h).to_csv(SD / "Supplementary_Fig_new_3fgj_h.tsv", sep="\t", index=False, float_format="%.6g")

fig.savefig(H / "Supplementary_Fig_new_3fgj.pdf")
fig.savefig(H / "Supplementary_Fig_new_3fgj.png", dpi=600)
print("saved")
