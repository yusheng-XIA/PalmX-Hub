#!/usr/bin/env python3
"""Supplementary Fig. 11 (with confounding-control panels e, f; fix/enh_dsv_confound): candidate dSVs re-identified within E. guineensis (African35) under three derived-state
criteria and three carrier thresholds; Fig. 5e,f statistics under each definition.

Inputs : ../enh_B_integrate/sd/SF11a..SF11d (sf11_data.py), ../enh_B/res/coloc_syri/coloc_summary.tsv,
         sd/SF11e_correlations.tsv, sd/SF11f_ratios.tsv (sf11ef_data.py).
Panels a-d are unchanged from ../enh_B_integrate/make_sf11.py; the figure is taller (196 mm) and a-d are shifted up.
Outputs: SF11/Supplementary_Fig_11.pdf (vector, Arial, 180 mm wide), SF11/Supplementary_Fig_11.png (600 dpi).
Style: palette scheme A (fix/beautify/common/palA.py); b and c in neutral greys; layouts follow Fig. 5e,f (observed /
background half-violins with boxes and points on a square-root axis; accumulated profile in 10-kb bins with
bootstrap s.e. on a log axis on #F4F4F4). P-value text 6 pt (exponents 4.2 pt); P values in d to two significant digits.
"""
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from scipy import stats

H = Path(__file__).resolve().parent; SD = H.parent / "enh_B_integrate/sd"; SD2 = H / "sd"; OUT = H / "SF11"; OUT.mkdir(exist_ok=True)
A = pd.read_csv(SD / "SF11a_dsv_definitions.tsv", sep="\t")
Bt = pd.read_csv(SD / "SF11b_sample_chrom.tsv", sep="\t")
Cp = pd.read_csv(SD / "SF11c_dSV_pairs.tsv", sep="\t")
Cf = pd.read_csv(SD / "SF11c_profile.tsv", sep="\t")
D = pd.read_csv(SD / "SF11d_definition_summary.tsv", sep="\t")
S = pd.read_csv(H.parent / "enh_B/res/coloc_syri/coloc_summary.tsv", sep="\t").set_index("label")
MAIN = S.loc["EG35_c_Phoenix_and_Oleifera_k1"]

plt.rcParams.update({"font.family": "Arial", "font.size": 6, "axes.linewidth": 0.5, "xtick.major.width": 0.5,
                     "ytick.major.width": 0.5, "xtick.minor.width": 0.4, "ytick.minor.width": 0.4,
                     "xtick.major.size": 2, "ytick.major.size": 2, "xtick.minor.size": 1.2, "ytick.minor.size": 1.2,
                     "axes.labelsize": 6, "xtick.labelsize": 5.5, "ytick.labelsize": 5.5, "pdf.fonttype": 42,
                     "ps.fonttype": 42, "axes.spines.top": False, "axes.spines.right": False, "legend.frameon": False,
                     "legend.fontsize": 5.5, "axes.unicode_minus": True, "mathtext.fontset": "custom",
                     "mathtext.rm": "Arial", "mathtext.it": "Arial:italic", "mathtext.bf": "Arial:bold"})
# Scheme A (fix/beautify/common/palA.py; audit 2026-09-24): category colours for the three derived-state criteria
# in a/d; neutral greys for b and c, whose colours carry a different meaning (observed/background, fit).
C_PHX, C_OLE, C_BOTH = "#3F7FB5", "#2F6E45", "#C8504A"
PT, FIT = "#6F8394", "#333333"            # b: points, regression line
OBS, BGD = "#4D4D4D", "#A6A6A6"           # c: observed, matched background
OBS_F, BGD_F = "#B8B8B8", "#E0E0E0"       # c: half-violin fills
SCH = [("a", "Date palm", C_PHX), ("b", "E. oleifera", C_OLE), ("c", "Date palm and E. oleifera", C_BOTH)]
P_PT = 6                                  # P-value text; mathtext exponent = 0.7 x 6 = 4.2 pt
THR = ["1/35", "≤2/35", "≤3/35"]


def sci(p, digits=1):
    m, e = f"{p:.{digits}e}".split("e")
    return rf"{m} $\times$ 10$^{{{int(e)}}}$"


def italic_species(lab):
    return lab.replace("E. oleifera", r"$\it{E.\,oleifera}$")


H_OLD, H_NEW = 128.0, 196.0            # mm; the old 128-mm layout (a-d) is kept and shifted up by 68 mm
fig = plt.figure(figsize=(180 / 25.4, H_NEW / 25.4))
Y = lambda y: (y * H_OLD + (H_NEW - H_OLD)) / H_NEW
HH = lambda h: h * H_OLD / H_NEW
_add = fig.add_axes
fig.add_axes = lambda r: _add([r[0], Y(r[1]), r[2], HH(r[3])])
L, Rr = 0.065, 0.985
# ---------------- layout (figure fractions)
ax_a = fig.add_axes([L, 0.595, 0.235, 0.35])
ax_b = fig.add_axes([0.385, 0.595, 0.235, 0.35])
ax_c1 = fig.add_axes([0.705, 0.775, Rr - 0.705, 0.185])
ax_c2 = fig.add_axes([0.705, 0.475, Rr - 0.705, 0.235])
ax_d1 = fig.add_axes([0.2, 0.075, 0.19, 0.33])
ax_d2 = fig.add_axes([0.435, 0.075, 0.19, 0.33])

# ---------------- a: numbers of dSVs per definition
a = ax_a; x = np.arange(3); wd = 0.27
for i, (code, lab, col) in enumerate(SCH):
    y = [int(A[(A.Definition.str.startswith(f"({code})")) & (A.Max_carriers_of_35 == k)].dSV.item()) for k in (1, 2, 3)]
    xs = x + (i - 1) * wd
    a.bar(xs, y, wd * 0.92, color=col, lw=0, label=f"({code}) " + lab)
    for xx, yy in zip(xs, y):
        a.text(xx, yy + 25, f"{yy:,}", ha="center", va="bottom", fontsize=5, rotation=90,
               bbox=dict(boxstyle="square,pad=0.05", fc="white", ec="none"))
a.axhline(1480, color="0.35", ls=(0, (3, 2)), lw=0.6, zorder=0)
a.text(-0.46, 1480 + 25, "All38: 1,480", fontsize=5, color="0.3", va="bottom", ha="left")
a.set_xticks(x); a.set_xticklabels(THR); a.set_xlim(-0.5, 2.5); a.set_ylim(0, 2450)
a.set_xlabel(r"Maximum carriers among 35 $\it{E.\,guineensis}$ haplotypes")
a.set_ylabel("Candidate dSVs")
h = [Patch(color=c, lw=0) for _, _, c in SCH]
leg = a.legend(h, [italic_species(f"({c}) {l}") for c, l, _ in SCH], title="Derived state supported by",
               title_fontsize=5.5, loc="upper left", bbox_to_anchor=(0.0, 1.07), handlelength=0.9, handleheight=0.7,
               labelspacing=0.25, borderaxespad=0.1)
leg._legend_box.align = "left"

# ---------------- b: Fig. 5e for definition (c), 1/35
b = ax_b
xv, yv = Bt["dSV count"].values.astype(float), Bt["dSNP count"].values.astype(float)
sl = stats.linregress(xv, yv); n = len(xv)
xx = np.linspace(0, xv.max(), 100); yh = sl.intercept + sl.slope * xx
se = np.sqrt(np.sum((yv - (sl.intercept + sl.slope * xv)) ** 2) / (n - 2)) * \
    np.sqrt(1 / n + (xx - xv.mean()) ** 2 / np.sum((xv - xv.mean()) ** 2))
tq = stats.t.ppf(0.975, n - 2)
b.fill_between(xx, yh - tq * se, yh + tq * se, color=FIT, alpha=0.15, lw=0)
b.plot(xx, yh, color=FIT, lw=0.8, zorder=3)
b.scatter(xv, yv, s=7, color=PT, edgecolors="white", linewidths=0.25, alpha=0.9, zorder=2)
r, p = stats.pearsonr(xv, yv)
assert abs(r - MAIN.e_r) < 1e-5 and n == int(MAIN.e_N)
b.text(0.05, 0.97, f"$\\it{{r}}$ = {r:.2f}\n$\\it{{P}}$ = {sci(MAIN.e_P)}\n$\\it{{n}}$ = {n}", transform=b.transAxes, va="top", ha="left",
       fontsize=P_PT, linespacing=1.35)
b.set_xlabel("dSV count per haplotype–chromosome"); b.set_ylabel("dSNP count per haplotype–chromosome")
b.set_xlim(-0.8, xv.max() + 1); b.set_ylim(-25, yv.max() * 1.05)

# ---------------- c (top): per-focal totals, observed vs matched background (sqrt axis)
c1 = ax_c1
obs, bg = Cp["Observed dSNP count"].values, Cp["Background dSNP count"].values
fw, iv = lambda v: np.sqrt(v), lambda v: np.square(v)
rng = np.random.default_rng(7)
xmax = max(obs.max(), bg.max())
grid = np.linspace(0, np.sqrt(xmax), 300)
for k, (vals, fc, ec, y0) in enumerate([(obs, OBS_F, OBS, 2.0), (bg, BGD_F, BGD, 0.0)]):
    t = np.sqrt(vals)
    kde = stats.gaussian_kde(t, bw_method=0.25)(grid); kde = kde / kde.max() * 0.75
    c1.fill_between(grid, y0 + 0.55, y0 + 0.55 + kde, color=fc, lw=0.5, ec="0.4")
    q1, q2, q3 = np.percentile(t, [25, 50, 75]); iqr = q3 - q1
    lo, hi = max(t.min(), q1 - 1.5 * iqr), min(t.max(), q3 + 1.5 * iqr)
    c1.add_patch(plt.Rectangle((q1, y0 + 0.3), q3 - q1, 0.18, fc="white", ec="0.15", lw=0.45, zorder=3))
    c1.plot([q2, q2], [y0 + 0.3, y0 + 0.48], color="0.15", lw=0.6, zorder=4)
    c1.plot([lo, q1], [y0 + 0.39] * 2, color="0.15", lw=0.45); c1.plot([q3, hi], [y0 + 0.39] * 2, color="0.15", lw=0.45)
    c1.scatter(t, y0 + 0.08 + rng.uniform(-0.07, 0.07, len(t)), s=1.2, color=ec, alpha=0.45, lw=0)
xb = np.sqrt(xmax) * 1.03
c1.plot([xb - 0.12, xb, xb, xb - 0.12], [2.39, 2.39, 0.39, 0.39], color="0.15", lw=0.5)
c1.text(xb + 0.08, 1.39, "***", rotation=90, va="center", ha="left", fontsize=6)
ticks = [0, 5, 10, 20, 40, 80, 140]
c1.set_xticks(np.sqrt(ticks)); c1.set_xticklabels([str(t) for t in ticks])
c1.set_xlim(-0.15, np.sqrt(xmax) * 1.075); c1.set_ylim(-0.1, 3.4); c1.set_yticks([])
for s in ("left", "top", "right"):
    c1.spines[s].set_visible(True)
c1.set_xlabel("Total dSNPs per focal dSV (±1 Mb)", labelpad=1.5)
c1.text(np.sqrt(xmax) * 0.99, 1.55, f"$\\it{{n}}$ = {int(MAIN.f_pairs)} dSV–carrier pairs\nobserved/background = {MAIN.f_obs_mean / MAIN.f_bg_mean:.1f}",
        ha="right", va="center", fontsize=5, linespacing=1.3)

# ---------------- c (bottom): accumulated profile
c2 = ax_c2
c2.set_facecolor("#F4F4F4")
xm = Cf.Distance_bin_mid_kb.values
for col, ycol, secol, lab in [(OBS, "Observed_accumulated_dSNP", "Observed_se", "Observed"),
                              (BGD, "Background_accumulated_dSNP", "Background_se", "Background")]:
    yv2, e2 = Cf[ycol].values, Cf[secol].values
    c2.errorbar(xm, yv2, yerr=[np.minimum(e2, yv2 - 1.05), e2], fmt="none", ecolor=col, alpha=0.45, elinewidth=0.35,
                capsize=0, zorder=2)
    c2.plot(xm, yv2, "-o", color=col, lw=0.6, ms=0.9, mfc=col, mec=col, mew=0, zorder=3, label=lab)
c2.axvline(0, color="0.35", ls=(0, (2.5, 2)), lw=0.5, zorder=1)
c2.set_yscale("log"); c2.set_ylim(3, 1000); c2.set_xlim(-1000, 1000)
c2.set_yticks([5, 10, 25, 50, 100, 250, 500]); c2.set_yticklabels(["5", "10", "25", "50", "100", "250", "500"])
c2.yaxis.set_minor_locator(matplotlib.ticker.NullLocator())
c2.set_xticks([-1000, -500, 0, 500, 1000]); c2.set_xticklabels(["−1000", "−500", "0", "500", "1000"])
for s in ("top", "right"):
    c2.spines[s].set_visible(True)
c2.set_xlabel("Distance to dSV (kb)"); c2.set_ylabel("Accumulated dSNP count")
c2.legend(loc="upper right", handlelength=1.4, borderaxespad=0.3, markerscale=1.5)
c2.text(0.02, 0.97, f"$\\it{{P}}$ = {sci(MAIN.f_wilcoxon_P)}", transform=c2.transAxes, fontsize=P_PT, va="top", ha="left")

# ---------------- d: nine definitions
Dn = D[D.Code != "All38"].copy(); ref = D[D.Code == "All38"].iloc[0]
order = [f"{c}{k}" for k in (1, 2, 3) for c in "abc"]
Dn = Dn.set_index("Code").loc[order].reset_index()
ys = []
y = 0.0
for k in (1, 2, 3):
    for c in "abc":
        ys.append(y); y -= 1
    y -= 0.6
ys = np.array(ys)
cols = {c: col for c, _, col in SCH}
labs = [f"({r.Code[0]}) {THR[int(r.Code[1]) - 1]}  n = {r.dSV_sites:,}" for r in Dn.itertuples()]
for axd, col, refv, xl, lim in [(ax_d1, "Pearson_r", ref.Pearson_r, r"Pearson $\it{r}$, dSV vs dSNP (Fig. 5e)", (0.7, 1.0)),
                                (ax_d2, "Observed_over_background", ref.Observed_over_background,
                                 "Observed/background dSNPs, ±1 Mb (Fig. 5f)", (0, 7))]:
    for yy, rr in zip(ys, Dn.itertuples()):
        axd.plot([lim[0], getattr(rr, col)], [yy, yy], color="0.85", lw=0.5, zorder=1)
        axd.scatter(getattr(rr, col), yy, s=14, color=cols[rr.Code[0]], edgecolors="white", linewidths=0.3, zorder=3)
    axd.axvline(refv, color="0.35", ls=(0, (3, 2)), lw=0.6, zorder=2)
    axd.set_xlim(*lim); axd.set_ylim(ys.min() - 0.7, ys.max() + 0.9)
    axd.set_xlabel(xl)
    axd.text(refv, ys.max() + 0.75, f"All38: {refv:.2f}" if col == "Pearson_r" else f"All38: {refv:.1f}",
             fontsize=5, color="0.3", ha="center", va="bottom")
ax_d1.set_yticks(ys); ax_d1.set_yticklabels(labs)
for tl, rr in zip(ax_d1.get_yticklabels(), Dn.itertuples()):
    tl.set_color(cols[rr.Code[0]])
ax_d2.set_yticks(ys); ax_d2.tick_params(labelleft=False)
ax_d2.set_xticks([0, 1, 2, 3, 4, 5, 6, 7])
ax_d2.axvline(1, color="0.6", lw=0.4, zorder=1)
for yy, rr in zip(ys, Dn.itertuples()):
    ax_d1.text(0.992, yy, f"{rr.Pearson_r:.2f}", fontsize=5, va="center", ha="right", color="0.25")
    ax_d2.text(7.25, yy, f"$\\it{{P}}$ = {sci(float(rr.Wilcoxon_P), 1)}", fontsize=P_PT, va="center", ha="left", color="0.25",
               clip_on=False)
h = [Line2D([], [], marker="o", ls="", color=c, mec="white", mew=0.3, ms=4) for _, _, c in SCH]
fig.legend(h, [italic_species(f"({c}) {l}") for c, l, _ in SCH], title="Derived state supported by",
           title_fontsize=5.5, loc="upper left", bbox_to_anchor=(0.765, Y(0.33)), handletextpad=0.3, labelspacing=0.35)
fig.text(0.772, Y(0.19), f"All nine definitions:\nPearson $\\it{{P}}$ ≤ {sci(D[D.Code != 'All38'].Pearson_P.astype(float).max())}\n"
         f"Wilcoxon $\\it{{P}}$ ≤ {sci(D[D.Code != 'All38'].Wilcoxon_P.astype(float).max())}\n$\\it{{n}}$ = 560 haplotype–chromosome\ncombinations",
         fontsize=P_PT, va="top", ha="left", linespacing=1.35)

fig.add_axes = _add

# ---------------- e, f: confounding controls (fix/enh_dsv_confound)
E = pd.read_csv(SD2 / "SF11e_correlations.tsv", sep="\t")
F = pd.read_csv(SD2 / "SF11f_ratios.tsv", sep="\t")
SETS = [("All38 (1,480 dSVs; 38 haplotypes)", "All38; $\\it{n}$ = 608"),
        ("Definition (c), 1/35 (550 dSVs; 35 haplotypes)", "Definition (c), 1/35; $\\it{n}$ = 560")]
CLS = [("dSNP (missense, stop, start, splice or conserved-proxy)", "dSNP", dict(marker="o", color="#333333", mfc="#333333")),
       ("Derived non-coding SNV (intron, UTR, intergenic or +-2 kb; outside conserved proxy)", "Derived non-coding SNV",
        dict(marker="s", color="#6F8394", mfc="#6F8394")),
       ("Derived synonymous SNV (all)", "Derived synonymous SNV", dict(marker="^", color="#8C8C8C", mfc="white")),
       ("Rare non-coding SNV, polarity not required", "Rare non-coding SNV, polarity not required",
        dict(marker="D", color="#B0B0B0", mfc="#B0B0B0"))]
OFF = [0.24, 0.08, -0.08, -0.24]
EROWS = [("Counts", "Counts (Fig. 5e)"), ("Counts, haplotype mean removed", "Haplotype mean\nremoved"),
         ("Counts, haplotype and chromosome effects removed (two-way fixed effects)", "Haplotype and\nchromosome FE"),
         ("Per Mb of chromosome", "Per Mb of\nchromosome"), ("Per Mb of CDS or conserved sequence", "Per Mb of CDS or\nconserved sequence"),
         ("Two-way fixed effects plus derived non-coding and synonymous counts", "FE + neutral\ncounts (dSNP)")]
FROWS = [("B0", "Random position on\nthe chromosome (Fig. 5f)"), ("B1", "Matched CDS +\nconserved content"),
         ("B2", "Matched content,\ncentre in CDS or\nconserved sequence")]
yb, hb = 10 / H_NEW, 36 / H_NEW
axe = [fig.add_axes([0.155, yb, 0.14, hb]), fig.add_axes([0.305, yb, 0.14, hb])]
axf = [fig.add_axes([0.64, yb, 0.16, hb]), fig.add_axes([0.81, yb, 0.16, hb])]
for j, (ax, (sn, st)) in enumerate(zip(axe, SETS)):
    for i, (mk, _) in enumerate(EROWS):
        y0 = -i
        ax.axhspan(y0 - 0.45, y0 + 0.45, color="#F4F4F4" if i % 2 == 0 else "white", lw=0, zorder=0)
        for (cn, _, st_), off in zip(CLS, OFF):
            v = E[(E.Set == sn) & (E["SNV class"] == cn) & (E.Adjustment == mk)].Pearson_r
            if mk.startswith("Two-way fixed effects plus") and cn != CLS[0][0]:
                continue
            assert len(v) == 1, (sn, cn, mk)
            ax.plot([v.item()], [y0 + off], ls="", ms=3.2, mew=0.5, zorder=3, **st_)
    ax.set_xlim(0, 1); ax.set_ylim(-len(EROWS) + 0.5, 0.5); ax.set_xticks([0, 0.5, 1]); ax.set_xticklabels(["0", "0.5", "1"])
    ax.set_xticks([0.25, 0.75], minor=True)
    ax.set_yticks(-np.arange(len(EROWS)))
    ax.set_yticklabels([l for _, l in EROWS] if j == 0 else [])
    ax.tick_params(axis="y", length=0)
    ax.set_title(st, fontsize=5.5, pad=2)
    ax.set_xlabel("Pearson $\\it{r}$ with dSV count", labelpad=1.5)
for j, (ax, (sn, st)) in enumerate(zip(axf, SETS)):
    for i, (bg, _) in enumerate(FROWS):
        y0 = -i
        ax.axhspan(y0 - 0.45, y0 + 0.45, color="#F4F4F4" if i % 2 == 0 else "white", lw=0, zorder=0)
        for (cn, _, st_), off in zip(CLS, OFF):
            q = F[(F.Set == sn) & (F["SNV class"] == cn) & (F.Background.str.startswith(bg + ":"))]
            assert len(q) == 1, (sn, cn, bg)
            q = q.iloc[0]
            ax.plot([q["Ratio 95% CI low"], q["Ratio 95% CI high"]], [y0 + off] * 2, color=st_["color"], lw=0.6, zorder=2)
            ax.plot([q.Ratio], [y0 + off], ls="", ms=3.2, mew=0.5, zorder=3, **st_)
        if bg == "B2":
            rr = F[(F.Set == sn) & (F["SNV class"] == "dSNP ratio / " + CLS[1][0] + " ratio") & (F.Background.str.startswith("B2:"))].iloc[0]
            ax.text(7.7, y0, f"dSNP/non-coding\n{rr.Ratio:.2f} ($\\it{{P}}$ = {float(rr.P):.2f})", fontsize=5, ha="right",
                    va="center", color="0.25", linespacing=1.15)
    ax.axvline(1, color="0.6", lw=0.4, zorder=1)
    ax.set_xscale("log"); ax.set_xlim(0.9, 8); ax.set_ylim(-len(FROWS) + 0.5, 0.5)
    ax.set_xticks([1, 2, 4, 8]); ax.set_xticklabels(["1", "2", "4", "8"]); ax.xaxis.set_minor_locator(matplotlib.ticker.NullLocator())
    ax.set_yticks(-np.arange(len(FROWS)))
    ax.set_yticklabels([l for _, l in FROWS] if j == 0 else [])
    ax.tick_params(axis="y", length=0)
    ax.set_title(st, fontsize=5.5, pad=2)
    ax.set_xlabel("Observed/background, ±1 Mb", labelpad=1.5)
h = [Line2D([], [], ls="", ms=3.5, mew=0.5, **st_) for _, _, st_ in CLS]
fig.legend(h, [l for _, l, _ in CLS], loc="lower center", bbox_to_anchor=(0.52, 55.5 / H_NEW), ncol=4, handletextpad=0.2,
           columnspacing=1.2, title="Rare derived variants carried by the dSV haplotypes (same frequency and date palm polarization rules unless stated)",
           title_fontsize=5.5)

# ---------------- panel letters
for lab, (xp, yp) in {"a": (0.008, Y(0.985)), "b": (0.328, Y(0.985)), "c": (0.648, Y(0.985)), "d": (0.008, Y(0.455)),
                      "e": (0.008, 52.5 / H_NEW), "f": (0.49, 52.5 / H_NEW)}.items():
    fig.text(xp, yp, lab, fontsize=8, fontweight="bold", va="top", ha="left")
fig.savefig(OUT / "Supplementary_Fig_11.pdf")
fig.savefig(OUT / "Supplementary_Fig_11.png", dpi=600)
print("saved", OUT)
