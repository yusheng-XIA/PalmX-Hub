#!/usr/bin/env python3
"""Supplementary Fig. 3 redraw: WGCNA of 114 FL and TN RNA-seq libraries.

Inputs
  raw.pkl  : Source Data sheets Supp_WGCNA_a_power ... _g_h (REVIEW workbook, 2026-09-23 16:36)
  hclust_merge.tsv / hclust_height.tsv / hclust_order.tsv / gene_colors.tsv :
             exported from RNA_WGCNA114_model.rds (net$dendrograms[[1]], net$colors)
PCA variance (PC1 36.455%, PC2 30.947%) recomputed from the non-grey module
eigengenes in the same .rds (prcomp, centred and scaled); PC scores match the sheet.
"""
import pickle
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.cluster import hierarchy as sch

HERE = Path(__file__).resolve().parent
OUT = HERE / "out"
OUT.mkdir(exist_ok=True)
MM = 1 / 25.4
TEAL, RED, DARK, GREY = "#2A9D8F", "#D95F5F", "#242A30", "#B9C0C7"
TKC, NSC = "#8C8C8C", "#756BB1"

mpl.rcParams.update({
    "font.family": "Arial", "font.size": 6, "axes.labelsize": 6, "axes.titlesize": 6.5,
    "xtick.labelsize": 5.5, "ytick.labelsize": 5.5, "legend.fontsize": 5.5,
    "axes.linewidth": 0.5, "xtick.major.width": 0.5, "ytick.major.width": 0.5,
    "xtick.major.size": 2, "ytick.major.size": 2, "axes.unicode_minus": True,
    "mathtext.fontset": "custom", "mathtext.rm": "Arial", "mathtext.it": "Arial:italic",
    "mathtext.bf": "Arial:bold", "pdf.fonttype": 42,
})
MODCOL = {"turquoise": "#40E0D0", "blue": "#0000FF", "brown": "#A52A2A", "yellow": "#FFFF00",
          "green": "#00FF00", "red": "#FF0000", "black": "#000000", "pink": "#FFC0CB",
          "magenta": "#FF00FF", "purple": "#A020F0", "greenyellow": "#ADFF2F", "grey": "#BEBEBE"}
# darker versions for points/lines where the pure WGCNA colour is too pale on white
MODDOT = dict(MODCOL, yellow="#D4C400", pink="#E7849A", greenyellow="#8CC21F", turquoise="#1FB8A8")
KEY = ["blue", "magenta", "purple", "pink", "turquoise"]


def clean(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def block(df, header_row, stop=None):
    t = df.iloc[header_row + 1:stop].copy()
    t.columns = df.iloc[header_row].tolist()
    t = t.dropna(how="all")
    return t.loc[:, [c for c in t.columns if isinstance(c, str)]]


def stage_label(s):
    return s.replace("d", " d").replace("h", " h")


d = pickle.load(open(HERE / "raw.pkl", "rb"))
pw = block(d["Supp_WGCNA_a_power"], 3).astype({"Power": int})
for c in ["Signed R²", "Mean connectivity"]:
    pw[c] = pw[c].astype(float)
tr = block(d["Supp_WGCNA_c_traits"], 3)
pca = block(d["Supp_WGCNA_d_PCA"], 3, 118)
cen = block(d["Supp_WGCNA_d_PCA"], 120)
ef = d["Supp_WGCNA_e_f"]
corr = ef.iloc[4:15, 0:12].copy()
corr.index = corr.iloc[:, 0]
corr = corr.iloc[:, 1:].astype(float)
corr.columns = ef.iloc[3, 1:12].tolist()
sizes = block(ef, 17).astype({"Genes": int})
gh = d["Supp_WGCNA_g_h"]
traj = block(gh, 3, 1676)
summ = block(gh, 1678, 2515)
hh = block(gh, 2517)

fig = plt.figure(figsize=(180 * MM, 230 * MM))


def letter(x, y, s):
    fig.text(x, y, s, fontsize=9, fontweight="bold", va="top", ha="left")


# ---------------- a: soft threshold ----------------
axa1 = fig.add_axes([0.065, 0.84, 0.155, 0.125])
axa2 = fig.add_axes([0.285, 0.84, 0.155, 0.125])
sel = pw[pw.Power == 22]
axa1.plot(pw.Power, pw["Signed R²"], color=DARK, lw=0.7, marker="o", ms=2.2, zorder=2)
axa1.plot(sel.Power, sel["Signed R²"], "o", ms=4.5, color=RED, zorder=3)
axa1.axhline(0.80, color=GREY, ls="--", lw=0.6)
axa1.text(1, 0.815, "0.80", fontsize=5.2, color="#777777", va="bottom")
axa1.annotate(f"power 22\n$R^2$ = {sel['Signed R²'].iloc[0]:.3f}", xy=(22, sel["Signed R²"].iloc[0]),
              xytext=(13.5, 0.30), fontsize=5.5, color=RED,
              arrowprops=dict(arrowstyle="-", color=RED, lw=0.5))
axa1.set_xlabel("Soft-threshold power")
axa1.set_ylabel("Scale-free fit, signed $R^2$")
axa1.set_ylim(-0.3, 0.9)
axa1.set_xticks([1, 10, 20, 30])
axa1.set_title("Scale-free topology fit", pad=2)
axa2.plot(pw.Power, pw["Mean connectivity"], color=DARK, lw=0.7, marker="o", ms=2.2)
axa2.plot(sel.Power, sel["Mean connectivity"], "o", ms=4.5, color=RED)
axa2.set_yscale("log")
axa2.set_xlabel("Soft-threshold power")
axa2.set_ylabel("Mean connectivity")
axa2.set_xticks([1, 10, 20, 30])
axa2.set_title("Network connectivity", pad=2)
for ax in (axa1, axa2):
    clean(ax)
letter(0.005, 0.985, "a")

# ---------------- b: dendrogram + module colours ----------------
merge = pd.read_csv(HERE / "hclust_merge.tsv", sep="\t", header=None).values
height = pd.read_csv(HERE / "hclust_height.tsv", sep="\t").height.values
colors = pd.read_csv(HERE / "gene_colors.tsv", sep="\t")
n = len(colors)
Z = np.zeros((n - 1, 4))
size = {}
for i, (a, b) in enumerate(merge):
    ia = -a - 1 if a < 0 else n + a - 1
    ib = -b - 1 if b < 0 else n + b - 1
    sa = 1 if a < 0 else size[ia]
    sb = 1 if b < 0 else size[ib]
    Z[i] = [ia, ib, height[i], sa + sb]
    size[n + i] = sa + sb
axb = fig.add_axes([0.52, 0.865, 0.465, 0.1])
dn = sch.dendrogram(Z, ax=axb, no_labels=True, color_threshold=0, above_threshold_color=DARK,
                    link_color_func=lambda k: DARK)
for coll in axb.collections:
    coll.set_linewidth(0.25)
axb.set_ylim(height.min() - 0.02, 1.005)
axb.set_ylabel("Height (1 − TOM)")
axb.set_xticks([])
axb.spines[["top", "right", "bottom"]].set_visible(False)
axb.set_title("Gene dendrogram and merged modules (5,000 genes)", pad=2)
axbc = fig.add_axes([0.52, 0.845, 0.465, 0.016])
leaf_cols = [MODCOL[c] for c in colors.merged.values[dn["leaves"]]]
axbc.imshow(np.array([mpl.colors.to_rgb(c) for c in leaf_cols])[None, :, :], aspect="auto",
            interpolation="nearest")
axbc.set_xticks([])
axbc.set_yticks([0], ["Module"])
for s in axbc.spines.values():
    s.set_linewidth(0.4)
letter(0.49, 0.985, "b")

# ---------------- c: module–trait heat map ----------------
mods = [m for m in corr.index]  # non-grey order used in the sheet
row_order = ["blue", "magenta", "purple", "pink", "turquoise", "brown", "yellow", "red", "black",
             "greenyellow", "green"]
traits = ["TN_vs_FL", "Developmental_day", "Postharvest_hour"]
tlab = ["TN vs FL", "Developmental\nday", "Postharvest\nhour"]
axc = fig.add_axes([0.085, 0.53, 0.175, 0.225])
R = np.zeros((len(row_order), 3))
for i, m in enumerate(row_order):
    for j, t in enumerate(traits):
        r = tr[(tr.Module == m) & (tr.Trait == t)].iloc[0]
        R[i, j] = float(r["Pearson correlation"])
        fdr = float(r.FDR)
        ftxt = "<0.001" if fdr < 0.001 else f"{fdr:.3f}"
        col = "white" if abs(R[i, j]) > 0.6 else DARK
        axc.text(j, i - 0.17, f"{R[i, j]:.2f}".replace("-", "\u2212"), ha="center", va="center", fontsize=5.5, color=col)
        axc.text(j, i + 0.22, ftxt, ha="center", va="center", fontsize=4.9, color=col)
im = axc.imshow(R, cmap="RdBu_r", vmin=-1, vmax=1, aspect="auto")
axc.set_xticks(range(3), tlab)
axc.set_yticks(range(len(row_order)), row_order)
for k, lab in enumerate(axc.get_yticklabels()):
    lab.set_color(DARK)
axc.tick_params(length=0)
for i, m in enumerate(row_order):
    axc.add_patch(mpl.patches.Rectangle((-0.74, i - 0.42), 0.18, 0.84, color=MODCOL[m], clip_on=False,
                                        transform=axc.transData, lw=0.2, ec="#666666"))
axc.tick_params(axis="y", pad=10)
plt.setp(axc.get_xticklabels(), rotation=35, ha="right", rotation_mode="anchor")
axc.set_title("Module–trait correlation", pad=2)
cb = fig.colorbar(im, cax=fig.add_axes([0.085, 0.462, 0.175, 0.007]), orientation="horizontal")
cb.set_label("Pearson $r$ (upper value in each cell); FDR below", fontsize=5.3, labelpad=1)
cb.ax.tick_params(labelsize=5, length=1.5)
cb.outline.set_linewidth(0.4)
letter(0.005, 0.795, "c")

# ---------------- d: PCA of eigengenes ----------------
axd = fig.add_axes([0.33, 0.53, 0.19, 0.225])
pca[["PC1", "PC2"]] = pca[["PC1", "PC2"]].astype(float)
cen[["PC1 centroid", "PC2 centroid", "Stage index"]] = cen[["PC1 centroid", "PC2 centroid", "Stage index"]].astype(float)
for mat, col in (("FL", TEAL), ("TN", RED)):
    s = pca[pca.Material == mat]
    axd.scatter(s.PC1, s.PC2, s=2, color=col, alpha=0.25, lw=0, zorder=1)
    c = cen[cen.Material == mat].sort_values("Stage index")
    axd.plot(c["PC1 centroid"], c["PC2 centroid"], color=col, lw=0.8, zorder=2)
    dv = c[c.Phase == "Development"]
    ph = c[c.Phase == "Postharvest"]
    axd.scatter(dv["PC1 centroid"], dv["PC2 centroid"], s=10, color=col, marker="o", zorder=3,
                edgecolor="white", lw=0.3)
    axd.scatter(ph["PC1 centroid"], ph["PC2 centroid"], s=13, color=col, marker="^", zorder=3,
                edgecolor="white", lw=0.3)
    first = c.iloc[0]
    axd.annotate("0 d", (first["PC1 centroid"], first["PC2 centroid"]), xytext=(3, 3),
                 textcoords="offset points", fontsize=5.2, color=col)
    last = c.iloc[-1]
    axd.annotate("72 h", (last["PC1 centroid"], last["PC2 centroid"]), xytext=(-12, 4),
                 textcoords="offset points", fontsize=5.2, color=col)
axd.set_xlabel("PC1 (36.5%)")
axd.set_ylabel("PC2 (30.9%)")
axd.axhline(0, color="#EEEEEE", lw=0.5, zorder=0)
axd.axvline(0, color="#EEEEEE", lw=0.5, zorder=0)
h = [mpl.lines.Line2D([], [], color=TEAL, lw=0.8, label="FL"),
     mpl.lines.Line2D([], [], color=RED, lw=0.8, label="TN"),
     mpl.lines.Line2D([], [], color=DARK, marker="o", ls="", ms=3, label="Development"),
     mpl.lines.Line2D([], [], color=DARK, marker="^", ls="", ms=3.4, label="Postharvest")]
axd.legend(handles=h, frameon=False, loc="upper center", ncol=2, handlelength=1.4, borderaxespad=0.1, columnspacing=0.8, fontsize=5.2)
axd.set_ylim(axd.get_ylim()[0], axd.get_ylim()[1] + 0.9)
axd.set_title("Eigengene PCA (stage centroids)", pad=2)
clean(axd)
letter(0.285, 0.795, "d")

# ---------------- e: eigengene correlation with clustering ----------------
C = corr.loc[corr.columns, corr.columns]
L = sch.linkage(1 - C.values[np.triu_indices(len(C), 1)], method="average")
axe_d = fig.add_axes([0.62, 0.715, 0.15, 0.04])
dn2 = sch.dendrogram(L, ax=axe_d, no_labels=True, color_threshold=0, link_color_func=lambda k: DARK)
for coll in axe_d.collections:
    coll.set_linewidth(0.5)
axe_d.axis("off")
order = [C.columns[i] for i in dn2["leaves"]]
axe = fig.add_axes([0.62, 0.555, 0.15, 0.16])
im2 = axe.imshow(C.loc[order, order].values, cmap="RdBu_r", vmin=-1, vmax=1, aspect="auto")
axe.set_xticks(range(len(order)), order, rotation=90)
axe.set_yticks(range(len(order)), order)
axe.tick_params(length=0, pad=1)
axe.set_title("")
fig.text(0.695, 0.762, "Eigengene correlation", ha="center", va="bottom", fontsize=6.5)
cb2 = fig.colorbar(im2, cax=fig.add_axes([0.776, 0.555, 0.006, 0.16]))
cb2.ax.tick_params(labelsize=5, length=1.5)
cb2.set_label("Pearson $r$", fontsize=5.5, labelpad=1)
cb2.outline.set_linewidth(0.4)
letter(0.555, 0.795, "e")

# ---------------- f: module sizes ----------------
axf = fig.add_axes([0.905, 0.53, 0.08, 0.225])
sz = sizes[sizes.Module != "grey"].sort_values("Genes")
axf.barh(range(len(sz)), sz.Genes, color=[MODCOL[m] for m in sz.Module], edgecolor="#555555", lw=0.25)
axf.set_yticks(range(len(sz)), sz.Module)
for i, v in enumerate(sz.Genes):
    axf.text(v + 30, i, f"{v:,}", va="center", fontsize=5)
axf.set_xlim(0, 2100)
axf.set_xticks([0, 1000, 2000])
axf.set_xlabel("Genes")
axf.set_title("Module size", pad=2)
axf.tick_params(axis="y", length=0, pad=1)
clean(axf)
letter(0.83, 0.795, "f")

# ---------------- g: trajectories ----------------
summ = summ.astype({"Stage index": int})
for c in ["Mean", "s.e.m."]:
    summ[c] = summ[c].astype(float)
traj = traj.astype({"Stage index": int})
traj["Score"] = traj["Score"].astype(float)
stages = summ[["Stage index", "Stage"]].drop_duplicates().sort_values("Stage index")
tick_idx = [1, 3, 5, 7, 9, 11, 13, 15, 17, 19]
tick_lab = [stage_label(stages.set_index("Stage index").Stage[i]) for i in tick_idx]
gw, gx0, gy0, gh_ = 0.152, 0.065, 0.255, 0.14
axes_g = []
for k, m in enumerate(KEY):
    ax = fig.add_axes([gx0 + k * (gw + 0.0395), gy0, gw, gh_])
    ax.axvspan(13.5, 19.5, color="#F4F4F4", zorder=0, lw=0)
    for mat, col in (("FL", TEAL), ("TN", RED)):
        s = summ[(summ.Material == mat) & (summ.Module == m)].sort_values("Stage index")
        ax.fill_between(s["Stage index"], s.Mean - s["s.e.m."], s.Mean + s["s.e.m."], color=col,
                        alpha=0.3, lw=0, zorder=2)
        ax.plot(s["Stage index"], s.Mean, color=col, lw=0.9, marker="o", ms=1.6, zorder=3)
    for mat, col in (("TK", TKC), ("NS", NSC)):
        s = traj[(traj.Material == mat) & (traj.Module == m)].sort_values("Stage index")
        ax.plot(s["Stage index"], s.Score, color=col, lw=0.7, ls=(0, (2.5, 1.5)), zorder=2)
    ax.set_xlim(0.5, 19.5)
    ax.set_xticks(tick_idx, tick_lab, rotation=90)
    ax.set_title(m, pad=7, color=DARK)
    ax.add_patch(mpl.patches.Rectangle((0, 1.01), 1, 0.035, transform=ax.transAxes, color=MODCOL[m],
                                       clip_on=False, lw=0))
    clean(ax)
    axes_g.append(ax)
axes_g[0].set_ylabel("Standardized module score")
fig.text(0.5, 0.203, "Sampling stage (shaded: postharvest)", ha="center", fontsize=6)
hg = [mpl.lines.Line2D([], [], color=TEAL, lw=0.9, label="FL (mean ± s.e.m., n = 3)"),
      mpl.lines.Line2D([], [], color=RED, lw=0.9, label="TN (mean ± s.e.m., n = 3)"),
      mpl.lines.Line2D([], [], color=TKC, lw=0.7, ls=(0, (2.5, 1.5)), label="TK (projection only)"),
      mpl.lines.Line2D([], [], color=NSC, lw=0.7, ls=(0, (2.5, 1.5)), label="NS (projection only)")]
fig.legend(handles=hg, frameon=False, ncol=4, loc="lower center", bbox_to_anchor=(0.5, 0.415),
           handlelength=2.2, columnspacing=1.6)
letter(0.005, 0.44, "g")

# ---------------- h: kME vs GS ----------------
hh = hh.copy()
hh["show"] = hh["Shown in panel h"].astype(str).isin(["1", "True", "1.0"])
for c in ["Developmental GS", "Assigned-module kME"]:
    hh[c] = hh[c].astype(float)
hy0, hh_ = 0.04, 0.1
hstats = []
for k, m in enumerate(KEY):
    ax = fig.add_axes([gx0 + k * (gw + 0.0395), hy0, gw, hh_])
    s = hh[(hh.Module == m) & hh.show]
    x = s["Assigned-module kME"].abs().values
    y = s["Developmental GS"].abs().values
    ax.scatter(x, y, s=1.4, color=MODDOT[m], alpha=0.55, lw=0, rasterized=True)
    b1, b0 = np.polyfit(x, y, 1)
    xx = np.linspace(x.min(), x.max(), 10)
    ax.plot(xx, b0 + b1 * xx, color=DARK, lw=0.6, ls="--")
    r = np.corrcoef(x, y)[0, 1]
    hstats.append((m, len(s), round(r, 3)))
    ax.set_title(f"{m} (n = {len(s):,})", pad=2)
    ax.set_ylim(0, 1)
    clean(ax)
    if k == 0:
        ax.set_ylabel("|Developmental GS|")
fig.text(0.5, 0.004, "Module membership |kME|", ha="center", va="bottom", fontsize=6)
letter(0.005, 0.172, "h")

for ext in ("pdf", "png"):
    fig.savefig(OUT / f"Supplementary_Fig_03.{ext}", dpi=600, facecolor="white")

# source tables for the new/derived values
sem = summ[summ.Module.isin(KEY)][["Material", "Module", "Stage index", "Stage", "Mean", "s.d.", "Samples (n)", "s.e.m."]]
sem.to_csv(HERE / "SF3g_FL_TN_mean_sem.tsv", sep="\t", index=False)
chk = traj[traj.Material.isin(["FL", "TN"]) & traj.Module.isin(KEY)].groupby(["Material", "Module", "Stage index"]).Score.agg(["mean", "std", "count"])
chk["sem_recalc"] = chk["std"] / np.sqrt(chk["count"])
mg = sem.set_index(["Material", "Module", "Stage index"]).join(chk)
print("max |mean diff|", float((mg.Mean - mg["mean"]).abs().max()), "max |sem diff|", float((mg["s.e.m."] - mg["sem_recalc"]).abs().max()))
pd.DataFrame({"dendrogram_leaf_position": range(1, n + 1), "gene_id": colors.gene.values[dn["leaves"]],
              "merged_module": colors.merged.values[dn["leaves"]]}).to_csv(HERE / "SF3b_leaf_order_modules.tsv", sep="\t", index=False)
pd.DataFrame({"merge_1": merge[:, 0], "merge_2": merge[:, 1], "height": height}).to_csv(HERE / "SF3b_hclust_merge_height.tsv", sep="\t", index=False)
print("h stats", hstats)
print("R2 max", pw["Signed R²"].max(), "any>=0.80", bool((pw["Signed R²"] >= 0.80).any()), "powers", pw.Power.tolist())
print("genes", int(sizes.Genes.sum()), "non-grey", int((sizes.Module != "grey").sum()))
