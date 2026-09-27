#!/usr/bin/env python3
"""Supplementary Fig. 6 redraw (proteome and metabolome QC).

Every panel is drawn from tabular data rather than edited pixels:
  a  data/SF6a_enzyme_family_OAU_coverage.tsv (values recovered exactly from the
     vector panel F2_SUPP_08_protein_OAU_coverage_QC.svg; all implied counts are integers)
  b   SF6/SF6b_{pos,neg}_PCA_current_batch.tsv: current Waters batch (114 FL/TN biological
      injections + 14 pooled QC; same batch as Supp Fig. 4). Pipeline PCA: log2 after
      half-minimum imputation, Pareto scaling, on the final QC-filtered, drift-corrected matrix.
  c-f Source Data sheets SF12c-SF12f (OilPalm_Source_Data_submission_114RNAseq.xlsx)
Fonts: Arial with a true minus sign (U+2212), TrueType embedded (pdf.fonttype 42).
"""
import pickle
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
OUT = HERE / "out"
OUT.mkdir(exist_ok=True)

TEAL, RED, DARK, GREY, PURPLE, BLUE = "#2A9D8F", "#D95F5F", "#242A30", "#B9C0C7", "#756BB1", "#3C78A8"
MM = 1 / 25.4

mpl.rcParams.update({
    "font.family": "Arial",
    "font.size": 7,
    "axes.labelsize": 7,
    "axes.titlesize": 7,
    "xtick.labelsize": 6,
    "ytick.labelsize": 6,
    "legend.fontsize": 6,
    "axes.linewidth": 0.6,
    "xtick.major.width": 0.6,
    "ytick.major.width": 0.6,
    "xtick.major.size": 2.5,
    "ytick.major.size": 2.5,
    "axes.unicode_minus": True,
    "mathtext.fontset": "custom",
    "mathtext.rm": "Arial",
    "mathtext.it": "Arial:italic",
    "mathtext.bf": "Arial:bold",
    "pdf.fonttype": 42,
    "svg.fonttype": "none",
})


def clean(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_axisbelow(True)


def letter(fig, ax, s, dx=-0.0, dy=0.0):
    bb = ax.get_position()
    fig.text(bb.x0 + dx, bb.y1 + dy, s, fontsize=9, fontweight="bold", va="bottom", ha="left")


d = pickle.load(open(HERE / "data/sf12.pkl", "rb"))
a = pd.read_csv(HERE / "data/SF6a_enzyme_family_OAU_coverage.tsv", sep="\t")

fig = plt.figure(figsize=(180 * MM, 170 * MM))
gs0 = fig.add_gridspec(1, 6, left=0.135, right=0.86, top=0.955, bottom=0.715, wspace=1.3)
gs1 = fig.add_gridspec(1, 6, left=0.075, right=0.985, top=0.63, bottom=0.40, wspace=1.05)
gs2 = fig.add_gridspec(1, 6, left=0.075, right=0.985, top=0.29, bottom=0.07, wspace=1.05)

# ---- a: enzyme-family protein coverage (FL vs TN) --------------------------------
ax = fig.add_subplot(gs0[0, 0:2])
y = np.arange(len(a))[::-1]
for yi, (_, r) in zip(y, a.iterrows()):
    ax.plot([r.FL_quantified_pct, r.TN_quantified_pct], [yi, yi], color="#C8C8C8", lw=1.0, zorder=1)
ax.scatter(a.FL_quantified_pct, y, s=13, color=TEAL, zorder=3, label="FL", lw=0)
ax.scatter(a.TN_quantified_pct, y, s=11, color=RED, marker="s", zorder=2, label="TN", lw=0)
ax.set_yticks(y, [f"{n} ({q}/{m})" for n, q, m in zip(a.enzyme_family, a.FL_quantified_OAUs, a.mapped_OAUs)],
              fontsize=5.5)
ax.set_xlim(-4, 104)
ax.set_xticks([0, 25, 50, 75, 100])
ax.set_xlabel("Mapped OAUs quantified (%)")
ax.set_ylabel("Enzyme family (FL quantified/mapped OAUs)", fontsize=6)
ax.grid(axis="x", color="#EEEEEE", lw=0.5)
ax.legend(frameon=False, loc="lower right", handletextpad=0.2, borderaxespad=0.1)
ax.set_ylim(-0.7, len(a) - 0.3)
clean(ax)
ax_a = ax

# ---- b: metabolome PCA, positive and negative mode --------------------------------
styles = {
    "FL": dict(color=TEAL, marker="o", s=9, label="FL"),
    "TN": dict(color=RED, marker="o", s=9, label="TN"),
    "Analysis QC": dict(color=DARK, marker="^", s=12, label="Analysis QC"),
    "Conditioning QC": dict(color=GREY, marker="s", s=12, label="Conditioning QC"),
}
axes_b = []
for k, mode in enumerate(["pos", "neg"]):
    ax = fig.add_subplot(gs0[0, 2 + 2 * k:4 + 2 * k])
    p = pd.read_csv(HERE / f"SF6/SF6b_{mode}_PCA_current_batch.tsv", sep="\t")
    v = pd.read_csv(HERE / f"SF6/SF6b_{mode}_PCA_variance_current_batch.tsv", sep="\t").set_index(
        "component").variance_percent
    for cls, st in styles.items():
        s = p[p.sample_class == cls]
        st = dict(st, label=f"{st['label']} ({len(s)})")
        ax.scatter(s.PC1, s.PC2, edgecolor="white", linewidth=0.3, zorder={"Analysis QC": 5, "Conditioning QC": 4}.get(cls, 2), **st)
    ax.set_xlabel(f"PC1 ({v['PC1']:.1f}%)")
    ax.set_ylabel(f"PC2 ({v['PC2']:.1f}%)")
    ax.set_title("Positive-ion mode" if mode == "pos" else "Negative-ion mode", pad=3)
    ax.axhline(0, color="#EEEEEE", lw=0.5, zorder=0)
    ax.axvline(0, color="#EEEEEE", lw=0.5, zorder=0)
    clean(ax)
    axes_b.append(ax)
axes_b[1].legend(frameon=False, loc="upper left", bbox_to_anchor=(1.0, 1.0), handletextpad=0.1, borderaxespad=0.2,
                 labelspacing=0.35, fontsize=5.5)

# ---- c: protein groups identified per stage -----------------------------------
ax = fig.add_subplot(gs1[0, 0:4])
c = d["SF12c_identification"]
order = c[["Stage", "Stage index"]].drop_duplicates().sort_values("Stage index")
for g, col in (("FL", TEAL), ("TN", RED)):
    s = c[c.genotype == g]
    ax.scatter(s["Stage index"] + (-0.08 if g == "FL" else 0.08), s.protein_group_n, s=5, color=col,
               alpha=0.45, lw=0, zorder=2)
    m = s.groupby("Stage index").protein_group_n.mean()
    ax.plot(m.index, m.values, color=col, lw=1.1, marker="o", ms=2.6, label=g, zorder=3)
ax.axvline(13.5, color=GREY, lw=0.6, ls="--", zorder=0)
ax.text(7, 1.0, "Development", transform=ax.get_xaxis_transform(), ha="center", va="bottom",
        fontsize=6, color="#555555")
ax.text(16.5, 1.0, "Postharvest", transform=ax.get_xaxis_transform(), ha="center", va="bottom",
        fontsize=6, color="#555555")
ax.set_xticks(order["Stage index"], [s.replace("d", " d").replace("h", " h") for s in order.Stage],
              rotation=45, ha="right")
ax.set_xlim(0.4, 19.6)
ax.set_ylabel("Protein groups identified")
ax.set_xlabel("Sampling stage")
ax.legend(frameon=False, loc="upper right", ncol=2, bbox_to_anchor=(1.0, 0.99))
clean(ax)
ax_c = ax

# ---- d: DIA-NN corrected median mass error --------------------------------------
ax = fig.add_subplot(gs1[0, 4:6])
me = d["SF12d_mass_error"]
vals = [me.Median_Mass_Acc_MS1_Corrected_ppm.values, me.Median_Mass_Acc_MS2_Corrected_ppm.values]
bp = ax.boxplot(vals, positions=[1, 2], widths=0.5, showfliers=False, patch_artist=True,
                medianprops=dict(color=DARK, lw=1.0), whiskerprops=dict(lw=0.6), capprops=dict(lw=0.6),
                boxprops=dict(lw=0.6))
for patch, col in zip(bp["boxes"], ["#94B5CF", "#8AC9C2"]):
    patch.set_facecolor(col)
rng = np.random.default_rng(1)
for i, v in enumerate(vals, start=1):
    ax.scatter(i + rng.uniform(-0.16, 0.16, len(v)), v, s=3, color=DARK, alpha=0.35, lw=0, zorder=3)
ax.set_xticks([1, 2], ["MS1", "MS2"])
ax.set_xlim(0.4, 2.6)
ax.set_ylabel("Corrected median mass error (ppm)")
ax.set_title(f"DIA-NN, {len(me)} Astral acquisitions", pad=3)
clean(ax)
ax_d = ax

# ---- e: missing-value structure -------------------------------------------------
pa = d["SF12a_protein_detect"]
se = d["SF12e_sample_complete"]
ax = fig.add_subplot(gs2[0, 0:2])
n_total = int(se.total_protein_groups.iloc[0])
counts = pa.quantified_samples.value_counts().reindex(range(1, 115), fill_value=0)
ax.bar(counts.index, counts.values, width=1.0, color=BLUE, lw=0)
full = int(counts.loc[114])
ax.set_xlabel("Acquisitions quantifying the protein group")
ax.set_ylabel("Protein groups")
ax.set_title(f"{n_total:,} protein groups", pad=3)
ax.annotate(f"quantified in all 114\n({full:,} groups)", xy=(114, full), xytext=(88, full * 0.93),
            fontsize=5.5, color=RED, ha="right", va="top",
            arrowprops=dict(arrowstyle="-", color=RED, lw=0.5))
ax.set_xlim(0, 116)
clean(ax)
ax_e1 = ax

ax = fig.add_subplot(gs2[0, 2:4])
comp = se.completeness_fraction * 100
ax.hist(comp, bins=np.arange(38, 77, 2), color=TEAL, alpha=0.85, lw=0.3, edgecolor="white")
med = comp.median()
ax.axvline(med, color=DARK, ls="--", lw=0.7)
ax.set_ylim(0, ax.get_ylim()[1] * 1.12)
ax.text(med - 0.6, ax.get_ylim()[1] * 0.98, f"median {med:.0f}%", ha="right", va="top", fontsize=5.5)
ax.set_xlabel("Per-acquisition completeness (%)")
ax.set_ylabel("Acquisitions")
ax.yaxis.set_major_locator(mpl.ticker.MaxNLocator(integer=True))
ax.set_title(f"{len(se)} acquisitions", pad=3)
clean(ax)
ax_e2 = ax

# ---- f: within-group biological-replicate correlations -------------------------
ax = fig.add_subplot(gs2[0, 4:6])
f = d["SF12f_replicates"]
r = f[[c for c in f.columns if c.startswith("pearson")][0]]
ax.hist(r, bins=np.arange(0.975, 0.9905, 0.0015), color=PURPLE, alpha=0.85, lw=0.3, edgecolor="white")
ax.axvline(r.median(), color=DARK, ls="--", lw=0.7)
ax.set_ylim(0, ax.get_ylim()[1] * 1.15)
ax.text(r.median() - 0.0003, ax.get_ylim()[1] * 0.98, f"median $r$ = {r.median():.3f}", ha="right", va="top",
        fontsize=5.5)
ax.set_xlabel(r"Pearson $r$ (log$_2$ abundance)")
ax.set_ylabel("Replicate pairs")
ax.set_title(f"{len(r)} biological-replicate pairs", pad=3)
ax.xaxis.set_major_locator(mpl.ticker.MultipleLocator(0.005))
clean(ax)
ax_f = ax

for ax, s in [(ax_a, "a"), (axes_b[0], "b"), (ax_c, "c"), (ax_d, "d"), (ax_e1, "e"), (ax_f, "f")]:
    letter(fig, ax, s, dx=-0.125 if s == "a" else -0.055, dy=0.012)

for ext in ("pdf", "png"):
    fig.savefig(OUT / f"Supplementary_Fig_06.{ext}", dpi=600, facecolor="white")
print("median r", round(r.median(), 6), "min r", round(r.min(), 6), "n", len(r))
print("median completeness", round(med, 2), "full groups", full, "total", n_total)
