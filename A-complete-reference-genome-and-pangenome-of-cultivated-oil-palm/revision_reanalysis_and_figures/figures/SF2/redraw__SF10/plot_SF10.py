#!/usr/bin/env python3
"""Supplementary Fig. 10: core lipid-pathway gene dosage across palms.

a  curated gene-locus counts (REVIEW Source Data sheet Lipid_gene_dosage = old Fig.2d_copy_number, identical to
   03_V3/02_figure/05_Figure2_evolution_multiomics_panels_flat_20260727/
   Fig2c_core_pathway_enzyme_copy_number_long.tsv)
b  CAFE5-inferred gains/losses at the oil-palm ancestor node for the linked
   orthogroups (03_V3/02_figure/03_Figure2_panels_flat_20260727/
   Fig2b_core_lipid_copy_number_rigorous_v2_ancestral_changes.tsv)
The two panels use different units (curated loci vs OrthoFinder orthogroups) and are
kept separate on purpose.
"""
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap, BoundaryNorm

HERE = Path(__file__).resolve().parent
MM = 1 / 25.4
DARK, GREY = "#242A30", "#B9C0C7"
GAIN, LOSS = "#D95F5F", "#3C78A8"

mpl.rcParams.update({
    "font.family": "Arial", "font.size": 7, "axes.labelsize": 7, "axes.titlesize": 7,
    "xtick.labelsize": 6, "ytick.labelsize": 6, "legend.fontsize": 6,
    "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6,
    "xtick.major.size": 2.5, "ytick.major.size": 2.5, "axes.unicode_minus": True,
    "mathtext.fontset": "custom", "mathtext.rm": "Arial", "mathtext.it": "Arial:italic",
    "mathtext.bf": "Arial:bold", "pdf.fonttype": 42, "svg.fonttype": "none",
})

GENOMES = ["Calamus", "Daemonorops", "Nypa_fruticans", "Phoenix_dactylifera", "Areca_catechu",
           "Cocos_nucifera", "American_hap1", "Dura", "Pisifera"]
LABELS = [r"$\it{Calamus}$", r"$\it{Daemonorops}$", r"$\it{Nypa\ fruticans}$",
          r"$\it{Phoenix\ dactylifera}$", r"$\it{Areca\ catechu}$", r"$\it{Cocos\ nucifera}$",
          r"$\it{E.\ oleifera}$-derived hap1", r"$\it{E.\ guineensis}$ dura",
          r"$\it{E.\ guineensis}$ pisifera"]
FOCAL = ["Cocos_nucifera", "American_hap1", "Dura", "Pisifera"]
SECTION_SHORT = {
    "Plastid fatty-acid synthesis": "Plastid FA synthesis",
    "FA export and modification": "FA export",
    "VLCFA side branch": "VLCFA branch",
    "ER TAG assembly": "ER TAG assembly",
    "PC–DAG exchange and desaturation": "PC–DAG exchange,\ndesaturation",
}

w = pd.read_csv(HERE / "_raw_Lipid_gene_dosage_REVIEW.tsv", sep="\t")  # REVIEW (16:36) sheet Lipid_gene_dosage; identical to old Fig.2d_copy_number
anc = pd.read_csv(HERE / "ref/Fig2b_core_lipid_copy_number_rigorous_v2_ancestral_changes.tsv", sep="\t")

# Source Data export (panel a wide table + identity flag; panel b table)
out = w[["Pathway_section", "Enzyme", "Scope_note", "N_Arabidopsis_anchor_loci", "N_linked_Orthogroups"] + GENOMES].copy()
out["identical_Cocos_and_3_Elaeis"] = w[FOCAL].nunique(axis=1).eq(1)
out.to_csv(HERE / "SF10_copy_number.tsv", sep="\t", index=False)
anc.to_csv(HERE / "SF10b_ancestral_changes.tsv", sep="\t", index=False)

M = w[GENOMES].to_numpy()
n_rows, n_cols = M.shape
same = w[FOCAL].nunique(axis=1).eq(1).to_numpy()

fig = plt.figure(figsize=(180 * MM, 165 * MM))
ax = fig.add_axes([0.30, 0.19, 0.39, 0.745])
bounds = [-0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 6.5, 10.5, 30]
cmap = LinearSegmentedColormap.from_list("cn", ["#F4F6F7", "#DDEBF1", "#B9D5E4", "#8DB9D3",
                                               "#5E98C0", "#3C78A8", "#2B5C87", "#1D3F5E"])
norm = BoundaryNorm(bounds, cmap.N)
im = ax.imshow(M, cmap=cmap, norm=norm, aspect="auto")
for i in range(n_rows):
    for j in range(n_cols):
        v = M[i, j]
        ax.text(j, i, "ND" if v == 0 else str(v), ha="center", va="center", fontsize=5.8,
                color="white" if v >= 5 else DARK)
ax.set_xticks(range(n_cols))
ax.set_xticklabels(LABELS, rotation=40, ha="right", rotation_mode="anchor", fontsize=6)
ax.set_yticks(range(n_rows))
ax.set_yticklabels([e.replace(" (VLCFA side branch)", "").replace(" (omega-6 FAD family)", " (ω6)")
                    .replace(" (omega-3 FAD family)", " (ω3)") for e in w.Enzyme], fontsize=6)
ax.tick_params(length=0, pad=2)
for s in ax.spines.values():
    s.set_visible(False)
ax.set_xticks(np.arange(-0.5, n_cols), minor=True)
ax.set_yticks(np.arange(-0.5, n_rows), minor=True)
ax.grid(which="minor", color="white", lw=0.8)
ax.tick_params(which="minor", length=0)
# Elaeis bracket and focal outline
ax.plot([5.5, 5.5], [-0.5, n_rows - 0.5], color=DARK, lw=0.6, ls=(0, (2, 1.5)))
ax.annotate("", xy=(5.6, -0.95), xytext=(8.4, -0.95), xycoords="data",
            arrowprops=dict(arrowstyle="-", lw=0.7, color=DARK), annotation_clip=False)
ax.text(7.0, -1.15, r"$\it{Elaeis}$", ha="center", va="bottom", fontsize=6.5)
# identity marks on the right
for i in range(n_rows):
    ax.text(n_cols - 0.35, i, "=" if same[i] else "≠", ha="left", va="center", fontsize=6.5,
            color="#7A8288" if same[i] else "#B43A32", fontweight="bold")
# pathway section brackets on the left
sections = w.Pathway_section.tolist()
start = 0
for i in range(1, n_rows + 1):
    if i == n_rows or sections[i] != sections[start]:
        y0, y1 = start - 0.4, i - 0.6
        tr = mpl.transforms.blended_transform_factory(ax.transAxes, ax.transData)
        x = -0.325
        ax.plot([x, x], [y0, y1], color=DARK, lw=0.7, clip_on=False, transform=tr)
        ax.text(x - 0.012, (y0 + y1) / 2, SECTION_SHORT[sections[start]], ha="right", va="center",
                fontsize=5.8, color="#38434A", clip_on=False, transform=tr)
        if i < n_rows:
            ax.axhline(i - 0.5, color=DARK, lw=0.5)
        start = i
ax.set_xlim(-0.5, n_cols - 0.5)
ax.set_ylim(n_rows - 0.5, -0.5)
cax = fig.add_axes([0.335, 0.035, 0.30, 0.014])
cb = fig.colorbar(im, cax=cax, orientation="horizontal", ticks=[0, 1, 2, 3, 4, 5.5, 8.5, 20])
cb.ax.set_xticklabels(["ND", "1", "2", "3", "4", "5–6", "7–10", ">10"], fontsize=5.5)
cb.outline.set_linewidth(0.4)
cb.ax.tick_params(length=1.5, width=0.4)
cax.set_title("Curated gene loci per genome", fontsize=6, pad=2)

# ---- b: CAFE5 ancestral gains/losses -------------------------------------------
axb = fig.add_axes([0.835, 0.19, 0.15, 0.745])
yb = np.arange(len(anc))
axb.barh(yb, anc.Total_inferred_gains, color=GAIN, height=0.55, lw=0)
axb.barh(yb, -anc.Total_inferred_losses, color=LOSS, height=0.55, lw=0)
for yi, (_, r) in zip(yb, anc.iterrows()):
    if r.Total_inferred_gains:
        axb.text(r.Total_inferred_gains + 0.1, yi, f"+{r.Total_inferred_gains}", va="center", ha="left",
                 fontsize=5.5, color=GAIN)
    if r.Total_inferred_losses:
        axb.text(-r.Total_inferred_losses - 0.1, yi, f"−{r.Total_inferred_losses}", va="center",
                 ha="right", fontsize=5.5, color=LOSS)
axb.axvline(0, color=DARK, lw=0.6)
axb.set_yticks(yb, [f"{e} ({n})" for e, n in zip(anc.Enzyme, anc.N_linked_OGs)], fontsize=5.8)
axb.set_ylim(len(anc) - 0.5, -0.5)
axb.set_xlim(-3.2, 3.2)
axb.set_xticks([-2, 0, 2], ["−2", "0", "+2"])
axb.set_xlabel("Inferred copies at the\noil-palm ancestor node", fontsize=6)
axb.set_title("Losses | gains", fontsize=6, pad=3)
axb.tick_params(axis="y", length=0, pad=2)
for s in ("top", "right", "left"):
    axb.spines[s].set_visible(False)

fig.text(0.012, 0.965, "a", fontsize=9, fontweight="bold", va="top")
fig.text(0.735, 0.965, "b", fontsize=9, fontweight="bold", va="top")
for ext in ("pdf", "png"):
    fig.savefig(HERE / f"Supplementary_Fig_10.{ext}", dpi=600, facecolor="white")
print("identical rows (Cocos + 3 Elaeis):", int(same.sum()), "of", n_rows)
