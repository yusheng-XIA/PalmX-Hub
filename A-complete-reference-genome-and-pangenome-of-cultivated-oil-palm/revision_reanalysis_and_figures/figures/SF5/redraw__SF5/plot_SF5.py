#!/usr/bin/env python3
"""Supplementary Fig. 5 redraw (proteome differences, module preservation, cross-platform concordance).

Data: Source Data SF11a/b/c (REVIEW workbook, via ../data/f8_raw.pkl).
"""
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from f8_style import DARK, GREY, MM, RED, TEAL, clean, letter, sheet, stage_label  # noqa: E402

HERE = Path(__file__).resolve().parent
fig = plt.figure(figsize=(180 * MM, 175 * MM))
gs = fig.add_gridspec(3, 1, left=0.1, right=0.985, top=0.965, bottom=0.075, hspace=0.55,
                      height_ratios=[1, 1.15, 1])

# ---- a: differential proteins per stage ------------------------------------------
de = sheet("SF11a_protein_DE")
de["significant_proteins"] = de.significant_proteins.astype(int)
piv = de.pivot(index="Stage index", columns="direction", values="significant_proteins")
stages = de.drop_duplicates("Stage index").set_index("Stage index").Stage
x = piv.index.values.astype(int)
ax = fig.add_subplot(gs[0])
w = 0.4
ax.bar(x - w / 2, piv["Higher in FL"], width=w, color=TEAL, label="Higher in FL")
ax.bar(x + w / 2, piv["Higher in TN"], width=w, color=RED, label="Higher in TN")
ax.axvline(12.5, color=GREY, ls="--", lw=0.6)
ax.text(6, 1.0, "Development", transform=ax.get_xaxis_transform(), ha="center", va="bottom", fontsize=6,
        color="#555555")
ax.text(15.5, 1.0, "Postharvest", transform=ax.get_xaxis_transform(), ha="center", va="bottom", fontsize=6,
        color="#555555")
ax.set_xticks(x, [stage_label(s) for s in stages], rotation=45, ha="right")
ax.set_xlim(-0.6, 18.6)
ax.set_ylabel("Differentially abundant\nproteins")
ax.set_xlabel("Sampling stage", labelpad=1)
ax.legend(frameon=False, loc="upper right", ncol=2, bbox_to_anchor=(1, 0.97))
clean(ax)
ax_a = ax

# ---- b: module preservation Z-summary ---------------------------------------------
pr = sheet("SF11b_preservation")
pr["Zsummary"] = pr.Zsummary.astype(float)
p = pr.pivot(index="module", columns="direction", values="Zsummary")
p = p.loc[p.max(axis=1).sort_values().index]
ax = fig.add_subplot(gs[1])
y = np.arange(len(p))
for yi, (_, r) in zip(y, p.iterrows()):
    ax.plot([r.FL_to_TN, r.TN_to_FL], [yi, yi], color="#C8C8C8", lw=1.0, zorder=1)
ax.scatter(p.FL_to_TN, y, s=16, color=TEAL, zorder=3, label="FL modules tested in TN", lw=0)
ax.scatter(p.TN_to_FL, y, s=16, color=RED, marker="D", zorder=3, label="TN modules tested in FL", lw=0)
for z in (2, 10):
    ax.axvline(z, color=GREY, ls="--", lw=0.6, zorder=0)
    ax.text(z, len(p) - 0.35, f"Z = {z}", ha="center", va="bottom", fontsize=5.5, color="#555555")
ax.set_yticks(y, [m.capitalize() for m in p.index])
ax.set_ylim(-0.6, len(p) - 0.1)
ax.set_xlim(0, 43)
ax.set_xlabel(r"Module-preservation $Z_{\mathrm{summary}}$")
ax.set_ylabel("Protein module")
ax.legend(frameon=False, loc="lower right")
clean(ax)
ax_b = ax

# ---- c: Astral vs timsTOF concordance --------------------------------------------
cc = sheet("SF11c_concordance")
cc["rho"] = cc.spearman_correlation_across_OAU.astype(float)
cc["Stage index"] = cc["Stage index"].astype(int)
ax = fig.add_subplot(gs[2])
spec = [("Astral vs timsTOF, FL", "FL abundance", TEAL, "o"),
        ("Astral vs timsTOF, TN", "TN abundance", RED, "s"),
        ("Astral vs timsTOF, TN − FL", "TN − FL difference", DARK, "D")]
rows = []
for comp, lab, col, mk in spec:
    s = cc[cc.comparison == comp].sort_values("Stage index")
    ax.plot(s["Stage index"], s.rho, color=col, lw=1.0, marker=mk, ms=2.6, label=lab)
    rows.append((lab, round(s.rho.min(), 3), round(s.rho.median(), 3), round(s.rho.max(), 3),
                 int(s.shared_allele_units.astype(int).min()), int(s.shared_allele_units.astype(int).max())))
ax.axhline(0, color=GREY, lw=0.5)
ax.axvline(12.5, color=GREY, ls="--", lw=0.6)
st = cc[cc.comparison == spec[0][0]].sort_values("Stage index")
ax.set_xticks(st["Stage index"], [stage_label(s) for s in st.Stage], rotation=45, ha="right")
ax.set_xlim(-0.6, 18.6)
ax.set_ylim(-0.1, 0.85)
ax.set_ylabel(r"Astral–timsTOF Spearman $\rho$")
ax.set_xlabel("Sampling stage", labelpad=1)
ax.legend(frameon=False, loc="center right", ncol=1, bbox_to_anchor=(1, 0.4))
clean(ax)

for a, s in [(ax_a, "a"), (ax_b, "b"), (ax, "c")]:
    letter(fig, a, s, dx=-0.09)
for ext in ("pdf", "png"):
    fig.savefig(HERE / f"Supplementary_Fig_05.{ext}", dpi=600, facecolor="white")

tot = piv.sum(axis=1)
print("DE FL range", piv["Higher in FL"].min(), piv["Higher in FL"].max(), "TN range", piv["Higher in TN"].min(),
      piv["Higher in TN"].max(), "stages TN>FL", int((piv["Higher in TN"] > piv["Higher in FL"]).sum()))
print(p.round(2).to_string())
print("concordance (label, min, median, max, OAU min, OAU max)", rows)
