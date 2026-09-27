#!/usr/bin/env python3
"""Supplementary Fig. 5 with OLE16 locus panels (d-i) added below the existing a-c.

a-c: same data and encodings as fix/beautify/work/SF5/sf5/plot_SF5.py (scheme A colours), re-laid out
(b and c side by side) to free space; d-i: OLE16a/OLE16b loci from build_sf5_data.py (sd/*.tsv).
Outputs out/Supplementary_Fig_05.{pdf,png} (600 dpi) and out/docx_2400/Supplementary_Fig_05.png."""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from Bio import Phylo
from PIL import Image

HERE = Path(__file__).resolve().parent
S = HERE.parents[2]
sys.path.insert(0, str(S / "fix/beautify/work/SF5"))
import f8_style as st  # noqa: E402  (scheme-A TEAL/RED, rcParams, sheet(), stage_label())
from f8_style import DARK, GREY, MM, clean, sheet, stage_label  # noqa: E402

FL, TN = st.TEAL, st.RED
SD = HERE / "sd"
OUT = HERE / "out"
(OUT / "docx_2400").mkdir(parents=True, exist_ok=True)
ST = ["0d", "15d", "35d", "50d", "65d", "80d", "95d", "110d", "125d", "140d", "155d", "170d", "185d",
      "12h", "24h", "36h", "48h", "60h", "72h"]

fig = plt.figure(figsize=(180 * MM, 232 * MM))
LETTERS = []


def letter(ax, s, dx=-0.06, dy=0.008):
    LETTERS.append((ax, s, dx, dy))


def stage_axis(ax, stages, rot=45):
    ax.set_xticks(range(len(stages)), [stage_label(s) for s in stages], rotation=rot, ha="right")


# ================================ a-c (unchanged content) ================================
top = fig.add_gridspec(2, 2, left=0.085, right=0.985, top=0.975, bottom=0.64, hspace=0.75, wspace=0.42,
                       width_ratios=[0.8, 1.45], height_ratios=[1, 1.15])
de = sheet("SF11a_protein_DE")
de["significant_proteins"] = de.significant_proteins.astype(int)
piv = de.pivot(index="Stage index", columns="direction", values="significant_proteins")
stages = de.drop_duplicates("Stage index").set_index("Stage index").Stage
x = piv.index.values.astype(int)
ax = fig.add_subplot(top[0, :])
w = 0.4
ax.bar(x - w / 2, piv["Higher in FL"], width=w, color=FL, label="Higher in FL")
ax.bar(x + w / 2, piv["Higher in TN"], width=w, color=TN, label="Higher in TN")
ax.axvline(12.5, color=GREY, ls="--", lw=0.6)
ax.text(6, 1.0, "Development", transform=ax.get_xaxis_transform(), ha="center", va="bottom", fontsize=6, color="#555555")
ax.text(15.5, 1.0, "Post-harvest", transform=ax.get_xaxis_transform(), ha="center", va="bottom", fontsize=6,
        color="#555555")
ax.set_xticks(x, [stage_label(s) for s in stages], rotation=45, ha="right")
ax.set_xlim(-0.6, 18.6)
ax.set_ylim(0, 2900)
ax.set_ylabel("Differentially\nabundant proteins")
ax.legend(frameon=False, loc="upper right", ncol=2, bbox_to_anchor=(1, 0.97))
clean(ax)
letter(ax, "a", dx=-0.075)

pr = sheet("SF11b_preservation")
pr["Zsummary"] = pr.Zsummary.astype(float)
p = pr.pivot(index="module", columns="direction", values="Zsummary")
p = p.loc[p.max(axis=1).sort_values().index]
ax = fig.add_subplot(top[1, 0])
y = np.arange(len(p))
for yi, (_, r) in zip(y, p.iterrows()):
    ax.plot([r.FL_to_TN, r.TN_to_FL], [yi, yi], color="#C8C8C8", lw=1.0, zorder=1)
ax.scatter(p.FL_to_TN, y, s=12, color=FL, zorder=3, label="FL reference, TN test", lw=0)
ax.scatter(p.TN_to_FL, y, s=12, color=TN, marker="D", zorder=3, label="TN reference, FL test", lw=0)
for z in (2, 10):
    ax.axvline(z, color=GREY, ls="--", lw=0.6, zorder=0)
    ax.text(z, len(p) - 0.35, f"$Z$ = {z}", ha="center", va="bottom", fontsize=5.5, color="#555555")
ax.set_yticks(y, [m.capitalize() for m in p.index])
ax.set_ylim(-0.6, len(p) - 0.1)
ax.set_xlim(0, 43)
ax.set_xlabel(r"Module-preservation $Z_{\mathrm{summary}}$")
ax.legend(frameon=False, loc="lower right", fontsize=5.5, handletextpad=0.2, borderaxespad=0.1)
clean(ax)
letter(ax, "b", dx=-0.075)

cc = sheet("SF11c_concordance")
cc["rho"] = cc.spearman_correlation_across_OAU.astype(float)
cc["Stage index"] = cc["Stage index"].astype(int)
ax = fig.add_subplot(top[1, 1])
for comp, lab, col, mk in [("Astral vs timsTOF, FL", "FL abundance", FL, "o"),
                           ("Astral vs timsTOF, TN", "TN abundance", TN, "s"),
                           ("Astral vs timsTOF, TN − FL", "TN − FL difference", DARK, "D")]:
    s = cc[cc.comparison == comp].sort_values("Stage index")
    ax.plot(s["Stage index"], s.rho, color=col, lw=0.9, marker=mk, ms=2.2, label=lab)
ax.axhline(0, color=GREY, lw=0.5)
ax.axvline(12.5, color=GREY, ls="--", lw=0.6)
s0 = cc[cc.comparison == "Astral vs timsTOF, FL"].sort_values("Stage index")
ax.set_xticks(s0["Stage index"], [stage_label(s) for s in s0.Stage], rotation=45, ha="right")
ax.set_xlim(-0.6, 18.6)
ax.set_ylim(-0.1, 0.85)
ax.set_ylabel(r"Astral–timsTOF Spearman $\rho$")
ax.legend(frameon=False, loc="center right", bbox_to_anchor=(1, 0.42), fontsize=5.5)
clean(ax)
letter(ax, "c", dx=-0.055)

# ================================ d-i (OLE16 loci) ================================
bot = fig.add_gridspec(3, 3, left=0.085, right=0.985, top=0.585, bottom=0.045, hspace=0.62, wspace=0.5,
                       width_ratios=[1.25, 1, 1])
# ---- d: pruned phylogeny ----
ax = fig.add_subplot(bot[:, 0])
t = Phylo.read(str(HERE / "tree_pruned.nwk"), "newick")
for c in t.get_nonterminals():
    c.name = None
col = {}
for c in t.get_terminals():
    col[c.name] = "#111111" if c.name.startswith("OLE16a") else ("#555555" if c.name.startswith("OLE16b") else "#333333")
depths = t.depths()
if not max(depths.values()):
    depths = t.depths(unit_branch_lengths=True)
terms = t.get_terminals()
ypos = {c: i for i, c in enumerate(terms)}


def ycoord(c):
    if c in ypos:
        return ypos[c]
    v = np.mean([ycoord(k) for k in c.clades])
    ypos[c] = v
    return v


ycoord(t.root)
for c in t.find_clades(order="preorder"):
    for k in c.clades:
        ax.plot([depths[c], depths[c]], [ypos[c], ypos[k]], color="#444444", lw=0.6)
        ax.plot([depths[c], depths[k]], [ypos[k], ypos[k]], color="#444444", lw=0.6)
    if c.clades and c.confidence is not None and c is not t.root:
        if c.confidence >= 95:
            ax.plot(depths[c], ypos[c], "o", ms=2.0, color="#222222", zorder=4)
        elif c.confidence >= 70:
            ax.plot(depths[c], ypos[c], "o", ms=2.0, mfc="white", mec="#222222", mew=0.5, zorder=4)
xmax = max(depths.values())
for c in terms:
    ax.text(depths[c] + xmax * 0.02, ypos[c], c.name, fontsize=5, va="center", color=col[c.name],
            fontweight="bold" if c.name.startswith("OLE16") else "normal")
ax.set_ylim(len(terms) + 1.6, -0.8)
ax.set_xlim(-xmax * 0.02, xmax * 1.75)
sb = 0.5
yb = len(terms) + 0.2
ax.plot([0, sb], [yb, yb], color="#222222", lw=0.8)
ax.text(sb / 2, yb - 0.25, f"{sb}", ha="center", va="bottom", fontsize=5)
ax.plot([xmax * 0.35], [yb], "o", ms=2.0, color="#222222")
ax.text(xmax * 0.38, yb, "UFBoot ≥ 95", fontsize=5, va="center")
ax.plot([xmax * 0.9], [yb], "o", ms=2.0, mfc="white", mec="#222222", mew=0.5)
ax.text(xmax * 0.93, yb, "UFBoot 70–94", fontsize=5, va="center")
ax.axis("off")
letter(ax, "d", dx=-0.075, dy=0.004)

# ---- e: bulk RNA ----
e = pd.read_csv(SD / "SF5e_OLE16_RNA.tsv", sep="\t")
vcol = "Normalized count (DESeq2 size factors, 114 libraries)"
ax = fig.add_subplot(bot[0, 1:])
for loc, ls, mk in [("OLE16a", "-", "o"), ("OLE16b", "--", "^")]:
    for m, c in [("FL", FL), ("TN", TN)]:
        d = e[(e.Locus == loc) & (e.Material == m)]
        mm = d.groupby("Stage index")[vcol].mean()
        ax.plot(mm.index - 1, np.log10(mm + 1), ls=ls, color=c, lw=0.9, label=f"{loc} {m}")
        ax.scatter(d["Stage index"] - 1, np.log10(d[vcol] + 1), s=3, marker=mk, color=c, lw=0, alpha=0.6)
ax.axvline(12.5, color=GREY, ls="--", lw=0.6)
stage_axis(ax, ST)
ax.set_xlim(-0.6, 18.6)
ax.set_ylabel("log$_{10}$(normalized\ncount + 1)")
ax.legend(frameon=False, ncol=2, fontsize=5.5, loc="upper left", handlelength=2.2, columnspacing=1)
clean(ax)
letter(ax, "e", dx=-0.07)

# ---- f: shared peptides ----
f = pd.read_csv(SD / "SF5f_OLE16_shared_peptides.tsv", sep="\t")
Q = pd.read_csv(SD / "SF5f_run_precursor_quantiles.tsv", sep="\t")
PS = ST[10:]
ax = fig.add_subplot(bot[1, 1])
for prec, loc, ls, mk in [("RVPGSEQLEQAR3", "OLE16a", "-", "o"), ("VPGSEQLEQAR2", "OLE16a", ":", "s"),
                          ("RPPGFEQLEQAR3", "OLE16b", "--", "^")]:
    for m, c in [("FL", FL), ("TN", TN)]:
        d = f[(f["Precursor.Id"] == prec) & (f.Material == m)]
        xi = d.Stage.map({s: i for i, s in enumerate(PS)})
        ax.scatter(xi, np.log10(d["Precursor.Quantity"]), s=3, marker=mk, color=c, lw=0, alpha=0.7)
        mm = d.groupby("Stage")["Precursor.Quantity"].mean().reindex(PS)
        ax.plot(range(len(PS)), np.log10(mm), ls=ls, color=c, lw=0.8)
med = np.log10(Q.p50.median())
ax.axhline(med, color=GREY, lw=0.5, ls="-.")
stage_axis(ax, PS)
ax.set_ylabel("log$_{10}$ precursor\nquantity")
from matplotlib.lines import Line2D  # noqa: E402
hd = [Line2D([], [], color="#555555", ls="-", marker="o", ms=2, lw=0.8, label="OLE16a RVPGSEQLEQAR (3+)"),
      Line2D([], [], color="#555555", ls=":", marker="s", ms=2, lw=0.8, label="OLE16a VPGSEQLEQAR (2+)"),
      Line2D([], [], color="#555555", ls="--", marker="^", ms=2, lw=0.8, label="OLE16b RPPGFEQLEQAR (3+)"),
      Line2D([], [], color=GREY, ls="-.", lw=0.6, label="Median precursor (per run)")]
ax.legend(handles=hd, frameon=False, fontsize=4.8, loc="upper left", handlelength=2.4, borderaxespad=0.1, labelspacing=0.25)
ax.set_ylim(5.2, 11.6)
clean(ax)
letter(ax, "f", dx=-0.1)

# ---- g: LD-coat families ----
g = pd.read_csv(SD / "SF5g_LD_coat_families.tsv", sep="\t")
gv = "Summed directLFQ abundance (mean of 3 replicates)"
PH = ["185d", "12h", "24h", "36h", "48h", "60h", "72h"]
ax = fig.add_subplot(bot[1, 2])
for fam, mk in [("Oleosin", "o"), ("LDAP (REF/SRPP)", "s"), ("Caleosin", "^")]:
    for m, c in [("FL", FL), ("TN", TN)]:
        d = g[(g.Family == fam) & (g.Material == m)].set_index("Stage").reindex(PH)
        ax.plot(range(len(PH)), np.log10(d[gv]), color=c, lw=0.8, marker=mk, ms=2.4,
                ls="-" if fam == "Oleosin" else ("--" if fam.startswith("LDAP") else ":"))
stage_axis(ax, PH)
ax.set_ylabel("log$_{10}$ summed\nprotein abundance")
hd = [Line2D([], [], color="#555555", ls="-", marker="o", ms=2, lw=0.8, label="Oleosin"),
      Line2D([], [], color="#555555", ls="--", marker="s", ms=2, lw=0.8, label="LDAP"),
      Line2D([], [], color="#555555", ls=":", marker="^", ms=2, lw=0.8, label="Caleosin")]
ax.legend(handles=hd, frameon=False, fontsize=5, loc="upper center", ncol=3, handlelength=2.0, columnspacing=0.8, borderaxespad=0.1)
ax.set_ylim(7.0, 10.6)
clean(ax)
letter(ax, "g", dx=-0.1)

# ---- h: FL allele origin ----
h = pd.read_csv(SD / "SF5h_OLE16a_FL_allele_origin.tsv", sep="\t")
h = h[h["Allele-informative fragments"] >= 50]
ax = fig.add_subplot(bot[2, 1])
xi = h["Stage index"] - 1
ax.bar(xi, h["FL-Hap1 fraction"] * 100, width=0.7, color=FL)
for xx, n in zip(xi, h["Allele-informative fragments"]):
    ax.text(xx, 101.5, f"{n:,}", rotation=90, fontsize=4.6, ha="center", va="bottom", color="#555555")
ax.axhline(50, color=GREY, lw=0.5, ls="--")
stage_axis(ax, ST)
ax.set_xlim(6.4, 18.6)
ax.set_ylim(0, 100)
ax.set_ylabel("FL-Hap1 fragments (%)")
clean(ax)
letter(ax, "h", dx=-0.1, dy=0.03)

# ---- i: snRNA clusters ----
r = pd.read_csv(SD / "SF5i_OLE16a_snRNA_clusters.tsv", sep="\t")
fl = r[r.Library == "FL_185"].set_index("Cluster")
tn = r[r.Library == "TN_185"].set_index("Cluster")
cls = [c for c in fl.index if fl.loc[c, "Nuclei"] >= 30]
cls = sorted(cls, key=lambda c: -fl.loc[c, "% nuclei with >=1 OLE16a UMI"])
ax = fig.add_subplot(bot[2, 2])
xx = np.arange(len(cls))
ax.bar(xx - 0.2, [fl.loc[c, "% nuclei with >=1 OLE16a UMI"] for c in cls], 0.4, color=FL, label="FL 185 d (188 of 6,836)")
ax.bar(xx + 0.2, [tn.loc[c, "% nuclei with >=1 OLE16a UMI"] if c in tn.index else 0 for c in cls], 0.4,
       color=TN, label="TN 185 d (0 of 12,501)")
ax.set_xticks(xx, cls, rotation=90)
ax.set_ylabel("Nuclei with OLE16a\ntranscripts (%)")
ax.legend(frameon=False, fontsize=5.5, loc="upper right")
clean(ax)
letter(ax, "i", dx=-0.1, dy=0.03)

fig.canvas.draw()
for a_, s_, dx, dy in LETTERS:
    bb = a_.get_position()
    fig.text(bb.x0 + dx, bb.y1 + dy, s_, fontsize=8, fontweight="bold", va="bottom", ha="left")
fig.savefig(OUT / "Supplementary_Fig_05.pdf", facecolor="white")
fig.savefig(OUT / "Supplementary_Fig_05.png", dpi=600, facecolor="white")
im = Image.open(OUT / "Supplementary_Fig_05.png")
im.resize((2400, round(im.size[1] * 2400 / im.size[0])), Image.LANCZOS).save(OUT / "docx_2400/Supplementary_Fig_05.png")
print("saved", im.size, "clusters", cls)
