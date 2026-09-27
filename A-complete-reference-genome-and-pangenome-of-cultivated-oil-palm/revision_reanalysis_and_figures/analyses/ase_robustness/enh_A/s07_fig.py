#!/usr/bin/env python3
"""Candidate Supplementary Figure: genome-wide vs trait-module ASE direction in FL (TN control) and DNA mapping-bias control."""
import numpy as np, pandas as pd, matplotlib as mpl
mpl.use("Agg"); import matplotlib.pyplot as plt
from pathlib import Path
O = Path(__file__).parent/"out"; MM = 1/25.4
mpl.rcParams.update({"font.family":"sans-serif","font.sans-serif":["Arial"],"font.size":6,"axes.labelsize":6,
    "xtick.labelsize":5.5,"ytick.labelsize":5.5,"legend.fontsize":5.5,"axes.linewidth":.6,"xtick.major.width":.6,
    "ytick.major.width":.6,"pdf.fonttype":42,"ps.fonttype":42})
MODS = ["Oil biosynthesis & storage","TAG assembly & oil body","De-novo / saturated FA","Unsaturated FA","Lipid oxidation / antioxidant","Shell / cell wall / lignin"]
AB = ["OBS","TOF","DSF","UFA","LOD","SCL"]; MC = ["#E07B39","#C9A227","#8C6BB1","#D6604D","#4D9221","#7F7F7F"]
g = pd.read_csv(O/"gene_window_ASE_unified.tsv.gz", sep="\t")
mem = pd.read_csv(O/"module_membership.tsv", sep="\t")
p1 = pd.read_csv(O/"part1_genome_vs_module_unified.tsv", sep="\t")
p2 = pd.read_csv(O/"part2_corrected_genome_vs_module.tsv", sep="\t")
G = pd.read_csv(O/"gene_DNA_ratios.tsv.gz", sep="\t", index_col=0)
k3 = pd.read_csv(O/"part3_key_FA_genes.tsv", sep="\t")
fig = plt.figure(figsize=(180*MM, 128*MM))
def letter(x, y, s): fig.text(x, y, s, fontsize=8, fontweight="bold", va="top")

# a: per-gene all-stage median log2(A/B), robust genes, genome-wide vs modules, FL and TN
for j, an in enumerate(["FL", "TN"]):
    ax = fig.add_axes([0.06 + j*0.25, 0.60, 0.20, 0.33])
    x = g[(g.analysis == an) & (g.stage_group == "All stages") & g.robust_any]
    sets = [("Genome", x.med_robust.values, "#BBBBBB")]
    for m, a, c in zip(MODS, AB, MC):
        gs = set(mem[(mem.analysis == an) & (mem.trait_module == m)].gene_id)
        sets.append((a, x[x.gene_id.isin(gs)].med_robust.values, c))
    for i, (lab, v, c) in enumerate(sets):
        vv = np.clip(v, -6, 6)
        if i == 0:
            parts = ax.violinplot(vv, positions=[i], widths=.8, showextrema=False)
            for b in parts["bodies"]: b.set_facecolor(c); b.set_alpha(.8); b.set_edgecolor("none")
        else:
            jit = (np.random.default_rng(i).random(len(vv)) - .5) * .5
            ax.scatter(i + jit, vv, s=2.5, color=c, lw=0, alpha=.85)
        ax.hlines(np.median(v), i - .35, i + .35, color="k", lw=.9)
        ax.text(i, 6.6, f"{100*(v<0).mean():.0f}%", ha="center", va="bottom", fontsize=5)
    ax.axhline(0, color="#555", lw=.5, ls=(0, (2, 2)))
    ax.set_xticks(range(len(sets)), [s[0] for s in sets], rotation=45, ha="right")
    ax.set_ylim(-6.5, 7.6); ax.set_yticks([-6, -3, 0, 3, 6])
    ax.set_ylabel("Per-gene median log2(A/B)\n(robust ASE stages)" if j == 0 else "")
    ax.spines[["top", "right"]].set_visible(False)
    lab = "FL (A = FL-Hap2, B = FL-Hap1)" if an == "FL" else "TN (A = TK-like, B = NS-like)"
    ax.set_title(lab, fontsize=6, pad=10)
fig.text(0.31, 0.985, "Numbers above each set: % B-biased genes", fontsize=5, color="#444", ha="center", va="top")
letter(0.005, 0.99, "a")

# b: FL % B-biased genes per window, genome-wide vs modules, before and after DNA correction
ax = fig.add_axes([0.60, 0.60, 0.38, 0.33]); letter(0.555, 0.99, "b")
WIN = ["Days 0–65", "Days 80–140", "Days 155–185", "Hours 12–72", "All stages"]
WL = ["0–65 d", "80–140 d", "155–185 d", "12–72 h", "All"]
vu = "RNA uncorrected (genes with Short-read reciprocal DNA)"; vc = "RNA - DNA (Short-read reciprocal)"
for wi, w in enumerate(WIN):
    for off, v, filled in [(-0.17, vu, True), (0.17, vc, False)]:
        d = p2[(p2.version == v) & (p2.window == w)]
        gw = d[d.set == "Genome-wide"].pct_B.iloc[0]
        ax.hlines(gw, wi + off - .14, wi + off + .14, color="k", lw=1.2, zorder=3)
        for m, c in zip(MODS, MC):
            r = d[d.set == m]
            if len(r):
                ax.scatter(wi + off, r.pct_B.iloc[0], s=10, marker="o", facecolor=c if filled else "white", edgecolor=c, lw=.7, zorder=4)
ax.axhline(50, color="#555", lw=.5, ls=(0, (2, 2)))
ax.set_xticks(range(5), WL); ax.set_ylabel("B-biased genes among robust ASE genes (%)"); ax.set_ylim(35, 85)
ax.spines[["top", "right"]].set_visible(False)
from matplotlib.lines import Line2D
h = [Line2D([0], [0], color="k", lw=1.2, label="Genome-wide")] + [Line2D([0], [0], marker="o", ls="", mfc=c, mec=c, ms=3.5, label=a) for a, c in zip(AB, MC)]
h += [Line2D([0], [0], marker="o", ls="", mfc="#999", mec="#999", ms=3.5, label="RNA"), Line2D([0], [0], marker="o", ls="", mfc="white", mec="#999", ms=3.5, label="RNA − DNA")]
ax.legend(handles=h, ncol=5, frameon=False, loc="upper center", bbox_to_anchor=(.5, 1.13), columnspacing=.8, handletextpad=.2)

# c: DNA allele ratios at the same diagnostic sites (gene level)
ax = fig.add_axes([0.06, 0.09, 0.36, 0.36]); letter(0.005, 0.47, "c")
bins = np.linspace(-2, 2, 81)
rna = g[(g.analysis == "FL") & (g.stage_group == "All stages")].med_all
for v, lab, c, ls in [(G.l2_srA, "DNA, reads on FL-Hap2 (A) reference", "#2166AC", "-"),
                      (G.l2_srB, "DNA, reads on FL-Hap1 (B) reference", "#B2182B", "-"),
                      (G.l2_sr_sym, "DNA, reciprocal mean", "k", "-"),
                      (rna, "RNA, per-gene median (all eligible stages)", "#E08214", "--")]:
    v = v.dropna(); hh, _ = np.histogram(np.clip(v, -2, 2), bins=bins); hh = hh / hh.sum() * 100
    ax.step(bins[:-1], hh, where="post", color=c, lw=.8, ls=ls, label=f"{lab} (median {np.median(v):+.2f})")
ax.axvline(0, color="#555", lw=.5, ls=(0, (2, 2)))
ax.set_xlabel("Gene-level log2(A/B) (values beyond ±2 pooled at the edges)"); ax.set_ylabel("Genes (%)"); ax.set_xlim(-2, 2); ax.set_ylim(0, 16)
ax.legend(frameon=False, loc="upper left", fontsize=5, bbox_to_anchor=(0, 1.02))
ax.spines[["top", "right"]].set_visible(False)

# d: key FA genes, DNA-corrected log2(A/B)
ax = fig.add_axes([0.58, 0.09, 0.40, 0.40]); letter(0.33, 0.50, "d")
k = k3[k3.corr_median_log2AB.notna() & (k3.FL_TPM >= 5)].copy()
order = ["SAD", "FAD2", "FAD3", "FATA/B", "KAS I/II", "KASIII", "LPAT", "DGAT", "PDAT"]
k["o"] = k.enzyme.map({e: i for i, e in enumerate(order)}); k = k.sort_values(["o", "FL_TPM"], ascending=[True, False]).reset_index(drop=True)
col = {"E. oleifera (B) allele higher": "#B2182B", "African (A) allele higher": "#2166AC"}
for i, r in k.iterrows():
    c = col.get(r.verdict_corrected, "#999999")
    v = np.clip(r.corr_median_log2AB, -4, 4)
    ax.hlines(i, 0, v, color=c, lw=1.4)
    ax.scatter(v, i, s=10, color=c, zorder=3)
    ax.scatter(r.DNA_sr_recip, i, marker="|", s=18, color="k", lw=.7, zorder=4)
    if abs(r.corr_median_log2AB) > 4: ax.text(v + (0.15 if v > 0 else -0.15), i, f"{r.corr_median_log2AB:.0f}", va="center", ha="left" if v > 0 else "right", fontsize=5)
ax.set_yticks(range(len(k)), [f"{e} {gid.replace('evm.TU.', '')} ({tpm:.0f} TPM)" for e, gid, tpm in zip(k.enzyme, k.gene_A_FLHap2, k.FL_TPM)], fontsize=5)
ax.invert_yaxis(); ax.axvline(0, color="#555", lw=.5)
ax.set_xlim(-4.6, 4.6); ax.set_xlabel("DNA-corrected median log2(A/B) across stages")
ax.spines[["top", "right"]].set_visible(False)
h = [Line2D([0], [0], color=col["E. oleifera (B) allele higher"], lw=1.4, label="B (E. oleifera-derived) higher"),
     Line2D([0], [0], color=col["African (A) allele higher"], lw=1.4, label="A (African-derived) higher"),
     Line2D([0], [0], color="#999999", lw=1.4, label="Mixed / stage-dependent"),
     Line2D([0], [0], marker="|", ls="", color="k", ms=5, label="DNA log2(A/B)")]
ax.legend(handles=h, frameon=False, loc="lower left", fontsize=5)
fig.savefig(O.parent/"SuppFig_ASE_genome_vs_module_DNAcontrol.pdf")
fig.savefig(O.parent/"SuppFig_ASE_genome_vs_module_DNAcontrol.png", dpi=600)
