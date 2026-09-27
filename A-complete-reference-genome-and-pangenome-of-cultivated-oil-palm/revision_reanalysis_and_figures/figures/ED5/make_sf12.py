#!/usr/bin/env python3
"""Supplementary Fig. 12: genome-wide allelic bias in FL and TN versus trait modules, and the FL DNA allele-ratio
control.

Inputs : sd/SF12a..SF12d (sf12_data.py), work/sf12_stats.json.
Outputs: SF12/Supplementary_Fig_12.pdf (vector, Arial, 180 mm wide), SF12/Supplementary_Fig_12.png (600 dpi).
Colours: palette scheme A (FL #23897D, TN #DE6B63; A-biased red #C8504A, B-biased blue #3F7FB5, direction as in
Fig. 3k).
"""
import json
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

H = Path(__file__).resolve().parent; SD = H / "sd"; OUT = H / "SF12"; OUT.mkdir(exist_ok=True)
A = pd.read_csv(SD / "SF12a_genome_vs_module.tsv", sep="\t")
Gv = pd.read_csv(SD / "SF12a_gene_values.tsv", sep="\t")
B = pd.read_csv(SD / "SF12b_window_DNAcorrected.tsv", sep="\t")
C = pd.read_csv(SD / "SF12c_gene_DNA_RNA.tsv", sep="\t")
D = pd.read_csv(SD / "SF12d_key_FA_genes.tsv", sep="\t")
ST = json.load(open(H / "work/sf12_stats.json"))

plt.rcParams.update({"font.family": "Arial", "font.size": 6, "axes.linewidth": 0.5, "xtick.major.width": 0.5,
                     "ytick.major.width": 0.5, "xtick.major.size": 2, "ytick.major.size": 2, "axes.labelsize": 6,
                     "xtick.labelsize": 5.5, "ytick.labelsize": 5.5, "pdf.fonttype": 42, "ps.fonttype": 42,
                     "axes.spines.top": False, "axes.spines.right": False, "legend.frameon": False,
                     "legend.fontsize": 5.5, "axes.unicode_minus": True, "mathtext.fontset": "custom",
                     "mathtext.rm": "Arial", "mathtext.it": "Arial:italic", "mathtext.bf": "Arial:bold"})
MAT = {"FL": "#23897D", "TN": "#DE6B63"}   # palette scheme A (fix/beautify/common/palA.py)
CA, CB = "#C8504A", "#3F7FB5"   # A-biased red, B-biased blue (scheme A; Fig. 3k direction)
MODS = ["OBS", "TOF", "DSF", "UFA", "LOD", "SCL"]
MK = dict(zip(MODS, ["o", "s", "^", "v", "D", "P"]))
WINS = ["0–65 d", "80–140 d", "155–185 d", "12–72 h", "All stages"]
LOG2 = r"log$_2$(A/B)"


def sci(p, digits=0):
    m, e = f"{p:.{digits}e}".split("e")
    return rf"{m} $\times$ 10$^{{{int(e)}}}$"


fig = plt.figure(figsize=(180 / 25.4, 168 / 25.4))

# ================= a: per-gene median log2(A/B), genome-wide vs modules (all stages)
ALAB = {"FL": ("FL-Hap2", r"FL-Hap1 ($\it{E.\,oleifera}$-derived)"), "TN": ("TK-like", "NS-like")}
for j, an in enumerate(["FL", "TN"]):
    x0 = 0.075 + j * 0.49
    ax = fig.add_axes([x0, 0.695, 0.36, 0.225]); axb = fig.add_axes([x0, 0.655, 0.36, 0.03])
    gv = Gv[Gv.Hybrid == an]
    sets = [("Genome", gv.Median_log2AB_robust_stages.values)] + \
           [(m, gv[gv[m]].Median_log2AB_robust_stages.values) for m in MODS]
    rows = A[(A.Hybrid == an) & (A.Window == "All stages")]
    for i, (lab, v) in enumerate(sets):
        vv = np.clip(v, -6, 6)
        if i == 0:
            vp = ax.violinplot(vv, positions=[i], widths=0.8, showextrema=False, bw_method=0.08)
            for b in vp["bodies"]:
                b.set_facecolor(MAT[an]); b.set_alpha(0.35); b.set_edgecolor(MAT[an]); b.set_linewidth(0.4)
        else:
            jit = (np.random.default_rng(i).random(len(vv)) - 0.5) * 0.5
            ax.scatter(i + jit, vv, s=2.2, color=MAT[an], lw=0, alpha=0.8, zorder=2)
        ax.hlines(np.median(v), i - 0.33, i + 0.33, color="k", lw=0.8, zorder=3)
        # stacked A/B bar
        r = rows[rows.Set == "Genome-wide"].iloc[0] if i == 0 else rows[rows.Set.str.startswith(lab + " (")].iloc[0]
        pb, pa = r.Pct_B_biased, r.Pct_A_biased
        axb.bar(i, pb, 0.8, color=CB, lw=0); axb.bar(i, pa, 0.8, bottom=pb, color=CA, lw=0)
        axb.text(i, pb / 2, f"{pb:.0f}", ha="center", va="center", fontsize=5, color="white")
    ax.axhline(0, color="0.4", lw=0.4, ls=(0, (2, 2)), zorder=1)
    ax.set_xlim(-0.6, 6.6); ax.set_ylim(-6.6, 6.6); ax.set_yticks([-6, -3, 0, 3, 6])
    ax.set_xticks([]); ax.spines["bottom"].set_visible(False)
    ax.set_ylabel(f"Per-gene median {LOG2}")
    ax.text(6.75, 6.3, f"A-biased\n({ALAB[an][0]})", color=CA, fontsize=5, va="top", ha="left", linespacing=1.1)
    ax.text(6.75, -6.3, f"B-biased\n({ALAB[an][1]})" if an == "TN" else "B-biased\n(FL-Hap1)", color=CB,
            fontsize=5, va="bottom", ha="left", linespacing=1.1)
    axb.set_xlim(-0.6, 6.6); axb.set_ylim(0, 100); axb.set_yticks([]); axb.spines["left"].set_visible(False)
    axb.set_xticks(range(7)); axb.tick_params(axis="x", length=0, pad=1.5)
    ns = [int(rows[rows.Set == "Genome-wide"].Robust_ASE_genes.item())] + \
         [int(rows[rows.Set.str.startswith(m + " (")].Robust_ASE_genes.item()) for m in MODS]
    axb.set_xticklabels([f"{l}\n{n:,}" for l, n in zip(["Genome"] + MODS, ns)], fontsize=5.5)
    axb.text(-0.75, 50, "% B", ha="right", va="center", fontsize=5, color=CB)
    st = ST[an]
    if an == "FL":
        ttl = (f"FL: {st['pct_B']:.1f}% of {st['n_robust']:,} genes favour FL-Hap1, $\\it{{P}}$ = {sci(st['binom_P'])}")
    else:
        ttl = (f"TN: {st['pct_A']:.1f}% of {st['n_robust']:,} genes favour the TK-like allele, "
               f"$\\it{{P}}$ = {sci(st['binom_P'])}")
    ax.text(0, 1.105, ttl, transform=ax.transAxes, fontsize=6, color=MAT[an], va="bottom", ha="left")
    ax.text(0, 1.035, f"Six modules vs genome-wide, 5 windows: permutation $\\it{{P}}$ ≥ "
            f"{np.floor(st['perm_two_min'] * 100) / 100:.2f} (two-sided)", transform=ax.transAxes, fontsize=5,
            va="bottom", ha="left", color="0.3")

# ================= b: FL % B-biased genes per window, RNA vs RNA - DNA
ax = fig.add_axes([0.075, 0.355, 0.40, 0.215])
for wi, w in enumerate(WINS):
    for off, v in [(-0.18, "RNA"), (0.18, "RNA − DNA")]:
        d = B[(B.Version == v) & (B.Window == w)]
        gw = d[d.Set == "Genome-wide"].Pct_B_biased.item()
        ax.hlines(gw, wi + off - 0.15, wi + off + 0.15, color=MAT["FL"], lw=1.6, zorder=3)
        for m in MODS:
            r = d[d.Set.str.startswith(m + " (")]
            if len(r):
                ax.scatter(wi + off, r.Pct_B_biased.item(), s=9, marker=MK[m], facecolor="0.3" if v == "RNA" else "white",
                           edgecolor="0.3", lw=0.5, zorder=4)
ax.axhline(50, color="0.4", lw=0.4, ls=(0, (2, 2)))
ax.axvspan(3.5, 4.5, color="0.94", zorder=0, lw=0)
ax.set_xticks(range(5)); ax.set_xticklabels(WINS); ax.set_xlim(-0.5, 4.5)
ax.set_ylim(33, 90); ax.set_yticks([40, 50, 60, 70, 80, 90])
ax.set_ylabel("FL genes favouring FL-Hap1 (%)")
cb = ST["b_corr"]
ax.text(0.01, 0.095, f"Genome-wide, all stages: RNA {ST['b_rna']['pct_B_all']:.1f}%, RNA − DNA {cb['pct_B_all']:.1f}%",
        transform=ax.transAxes, ha="left", va="bottom", fontsize=5, color=MAT["FL"])
ax.text(0.01, 0.025, "Filled, RNA; open, RNA − DNA (DNA-corrected)", transform=ax.transAxes, ha="left", va="bottom",
        fontsize=5, color="0.3")
h = [Line2D([], [], color=MAT["FL"], lw=1.6)] + [Line2D([], [], marker=MK[m], ls="", mfc="0.3", mec="0.3", ms=3) for m in MODS]
ax.legend(h, ["Genome-wide"] + MODS, ncol=7, loc="lower left", bbox_to_anchor=(-0.01, 1.0),
          handletextpad=0.15, columnspacing=0.7, handlelength=1.1, borderaxespad=0.2)

# ================= c: gene-level DNA vs RNA log2(A/B)
ax = fig.add_axes([0.075, 0.065, 0.40, 0.2])
bins = np.linspace(-3, 3, 61); ctr = (bins[:-1] + bins[1:]) / 2
for col, lab, c, ls, lw in [("DNA_log2AB_FL-Hap2_ref", "DNA, FL-Hap2 reference", "0.55", (0, (1, 1)), 0.8),
                            ("DNA_log2AB_FL-Hap1_ref", "DNA, FL-Hap1 reference", "0.55", (0, (3, 1.5)), 0.8),
                            ("DNA_log2AB_reciprocal_mean", "DNA, reciprocal mean", "k", "-", 0.9),
                            ("RNA_median_log2AB_all_eligible_stages", "RNA, eligible stages", MAT["FL"], "-", 0.9)]:
    v = C[col].dropna().values
    hh, _ = np.histogram(np.clip(v, -2.999, 2.999), bins=bins); hh = hh / len(v) * 100
    ax.step(bins, np.r_[hh, hh[-1]], where="post", color=c, lw=lw, ls=ls,
            label=f"{lab} (median {np.median(v):.2f}; $\\it{{n}}$ = {len(v):,})".replace("-", "−"))
ax.plot([0, 0], [0, 20.5], color="0.4", lw=0.4, ls=(0, (2, 2)))
ax.set_xlim(-3, 3); ax.set_ylim(0, 30)
ax.set_xticks([-3, -2, -1, 0, 1, 2, 3]); ax.set_xticklabels(["≤−3", "−2", "−1", "0", "1", "2", "≥3"])
ax.set_xlabel(f"Gene-level {LOG2} in FL"); ax.set_ylabel("Genes (%)")
ax.legend(loc="upper left", fontsize=5, handlelength=1.8, bbox_to_anchor=(0.0, 1.02), labelspacing=0.3)
ax.text(0.99, 0.55, f"DNA reciprocal mean:\n|{LOG2}| < 0.5\nfor {100 * ST['c']['frac_lt05']:.1f}% of genes",
        transform=ax.transAxes, ha="right", va="top", fontsize=6, linespacing=1.3)

# ================= d: key fatty-acid genes, DNA-corrected RNA log2(A/B) with DNA log2(A/B)
k = D[D.Plotted.astype(str) == "True"].copy()
order = ["SAD", "FAD2", "FAD3", "FATA/B", "KAS I/II", "KASIII", "LPAT", "DGAT", "PDAT"]
k["o"] = k.Enzyme.map({e: i for i, e in enumerate(order)})
k = k.sort_values(["o", "FL_mean_TPM"], ascending=[True, False]).reset_index(drop=True)
ax = fig.add_axes([0.745, 0.065, 0.235, 0.505])
col = {"FL-Hap1 (B) higher": CB, "FL-Hap2 (A) higher": CA}
for i, r in k.iterrows():
    c = col.get(r.Call_after_DNA_correction, "0.6")
    v = float(np.clip(r.DNA_corrected_median_log2AB, -3, 3))
    ax.hlines(i, 0, v, color=c, lw=1.2)
    ax.scatter(v, i, s=9, color=c, zorder=3, lw=0)
    ax.scatter(r["DNA_log2AB_reciprocal_mean"], i, marker="|", s=16, color="k", lw=0.7, zorder=4)
    if abs(r.DNA_corrected_median_log2AB) > 3:
        ax.text(v + 0.15, i, f"{r.DNA_corrected_median_log2AB:.1f}", va="center", ha="left", fontsize=5)
lab = [f"{e}  {g.replace('evm.TU.', '')}  ({t:,.0f})" for e, g, t in zip(k.Enzyme, k["Gene_ID_FL-Hap2"], k.FL_mean_TPM)]
ax.set_yticks(range(len(k))); ax.set_yticklabels(lab, fontsize=5)
ax.set_ylim(len(k) - 0.4, -0.6)
ax.axvline(0, color="0.4", lw=0.4)
ax.set_xlim(-3.3, 3.9); ax.set_xticks([-3, -2, -1, 0, 1, 2, 3])
ax.set_xlabel(f"DNA-corrected median {LOG2}")
ax.text(-0.02, 1.005, "Gene (FL TPM)", transform=ax.transAxes, ha="right", va="bottom", fontsize=5)
h = [Line2D([], [], color=CA, lw=1.2, marker="o", ms=2.5), Line2D([], [], color=CB, lw=1.2, marker="o", ms=2.5),
     Line2D([], [], color="0.6", lw=1.2, marker="o", ms=2.5), Line2D([], [], marker="|", ls="", color="k", ms=4, mew=0.7)]
ax.legend(h, ["FL-Hap2 higher", "FL-Hap1 higher", "Stage-dependent", f"DNA {LOG2}"], loc="upper right", bbox_to_anchor=(1.03, 0.947),
          fontsize=6, handlelength=1.0, handletextpad=0.5, borderaxespad=0.1, labelspacing=0.3)

for t, (xp, yp) in {"a": (0.008, 0.985), "b": (0.008, 0.615), "c": (0.008, 0.29), "d": (0.535, 0.615)}.items():
    fig.text(xp, yp, t, fontsize=8, fontweight="bold", va="top", ha="left")
fig.savefig(OUT / "Supplementary_Fig_12.pdf")
fig.savefig(OUT / "Supplementary_Fig_12.png", dpi=600)
print("saved", OUT)
