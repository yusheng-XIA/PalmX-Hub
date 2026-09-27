#!/usr/bin/env python3
"""Vector panels for the Figure 2 re-layout candidate (fix/fig2_ole16/).

Outputs (all vector PDF, Arial, TrueType-embedded):
  panels/c_axes.pdf        new panel c = former 2d (four composite axes), redrawn taller from Source Data Fig2d_axis_scores
  panels/row_ole16.pdf     OLE16a row, panels f-h (no promoter placeholder)
  panels/row_ole16_ph.pdf  same row with a placeholder panel i for the promoter comparison
  panels/sf4_metab.pdf     former 2c (Delta consensus z heatmap) redrawn for Supplementary Fig. 4
Data: deliver/Source_Data/Source_Data.xlsx (Fig2c_metabolites, Fig2d_axis_scores) and
      fix/idea_A/integrate/sd/SF5{d,e,g,h}_*.tsv (identical to the SF5d-i Source Data sheets).
"""
import io
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, BoundaryNorm
import numpy as np
import pandas as pd
from Bio import Phylo

HERE = Path(__file__).resolve().parent
S = HERE.parents[1]
sys.path.insert(0, str(S / "fix/beautify/common"))
import palA  # noqa: E402

FL, TN, EOL = palA.FL, palA.TN, palA.EOL
GREY, DARK, STRIP, SHADE = "#9A9A9A", "#222222", "#E6E7EA", "#EDEEF0"
PT = 1 / 72
plt.rcParams.update({
    "font.family": "Arial", "font.size": 6.5, "pdf.fonttype": 42, "svg.fonttype": "none",
    "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6,
    "xtick.major.size": 2.2, "ytick.major.size": 2.2, "axes.labelsize": 6.5,
    "xtick.labelsize": 6, "ytick.labelsize": 6, "legend.fontsize": 6, "axes.edgecolor": DARK,
    "xtick.color": DARK, "ytick.color": DARK, "text.color": DARK, "axes.labelcolor": DARK,
    "mathtext.fontset": "custom", "mathtext.rm": "Arial", "mathtext.it": "Arial:italic", "mathtext.bf": "Arial:bold",
})
EO = r"$\it{E.\ oleifera}$"
ST = ["0d", "15d", "35d", "50d", "65d", "80d", "95d", "110d", "125d", "140d", "155d", "170d", "185d",
      "12h", "24h", "36h", "48h", "60h", "72h"]
SDX = S / "deliver/Source_Data/Source_Data.xlsx"
SD5 = S / "fix/idea_A/integrate/sd"
OUT = HERE / "panels"
OUT.mkdir(exist_ok=True)


def clean(ax):
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)


def letter(fig, x_pt, y_pt, s, H_pt):
    fig.text(x_pt / fig.get_figwidth() / 72, 1 - y_pt / H_pt, s, fontsize=8, fontweight="bold",
             va="top", ha="left")


# ---------------------------------------------------------------- panel c (former 2d)
def panel_c():
    W, H = 247.0, 190.0
    fig = plt.figure(figsize=(W * PT, H * PT))
    d = pd.read_excel(SDX, sheet_name="Fig2d_axis_scores")
    order = [("P02", "Storage lipids"), ("P01", "Oleic balance"),
             ("P03", "Hydrolytic deterioration"), ("P04", "Oxidative deterioration")]
    gs = fig.add_gridspec(2, 2, left=0.13, right=0.985, top=0.80, bottom=0.10, hspace=0.62, wspace=0.28)
    for k, (aid, name) in enumerate(order):
        ax = fig.add_subplot(gs[k // 2, k % 2])
        ax.axvspan(12.5, 18.5, color=SHADE, lw=0, zorder=0)
        ax.axvline(12.5, color="#8C8C8C", ls=(0, (3, 2)), lw=0.6, zorder=1)
        ax.axhline(0, color="#BDBDBD", lw=0.5, zorder=1)
        for g, col in (("FL", FL), ("TN", TN)):
            x = d[(d.axis_id == aid) & (d.genotype == g)].set_index("stage").reindex(ST)
            xs = np.arange(19)
            ax.fill_between(xs, x["mean"] - x.se, x["mean"] + x.se, color=col, alpha=0.2, lw=0, zorder=2)
            ax.plot(xs, x["mean"], color=col, lw=0.9, marker="o", ms=1.8, zorder=3, label=g)
        ax.set_xlim(-0.6, 18.6)
        ax.set_xticks([0, 5, 12, 18], ["0d", "80d", "185d", "72h"])
        clean(ax)
        ax.add_patch(plt.Rectangle((0, 1.03), 1, 0.17, transform=ax.transAxes, fc=STRIP, ec="none", clip_on=False))
        ax.text(0.5, 1.115, name, transform=ax.transAxes, ha="center", va="center", fontsize=6.5, fontweight="bold")
        ax.yaxis.set_major_locator(matplotlib.ticker.MaxNLocator(3))
        if k % 2 == 0:
            ax.set_ylabel("Axis score" if k == 0 else "", labelpad=1)
    fig.text(0.02, 0.47, "Axis score", rotation=90, va="center", ha="left", fontsize=6.5)
    for a in fig.axes:
        a.set_ylabel("")
    fig.text(0.13, 0.905, "(mean ± s.e.m.)", fontsize=6, va="center")
    h = [plt.Line2D([], [], color=FL, marker="o", ms=1.8, lw=0.9), plt.Line2D([], [], color=TN, marker="o", ms=1.8, lw=0.9)]
    fig.legend(h, ["FL", "TN"], loc="center right", bbox_to_anchor=(0.985, 0.905), ncol=2, frameon=False,
               handlelength=1.6, columnspacing=1.0, borderaxespad=0)
    letter(fig, 0.6, 0.3, "c", H)
    fig.savefig(OUT / "c_axes.pdf")
    plt.close(fig)


# ---------------------------------------------------------------- OLE16a row
TIPS = {  # tree tip -> (display label, style)
    "OLE16a FL-Hap1 (chr11A)": ("OLE16a FL-Hap1", "a"),
    "OLE16a FL-Hap2 (chr11B)": ("OLE16a FL-Hap2", "a"),
    "OLE16a TN-Hap1": ("OLE16a TN", "a"),
    "Date palm LOC103704664": ("Date palm LOC103704664", "o"),
    "OLE16b FL-Hap2 (chr04B)": ("OLE16b FL-Hap2", "b"),
    "OLE16b TN-Hap1": ("OLE16b TN", "b"),
    "Date palm LOC103700664": ("Date palm LOC103700664", "o"),
    "Rice OLE16": ("Rice OLE16", "o"), "Maize OLE16": ("Maize OLE16", "o"), "Brome OLE16": ("Brome OLE16", "o"),
    "Arabidopsis OLE1": ("Arabidopsis OLE1", "o"), "Sesame OLE-L": ("Sesame OLE-L", "o"),
    "Rice OLE18": ("Rice OLE18", "o"), "Maize OLE18": ("Maize OLE18", "o"), "Sesame OLE-H1": ("Sesame OLE-H1", "o"),
}


def pruned_tree():
    t = Phylo.read(str(S / "fix/idea_A/integrate/tree_pruned.nwk"), "newick")
    for tip in list(t.get_terminals()):
        if tip.name not in TIPS:
            t.prune(tip)
    return t


def draw_tree(ax):
    t = pruned_tree()
    tips = t.get_terminals()
    y = {tip: i for i, tip in enumerate(tips)}
    depth = t.depths()
    if not max(depth.values()):
        depth = t.depths(unit_branch_lengths=True)

    def ypos(c):
        if c.is_terminal():
            return y[c]
        return np.mean([ypos(ch) for ch in c.clades])

    def draw(c):
        yc = ypos(c)
        for ch in c.clades:
            ych = ypos(ch)
            ax.plot([depth[c], depth[c]], [yc, ych], color=DARK, lw=0.55, solid_capstyle="butt")
            ax.plot([depth[c], depth[ch]], [ych, ych], color=DARK, lw=0.55, solid_capstyle="butt")
            draw(ch)
        if not c.is_terminal() and c.confidence is not None and c is not t.root:
            if c.confidence >= 95:
                ax.plot(depth[c], yc, "o", ms=2.0, color=DARK, zorder=5)
            elif c.confidence >= 70:
                ax.plot(depth[c], yc, "o", ms=2.0, mfc="white", mec=DARK, mew=0.5, zorder=5)

    draw(t.root)
    xmax = max(depth.values())
    for tip in tips:
        lab, sty = TIPS[tip.name]
        col = {"a": FL, "b": "#4D4D4D", "o": "#6B6B6B"}[sty]
        ax.text(depth[tip] + 0.02 * xmax, y[tip], lab, va="center", ha="left", fontsize=5.2, color=col,
                fontweight="bold" if sty in "ab" else "normal",
                fontstyle="italic" if False else "normal")
    # brackets: L-oleosin (OLE16/OLE1) vs H-oleosin (OLE18/OLE-H1)
    names = [tip.name for tip in tips]
    L = [i for i, n in enumerate(names) if not ("OLE18" in n or "OLE-H1" in n)]
    Hh = [i for i, n in enumerate(names) if ("OLE18" in n or "OLE-H1" in n)]
    import matplotlib.transforms as mt
    tr = mt.blended_transform_factory(ax.figure.transFigure, ax.transData)
    xb = 141 / 510
    for idx, lab in ((L, "L-oleosins"), (Hh, "H-oleosins")):
        ax.plot([xb, xb], [min(idx) - 0.3, max(idx) + 0.3], color=GREY, lw=0.6, clip_on=False, transform=tr)
        ax.text(xb + 3 / 510, np.mean(idx), lab, va="center", ha="left", fontsize=5.2, color="#555555",
                clip_on=False, transform=tr)
    ax.set_ylim(len(tips) - 0.4, -0.6)
    ax.set_xlim(-0.02 * xmax, xmax * 1.02)
    ax.axis("off")
    # scale bar
    ax.plot([0, 0.5], [len(tips) - 0.1] * 2, color=DARK, lw=0.6, clip_on=False)
    ax.text(0.25, len(tips) + 0.1, "0.5", ha="center", va="top", fontsize=5.2, clip_on=False)


def draw_rna(ax, axb):
    r = pd.read_csv(SD5 / "SF5e_OLE16_RNA.tsv", sep="\t")
    r = r[r.Locus == "OLE16a"]
    val = "Normalized count (DESeq2 size factors, 114 libraries)"
    r["l"] = np.log10(r[val] + 1)
    xs = np.arange(19)
    for g, col in (("TN", TN), ("FL", FL)):
        q = r[r.Material == g].groupby("Stage")["l"].agg(["mean", "sem"]).reindex(ST)
        ax.fill_between(xs, q["mean"] - q["sem"], q["mean"] + q["sem"], color=col, alpha=0.2, lw=0)
        ax.plot(xs, q["mean"], color=col, lw=0.9, marker="o", ms=1.8, label=g)
    for a in (ax, axb):
        a.axvspan(12.5, 18.6, color=SHADE, lw=0, zorder=0)
        a.axvline(12.5, color="#8C8C8C", ls=(0, (3, 2)), lw=0.6)
        a.set_xlim(-0.6, 18.6)
        clean(a)
    ax.set_ylabel("OLE16a RNA\nlog$_{10}$(count + 1)", labelpad=1)
    ax.set_ylim(-0.2, 4.6)
    ax.set_yticks([0, 2, 4])
    ax.set_xticks(xs, [])
    ax.legend(loc="upper left", frameon=False, handlelength=1.4, borderaxespad=0.1, ncol=2, columnspacing=0.8)
    a = pd.read_csv(SD5 / "SF5h_OLE16a_FL_allele_origin.tsv", sep="\t").set_index("Stage")
    for i, s in enumerate(ST):
        if s in a.index:
            f1 = a.loc[s, "FL-Hap1 fraction"] * 100
            axb.bar(i, f1, width=0.72, color=EOL, lw=0)
            axb.bar(i, 100 - f1, bottom=f1, width=0.72, color="#D9D9D9", lw=0)
    axb.set_ylim(0, 100)
    axb.set_yticks([0, 50, 100])
    axb.set_ylabel("FL allele\norigin (%)", labelpad=1)
    axb.set_xticks(xs, [s if s in ("0d", "80d", "140d", "185d", "72h") else "" for s in ST], rotation=0)
    axb.tick_params(axis="x", length=1.5)
    h = [plt.Rectangle((0, 0), 1, 1, color=EOL), plt.Rectangle((0, 0), 1, 1, color="#D9D9D9")]
    axb.legend(h, ["FL-Hap1 (" + EO + "-derived)", "FL-Hap2"], loc="lower left", bbox_to_anchor=(-0.01, -1.02),
               ncol=2, frameon=False, handlelength=0.9, handleheight=0.8, columnspacing=0.8, borderaxespad=0,
               fontsize=5.5)


def draw_ratio(ax):
    g = pd.read_csv(SD5 / "SF5g_LD_coat_families.tsv", sep="\t")
    v = "Summed directLFQ abundance (mean of 3 replicates)"
    pv = g.pivot_table(index=["Material", "Stage"], columns="Family", values=v)
    stg = ["185d", "12h", "24h", "36h", "48h", "60h", "72h"]
    xs = np.arange(len(stg))
    for m, col in (("FL", FL), ("TN", TN)):
        q = pv.loc[m].reindex(stg)
        rr = q["Oleosin"] / q["LDAP (REF/SRPP)"]
        ax.plot(xs, rr, color=col, lw=0.9, marker="o", ms=2.0, label=m)
        ax.annotate({"FL": "0.46", "TN": "0.005"}[m], (xs[0], rr.iloc[0]), xytext=(0, 4 if m == "FL" else -4),
                    textcoords="offset points", ha="center", va="bottom" if m == "FL" else "top", fontsize=5.5, color=col)
    ax.set_yscale("log")
    ax.yaxis.set_minor_locator(matplotlib.ticker.NullLocator())
    ax.set_ylim(1e-3, 20)
    ax.set_xlim(-0.6, len(stg) - 0.4)
    ax.set_xticks(xs, stg, rotation=45, ha="right")
    ax.set_ylabel("Oleosin : LDAP\n(protein abundance)", labelpad=1)
    clean(ax)
    ax.legend(loc="upper right", frameon=False, handlelength=1.4, borderaxespad=0.1, ncol=2, columnspacing=0.8)


def row(placeholder):
    W, H = 510.0, 108.0
    fig = plt.figure(figsize=(W * PT, H * PT))
    if placeholder:
        tx, gx, hx, ix = (0, 170), (207, 338), (374, 446), (458, 510)
    else:
        tx, gx, hx = (0, 170), (205, 368), (408, 506)
    fx = lambda a, b: (a / W, b / W)  # noqa: E731
    at = fig.add_axes([tx[0] / W + 0.004, 0.07, (tx[1] - tx[0]) / W * 0.52, 0.86])
    draw_tree(at)
    x0, x1 = fx(*gx)
    ag = fig.add_axes([x0, 0.47, x1 - x0, 0.46])
    agb = fig.add_axes([x0, 0.26, x1 - x0, 0.17])
    draw_rna(ag, agb)
    x0, x1 = fx(*hx)
    ah = fig.add_axes([x0, 0.30, x1 - x0, 0.63])
    draw_ratio(ah)
    letter(fig, 0, 1, "f", H)
    letter(fig, gx[0] - 30, 1, "g", H)
    letter(fig, hx[0] - 30, 1, "h", H)
    if placeholder:
        x0, x1 = fx(*ix)
        ai = fig.add_axes([x0, 0.12, x1 - x0, 0.80])
        ai.set_xticks([]); ai.set_yticks([])
        for s in ai.spines.values():
            s.set_linestyle((0, (2, 2))); s.set_color(GREY); s.set_linewidth(0.6)
        ai.text(0.5, 0.5, "OLE16a upstream\ncomparison\n(FL-Hap1, FL-Hap2,\nTN alleles)\n\nplaceholder", ha="center",
                va="center", fontsize=5.5, color=GREY, transform=ai.transAxes)
        letter(fig, ix[0] - 9, 1, "i", H)
    fig.savefig(OUT / ("row_ole16_ph.pdf" if placeholder else "row_ole16.pdf"))
    plt.close(fig)


# ---------------------------------------------------------------- SF4 candidate panel (former 2c)
def sf4_panel():
    d = pd.read_excel(SDX, sheet_name="Fig2c_metabolites")
    names = list(dict.fromkeys(d.display_name))
    piv = d.pivot_table(index=["display_name", "stage"], columns="genotype", values="mean_consensus_z")
    delta = (piv["FL"] - piv["TN"]).unstack("stage").reindex(index=names, columns=ST)
    W, H = 360.0, 118.0
    fig = plt.figure(figsize=(W * PT, H * PT))
    ax = fig.add_axes([0.20, 0.25, 0.68, 0.58])
    cols = [palA.FL_TN[4], palA.FL_TN[3], "#F7F7F7", palA.FL_TN[1], palA.FL_TN[0]]
    cmap = LinearSegmentedColormap.from_list("d", cols[::1][::-1][::-1])
    cmap = LinearSegmentedColormap.from_list("d", [palA.FL_TN[4], palA.FL_TN[3], "#F7F7F7", palA.FL_TN[1], palA.FL_TN[0]])
    cmap.set_bad("#BDC0C6")
    m = np.ma.masked_invalid(delta.values)
    im = ax.pcolormesh(np.arange(20) - 0.5, np.arange(len(names) + 1) - 0.5, m, cmap=cmap, vmin=-2, vmax=2,
                       edgecolors="white", linewidth=0.4)
    ax.set_ylim(len(names) - 0.5, -0.5)
    ax.set_yticks(range(len(names)), [("$\\it{p}$-Coumaric acid" if n.lower().startswith("p-coumaric") else n) for n in names])
    ax.set_xticks(range(19), ST, rotation=90)
    ax.axvline(12.5, color="#666666", ls=(0, (3, 2)), lw=0.6)
    ax.tick_params(length=0)
    for s in ax.spines.values():
        s.set_visible(False)
    ax.text(6, -0.9, "Development (d)", ha="center", va="bottom", fontsize=6)
    ax.text(15.5, -0.9, "Post-harvest (h)", ha="center", va="bottom", fontsize=6)
    fig.text(0.20, 0.97, "Δ consensus z (FL − TN)", fontsize=6.5, va="top")
    cax = fig.add_axes([0.91, 0.45, 0.018, 0.36])
    cb = fig.colorbar(im, cax=cax, ticks=[-2, 0, 2])
    cb.ax.set_yticklabels(["≤−2", "0", "≥2"])
    cb.outline.set_visible(False)
    cb.ax.tick_params(length=0, labelsize=6)
    fig.text(0.905, 0.33, "■", color="#BDC0C6", fontsize=7, va="center")
    fig.text(0.93, 0.33, "NA", fontsize=6, va="center")
    letter(fig, 0, 1, "d", H)
    fig.savefig(OUT / "sf4_metab.pdf")
    fig.savefig(OUT / "sf4_metab.png", dpi=300)
    plt.close(fig)


if __name__ == "__main__":
    panel_c()
    row(False)
    row(True)
    sf4_panel()
    print("panels written to", OUT)
