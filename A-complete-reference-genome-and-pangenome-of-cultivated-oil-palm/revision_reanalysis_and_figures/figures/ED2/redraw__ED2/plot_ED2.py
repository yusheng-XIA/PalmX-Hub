#!/usr/bin/env python3
"""Extended Data Fig. 2 redraw: multi-haplotype microsynteny at four annotation-based PAV candidates.

Inputs are the tables written by the original pipeline
(40_plot_lipid_priority_multihaplotype_microsynteny.py; 03_V3/.../08_gene_pav/EG11_FL_Africa/microsynteny/):
Tracks.tsv, Genes.tsv (with roles), Adjacent_Track_Links.tsv, Target_Presence.tsv.
Geometry (x mapping, Bezier ribbons, inferred PAV position, colours, alpha) follows the original
draw_figure(); only sizes, fonts and haplotype labels are changed.
"""
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import PathPatch, Rectangle
from matplotlib.path import Path as MplPath

HERE = Path(__file__).resolve().parent
MM = 1 / 25.4
mpl.rcParams.update({"font.family": "Arial", "font.size": 6, "pdf.fonttype": 42, "axes.unicode_minus": True,
                     "mathtext.fontset": "custom", "mathtext.rm": "Arial", "mathtext.it": "Arial:italic"})

TARGET_COLOR, CNV_COLOR, SYN_COLOR, OTHER_COLOR, TRACK_COLOR = "#E68445", "#4C78A8", "#B8B8B8", "#E6E6E6", "#333333"
LABEL = {"FL_Africa_hap2": "FL-Hap2", "EG11": "EG11", "dura_hap1": "TK-Hap1", "dura_hap2": "TK-Hap2",
         "pisifera_hap1": "NS-Hap1", "pisifera_hap2": "NS-Hap2", "nrly_hap1": "Nigerian-Hap1",
         "nrly_hap2": "Nigerian-Hap2", "BK_hap1": "TN-Hap1", "BK_hap2": "TN-Hap2", "FL_American_hap1": "FL-Hap1"}
PANELS = [("a", "chr03B.814", "PDAT1"), ("b", "chr01B.3982", "Polyketide cyclase-like"),
          ("c", "chr04B.1103", "GDSL lipase"), ("d", "chr05B.150", "HACD")]
X0, X1 = 0.20, 0.76


def load(gid):
    d = HERE / "src" / gid
    tr = pd.read_csv(d / "Tracks.tsv", sep="\t")
    ge = pd.read_csv(d / "Genes.tsv", sep="\t")
    li = pd.read_csv(d / "Adjacent_Track_Links.tsv", sep="\t")
    tp = pd.read_csv(d / "Target_Presence.tsv", sep="\t")
    return tr, ge, li, tp


def xpos(v, t):
    return X0 + (v - t.Region_Start) / max(1, t.Region_End - t.Region_Start) * (X1 - X0)


def ribbon(x1, y1, x2, y2, hw=0.0045):
    dl = (y2 - y1) * 0.45
    v = [(x1 - hw, y1), (x1 - hw, y1 + dl), (x2 - hw, y2 - dl), (x2 - hw, y2), (x2 + hw, y2),
         (x2 + hw, y2 - dl), (x1 + hw, y1 + dl), (x1 + hw, y1), (x1 - hw, y1)]
    c = [MplPath.MOVETO] + [MplPath.CURVE4] * 3 + [MplPath.LINETO] + [MplPath.CURVE4] * 3 + [MplPath.CLOSEPOLY]
    return MplPath(v, c)


def chrom_label(t):
    import re
    m = re.search(r"chr(1[0-6]|0?[1-9])", str(t.Source_Seqid)) or re.search(r"chr(1[0-6]|0?[1-9])", str(t.Seqid))
    return f"chr{int(m.group(1)):02d}" if m else str(t.Seqid)


def inferred_x(t, genes_t, fl_genes, target_og):
    idx = next(i for i, g in enumerate(fl_genes.itertuples()) if g.Orthogroup == target_og)
    up = [g.Orthogroup for g in list(fl_genes.itertuples())[:idx][::-1]]
    dn = [g.Orthogroup for g in list(fl_genes.itertuples())[idx + 1:]]
    ug = dg = None
    for og in up:
        c = genes_t[genes_t.Orthogroup == og]
        if len(c):
            ug = c.loc[c.End.idxmax()]; break
    for og in dn:
        c = genes_t[genes_t.Orthogroup == og]
        if len(c):
            dg = c.loc[c.Start.idxmin()]; break
    cx = lambda g: xpos((g.Start + g.End) / 2, t)
    if ug is not None and dg is not None:
        return (cx(ug) + cx(dg)) / 2
    if ug is not None:
        return min(X1 - 0.015, cx(ug) + 0.03)
    if dg is not None:
        return max(X0 + 0.015, cx(dg) - 0.03)
    return (X0 + X1) / 2


def draw(ax, gid, name, summary):
    tr, ge, li, tp = load(gid)
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")
    target_og = tp.Orthogroup.dropna().iloc[0]
    samples = list(tr.Sample)
    n = len(samples)
    ytop, ybot = 0.93, 0.03
    Y = {s: ytop - i * (ytop - ybot) / (n - 1) for i, s in enumerate(samples)}
    T = {r.Sample: r for r in tr.itertuples()}
    G = {s: ge[ge.Sample == s] for s in samples}
    gene_pos = {(r.Sample, r.Gene_ID): xpos((r.Start + r.End) / 2, T[r.Sample]) for r in ge.itertuples()}
    for r in li.itertuples():
        col, al = {"Target_PAV": (TARGET_COLOR, 0.31), "CNV": (CNV_COLOR, 0.30)}.get(r.Link_Type, (SYN_COLOR, 0.20))
        ax.add_patch(PathPatch(ribbon(gene_pos[(r.Upper_Sample, r.Upper_Gene)], Y[r.Upper_Sample],
                                      gene_pos[(r.Lower_Sample, r.Lower_Gene)], Y[r.Lower_Sample]),
                               fc=col, ec="none", alpha=al, zorder=1))
    fl_genes = G["FL_Africa_hap2"].sort_values("Start")
    fl_genes = fl_genes[fl_genes.Orthogroup.notna()]
    absent = []
    for s in samples:
        t, y = T[s], Y[s]
        ax.add_line(Line2D([X0 + 0.005, X1], [y, y], color=TRACK_COLOR, lw=0.5, zorder=2))
        ax.text(X0 - 0.015, y, LABEL[s], ha="right", va="center", fontsize=5.6,
                fontweight="bold" if s in ("FL_Africa_hap2", "EG11") else "normal")
        ax.text(X1 + 0.012, y, f"{chrom_label(t)} {t.Region_Start / 1e6:.2f}–{t.Region_End / 1e6:.2f} Mb",
                ha="left", va="center", fontsize=4.9, color="#666666")
        for g in G[s].itertuples():
            a, b = max(X0, xpos(g.Start, t)), min(X1, xpos(g.End, t))
            if b <= X0 or a >= X1:
                continue
            w = max(0.0035, b - a)
            col = {"Target_PAV": TARGET_COLOR, "CNV": CNV_COLOR, "Syntenic_ortholog": SYN_COLOR}.get(g.Role, OTHER_COLOR)
            ax.add_patch(Rectangle(((a + b) / 2 - w / 2, y - 0.014), w, 0.028, fc=col, ec="none", zorder=3))
        if not (G[s].Orthogroup == target_og).any():
            ptp = getattr(t, "Projected_Target_Position", None)
            x = xpos(float(ptp), t) if ptp is not None and pd.notna(ptp) and str(ptp) != "NA" else inferred_x(t, G[s], fl_genes, target_og)
            ax.add_patch(Rectangle((x - 0.007, y - 0.02), 0.014, 0.04, fill=False, ec=TARGET_COLOR, lw=0.7,
                                   ls=(0, (2.5, 1.5)), zorder=4))
            star = "" if getattr(t, "Anchor_Context", "bilateral") == "bilateral" else "*"
            ax.text(x, y + 0.026, "PAV" + star, ha="center", va="bottom", fontsize=4.9, color=TARGET_COLOR)
            absent.append(LABEL[s] + star)
    # E. guineensis / E. oleifera-derived separator above FL-Hap1 (last track)
    ysep = (Y[samples[-2]] + Y[samples[-1]]) / 2
    ax.add_line(Line2D([X0, X1], [ysep, ysep], color="#666666", lw=0.5, ls=(0, (4, 3))))
    ax.text(X1 + 0.012, ysep + 0.016, r"$\mathit{E.\ guineensis}$", fontsize=4.9, color="#777777", va="center")
    ax.text(X1 + 0.012, ysep - 0.016, r"$\mathit{E.\ oleifera}$-derived", fontsize=4.9, color="#777777", va="center")
    summary.append((gid, name, target_og, ";".join(absent)))


fig = plt.figure(figsize=(180 * MM, 168 * MM))
summary = []
pos = {"a": (0.0, 0.49), "b": (0.5, 0.49), "c": (0.0, 0.0), "d": (0.5, 0.0)}
for letter, gid, name in PANELS:
    x0, y0 = pos[letter]
    ax = fig.add_axes([x0 + 0.005, y0 + 0.012, 0.49, 0.40])
    draw(ax, gid, name, summary)
    fig.text(x0 + 0.008, y0 + 0.418, letter, fontsize=9, fontweight="bold", va="bottom")
    fig.text(x0 + 0.035, y0 + 0.419, f"{name} (evm.TU.{gid})", fontsize=6.5, va="bottom",
             fontstyle="normal")
handles = [Rectangle((0, 0), 1, 1, fc=TARGET_COLOR, ec="none", label="Focal PAV candidate"),
           Rectangle((0, 0), 1, 1, fc=CNV_COLOR, ec="none", label="Local copy-number-variable orthogroup"),
           Rectangle((0, 0), 1, 1, fc=SYN_COLOR, ec="none", label="Syntenic orthologue"),
           Rectangle((0, 0), 1, 1, fc=OTHER_COLOR, ec="none", label="Other local gene"),
           Rectangle((0, 0), 1, 1, fill=False, ec=TARGET_COLOR, ls=(0, (2.5, 1.5)), lw=0.7,
                     label="Inferred position of unannotated candidate")]
fig.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, 0.995), ncol=5, frameon=False, fontsize=5.6,
           handlelength=1.0, handleheight=0.8, columnspacing=1.2)
for ext in ("pdf", "png"):
    fig.savefig(HERE / f"Extended_Data_Fig_02.{ext}", dpi=600, facecolor="white")
pd.DataFrame(summary, columns=["focal_gene", "name", "target_orthogroup", "tracks_without_annotated_target"]).to_csv(
    HERE / "ED2_panel_summary.tsv", sep="\t", index=False)
print(summary)
