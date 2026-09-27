#!/usr/bin/env python3
"""Fig. 3k candidate redraw: unique-gene metric, drawn at 1:1 page scale.

Geometry, fonts, colours and line widths are measured from the visible 3k
panel in deliver/Main_Figures_revised/Figure3.pdf (page units, pt), so the new
panel can be placed back with show_pdf_page at an identity scale.
Style origin: final_ms/01_figure3/redraw_Fig3b_Fig3k_Fig3l_final.py::render_k
(coolwarm, TwoSlopeNorm +/-2.5, FL circles #E75F5F edge, TN squares #2F68A2 edge,
bands #F5F7F7 / #FFF7E5), after the 0.8189 figure-level scaling.

usage: render_fig3k.py {old|new} out.pdf
  old = legacy source_06A values + legacy labels/legend (placement validation only)
  new = Figure3k_plotdata_unique_genes.tsv (candidate)
"""
import sys
from pathlib import Path
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from matplotlib.patches import Rectangle
import pandas as pd

HERE = Path(__file__).resolve().parent
SCR = HERE.parent.parent
PHASES = ["Days 0–65", "Days 80–140", "Days 155–185", "Hours 12–72"]
MODULES = ["Oil biosynthesis & storage", "TAG assembly & oil body", "De-novo / saturated FA",
           "Unsaturated FA", "Lipid oxidation / antioxidant", "Shell / cell wall / lignin"]
ABBR = ["OBS", "TOF", "DSF", "UFA", "LOD", "SCL"]

# ---- page-coordinate frame of the replaced area (pt) ----
X0, X1, Y0, Y1 = 282.0, 424.8, 381.0, 500.0
SC = 0.8189                     # figure-level scale applied to the published panel
FS = 5.04                       # 6.15 pt * 0.8189
C_TICK, C_TXT = '#333333', '#222222'
COL0, DCOL = 303.73, 15.795     # x of data column 0, column pitch
ROW0, DROW = 408.61, 11.717     # y of data row 0, row pitch
AX_Y0, AX_Y1 = 400.41, 474.23   # band top / bottom
HEAD_X = [311.25, 342.70, 374.50, 405.30]; HEAD_BASE = 389.6
TICK_BASE = 398.2
ROW_RIGHT = 293.1; ROW_BASE = [410.4, 422.1, 433.8, 445.5, 457.3, 469.0]
CB = (361.84, 410.82, 478.98, 481.88); CB_TICK_LEN = 1.31; CB_LAB_BASE = 487.6
CB_TITLE = (386.35, 493.2)
LEG_TITLE = (298.2, 482.0); LEG_MK_X = [297.03, 312.93, 328.83]; LEG_MK_Y = 486.4
LEG_LAB_X = [301.1, 317.0, 332.9]; LEG_LAB_BASE = 488.2

mpl.rcParams.update({'font.family': 'sans-serif', 'font.sans-serif': ['Arial'],
                     'pdf.fonttype': 42, 'ps.fonttype': 42, 'axes.unicode_minus': True})


def load(mode):
    if mode == "old":
        d = pd.read_csv(SCR / "trace/Fig3B/src/source_06A_FL_TN_trait_complement.tsv", sep="\t")
        return d.rename(columns={"robust_ASE_percentage": "pct", "robust_median_log2_ratio": "med"})
    d = pd.read_csv(HERE / "Figure3k_plotdata_unique_genes.tsv", sep="\t")
    return d.rename(columns={"robust_ASE_unique_gene_percentage": "pct",
                             "median_of_within_gene_robust_log2_ratios": "med"})


def size_old(p):   # legacy: s = 7 + 0.21*pct (pt^2 at panel scale)
    return (7 + .21 * p) * SC ** 2


def size_new(p):   # unique-gene range 74.2-100 %: linear-in-area, spread over the new range
    return (8 + .9 * (p - 70)) * SC ** 2


def render(mode, path):
    new = mode == "new"
    d = load(mode)
    size = size_new if new else size_old
    W, H = X1 - X0, Y1 - Y0
    fig = plt.figure(figsize=(W / 72, H / 72))
    ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(X0, X1); ax.set_ylim(Y1, Y0); ax.axis('off')
    norm = TwoSlopeNorm(vmin=-2.5, vcenter=0, vmax=2.5)
    cmap = mpl.colormaps['coolwarm']
    # time-window bands
    for pi in range(4):
        c = '#F5F7F7' if pi < 3 else '#FFF7E5'
        xa = COL0 + (pi * 2 - .5) * DCOL
        ax.add_patch(Rectangle((xa, AX_Y0), 2 * DCOL, AX_Y1 - AX_Y0, fc=c, ec=c, lw=SC, zorder=0))
    # symbols
    for yi, mod in enumerate(MODULES):
        for pi, phase in enumerate(PHASES):
            for ai, an in enumerate(['FL', 'TN']):
                q = d[(d.trait_module == mod) & (d.stage_group == phase) & (d.analysis == an)]
                assert len(q) == 1, (mod, phase, an)
                r = q.iloc[0]
                ax.scatter(COL0 + (pi * 2 + ai) * DCOL, ROW0 + yi * DROW, s=size(float(r.pct)),
                           c=[float(r.med)], cmap=cmap, norm=norm, marker='o' if an == 'FL' else 's',
                           edgecolors='#E75F5F' if an == 'FL' else '#2F68A2', linewidths=.55 * SC, zorder=3)
    t = dict(fontsize=FS, va='baseline')
    heads = ['0–65 d', '80–140 d', '155–185 d', '12–72 h'] if new else ['Early', 'Middle', 'Late', 'Postharvest']
    for x, s in zip(HEAD_X, heads):
        ax.text(x, HEAD_BASE, s, ha='center', color=C_TXT, **t)
    for j in range(8):
        ax.text(COL0 + j * DCOL, TICK_BASE, ['FL', 'TN'][j % 2], ha='center', color=C_TICK, **t)
    for yb, s in zip(ROW_BASE, ABBR):
        ax.text(ROW_RIGHT, yb, s, ha='right', color=C_TICK, **t)
    # colour bar
    cx0, cx1, cy0, cy1 = CB
    cax = fig.add_axes([(cx0 - X0) / W, (Y1 - cy1) / H, (cx1 - cx0) / W, (cy1 - cy0) / H])
    cb = fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap), cax=cax, orientation='horizontal')
    cb.set_ticks([-2.5, 0, 2.5]); cb.outline.set_linewidth(.55 * SC); cb.outline.set_edgecolor(C_TICK)
    cax.tick_params(length=CB_TICK_LEN, width=.5 * SC, color=C_TICK, labelbottom=False)
    for v, s in zip([-2.5, 0, 2.5], ['−2.5', '0.0', '2.5']):
        ax.text(cx0 + (v + 2.5) / 5 * (cx1 - cx0), CB_LAB_BASE, s, ha='center', color=C_TICK, **t)
    ax.text(CB_TITLE[0], CB_TITLE[1], 'Median log2(A/B)', ha='center', color=C_TXT, **t)
    # size legend
    ax.text(LEG_TITLE[0], LEG_TITLE[1], 'Robust ASE genes (%)' if new else 'Robust ASE (%)', ha='left', color=C_TXT, **t)
    keys = [80, 90, 100] if new else [50, 70, 90]
    for xm, xl, v in zip(LEG_MK_X, LEG_LAB_X, keys):
        ax.scatter(xm, LEG_MK_Y, s=size(v), marker='o', c='#DDE2E5', edgecolors='#777777', linewidths=.45 * SC, zorder=3)
        ax.text(xl, LEG_LAB_BASE, str(v), ha='left', color=C_TXT, **t)
    fig.savefig(path, transparent=True, metadata={'Title': 'Fig3k_unique_genes', 'Creator': 'render_fig3k.py'})
    plt.close(fig)


if __name__ == '__main__':
    render(sys.argv[1], sys.argv[2])
