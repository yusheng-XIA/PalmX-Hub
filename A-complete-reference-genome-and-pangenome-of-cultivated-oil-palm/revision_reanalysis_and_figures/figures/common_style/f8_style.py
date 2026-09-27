"""Shared style for the F8 redraws (Supplementary Figs. 4, 5, 9); mirrors plot_SF6.py."""
import pickle
from pathlib import Path

import matplotlib as mpl

HERE = Path(__file__).resolve().parent
TEAL, RED, DARK, GREY, PURPLE, BLUE, ORANGE = ("#2A9D8F", "#D95F5F", "#242A30", "#B9C0C7",
                                               "#756BB1", "#3C78A8", "#E69F00")
MM = 1 / 25.4

mpl.rcParams.update({
    "font.family": "Arial", "font.size": 7, "axes.labelsize": 7, "axes.titlesize": 7,
    "xtick.labelsize": 6, "ytick.labelsize": 6, "legend.fontsize": 6,
    "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6,
    "xtick.major.size": 2.5, "ytick.major.size": 2.5, "axes.unicode_minus": True,
    "mathtext.fontset": "custom", "mathtext.rm": "Arial", "mathtext.it": "Arial:italic",
    "mathtext.bf": "Arial:bold", "pdf.fonttype": 42, "svg.fonttype": "none",
})


def clean(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_axisbelow(True)


def letter(fig, ax, s, dx=-0.055, dy=0.012):
    bb = ax.get_position()
    fig.text(bb.x0 + dx, bb.y1 + dy, s, fontsize=9, fontweight="bold", va="bottom", ha="left")


def sheet(name):
    d = pickle.load(open(HERE / "data/f8_raw.pkl", "rb"))
    v = d[name].copy()
    v.columns = v.iloc[0]
    v = v.iloc[1:].reset_index(drop=True)
    return v.infer_objects()


def stage_label(s):
    s = str(s)
    return s[:-1] + " " + s[-1] if s[-1] in "dh" else s


def sci(p, digits=2):
    """P value as 'm × 10^{e}' mathtext."""
    m, e = f"{p:.{digits}e}".split("e")
    return rf"{m} \times 10^{{{int(e)}}}"
