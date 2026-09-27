#!/usr/bin/env python3
"""Faithfully redraw Figure 3a, 3b and 3d from reviewed original sources.

Data sources
------------
Fig. 3b: exact counts in plot_allele_classification.R.
Fig. 3d: complete 4 x 4 FL/TN transition table in source_E_FL_TN_alluvial.tsv.
Fig. 3a: schematic (no quantitative values), matched to Figure3-A4_1.5fold.pdf.

Physical output sizes
---------------------
Fig3a  58 x 22 mm
Fig3b  58 x 30 mm
Fig3ab 58 x 52 mm (a over b, no scaling)
Fig3d  38 x 42 mm
All displayed text is 6--8 pt.
"""
from __future__ import annotations

import hashlib
import os
import tempfile
from pathlib import Path

import fitz
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch, PathPatch, Rectangle
from matplotlib.path import Path as MplPath
import numpy as np
import pandas as pd

OUT = Path(__file__).resolve().parent
MM = 1 / 25.4
D_SOURCE = Path(
    "${ANALYSIS_DIR}/"
    "22_answer_reviews/00_ms/03_V3/03_figure3/04_FL_TN_ASE_reference_panels/"
    "source_E_FL_TN_alluvial.tsv"
)
B_SOURCE = Path(
    "${ANALYSIS_DIR}/"
    "21_MS/02_result/01_figure/plot_allele_classification.R"
)
REFERENCE = Path(
    "${ANALYSIS_DIR}/"
    "22_answer_reviews/00_ms/05_MS/new_revision/Final_figures/"
    "Figure3-A4_1.5fold.pdf"
)

# Colors sampled/matched to the supplied previous Figure 3.
BIAL = "#76B2D4"
SAME = "#D58774"
HAPS = "#F1D992"
ABSENT_EDGE = "#A6A6A6"
TRACK1 = "#6F9DC0"
TRACK2 = "#C99169"

ORDER = ["NoDiff", "HapDom", "Sub", "NoASE"]
D_COLORS = {"NoDiff": "#F8766D", "HapDom": "#FFF19A", "Sub": "#82CBBB", "NoASE": "#58A6CF"}


def set_style() -> None:
    mpl.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Liberation Sans", "DejaVu Sans"],
        "font.size": 6.0,
        "axes.labelsize": 6.0,
        "xtick.labelsize": 6.0,
        "ytick.labelsize": 6.0,
        "legend.fontsize": 6.0,
        "legend.title_fontsize": 6.0,
        "axes.linewidth": 0.65,
        "xtick.major.width": 0.65,
        "ytick.major.width": 0.65,
        "xtick.major.size": 2.3,
        "ytick.major.size": 2.3,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "savefig.facecolor": "white",
    })


def panel_letter(fig, value: str) -> None:
    fig.text(0.012, 0.985, value, fontsize=8, fontweight="bold", ha="left", va="top")


def save_pdf(fig, path: Path) -> None:
    fig.savefig(
        path,
        facecolor="white",
        metadata={
            "Title": path.stem,
            "Subject": "Faithful layout redraw from reviewed original Figure 3 sources",
            "Creator": Path(__file__).name,
        },
    )
    plt.close(fig)


def draw_gene(ax, x: float, y: float, kind: str, track_color: str) -> None:
    colors = {"biallelic": BIAL, "same": SAME, "hap": HAPS, "absent": "white"}
    ax.add_patch(Rectangle(
        (x - 0.043, y - 0.065), 0.086, 0.13,
        facecolor=colors[kind],
        edgecolor=ABSENT_EDGE if kind == "absent" else track_color,
        linewidth=0.5,
        linestyle=(0, (4, 3)) if kind == "absent" else "solid",
        zorder=3,
    ))


def render_a(path: Path) -> None:
    fig = plt.figure(figsize=(58 * MM, 22 * MM), facecolor="white")
    panel_letter(fig, "a")
    ax = fig.add_axes([0.075, 0.08, 0.90, 0.84])
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")

    y1, y2 = 0.78, 0.19
    ax.text(0.005, y1, "Hap1", fontweight="bold", ha="left", va="center")
    ax.text(0.005, y2, "Hap2", fontweight="bold", ha="left", va="center")
    x0, x1 = 0.105, 0.985
    ax.add_patch(Rectangle((x0, y1 - 0.045), x1-x0, 0.09, fc="#EEF4F7", ec=TRACK1, lw=0.45, zorder=1))
    ax.add_patch(Rectangle((x0, y2 - 0.045), x1-x0, 0.09, fc="#FFF7EF", ec=TRACK2, lw=0.45, zorder=1))

    xs = [0.16, 0.30, 0.44, 0.58, 0.72, 0.86]
    top = ["same", "biallelic", "biallelic", "hap", "biallelic", "absent"]
    bottom = ["same", "biallelic", "biallelic", "absent", "biallelic", "hap"]
    pair_cols = {"same": SAME, "biallelic": BIAL}
    for x, t, b in zip(xs, top, bottom):
        draw_gene(ax, x, y1, t, TRACK1)
        draw_gene(ax, x, y2, b, TRACK2)
        # Keep paired-gene connectors in the left schematic field only; the
        # right middle field is reserved for the legend, so no text crosses art.
        if x < 0.50 and t == b and t in pair_cols:
            ax.add_patch(Rectangle((x-0.031, y2+0.075), 0.062, y1-y2-0.15,
                                   fc=mpl.colors.to_rgba(pair_cols[t], 0.33), ec="none", zorder=0))

    handles = [
        Patch(fc=BIAL, ec="none", label="Biallelic"),
        Patch(fc=SAME, ec="none", label="Same CDS"),
        Patch(fc=HAPS, ec="none", label="Hap-specific"),
        Patch(fc="white", ec=ABSENT_EDGE, lw=.5, linestyle=(0,(4,3)), label="Gene absent"),
    ]
    # Opaque legend background prevents any connector from crossing text.
    leg = ax.legend(handles=handles, frameon=True, facecolor="white", edgecolor="none",
                    framealpha=1, ncol=2, loc="center", bbox_to_anchor=(0.68, 0.49),
                    columnspacing=0.8, handlelength=1.25, handletextpad=0.35,
                    borderpad=0.25, labelspacing=0.35)
    leg.set_zorder(10)
    save_pdf(fig, path)


def render_b(path: Path) -> None:
    # Exact counts copied from the reviewed source script; no inferred values.
    labels = ["FL_hap1", "FL_hap2", "TN_hap1", "TN_hap2"]
    counts = np.array([
        [27106, 264, 3214],
        [27106, 264, 7430],
        [18439, 7747, 7720],
        [18439, 7747, 8486],
    ], dtype=float)
    pct = counts / counts.sum(axis=1, keepdims=True) * 100.0

    fig = plt.figure(figsize=(58 * MM, 30 * MM), facecolor="white")
    panel_letter(fig, "b")
    ax = fig.add_axes([0.17, 0.22, 0.80, 0.69])
    x = np.arange(4)
    bottom = np.zeros(4)
    for vals, color in zip(pct.T, [BIAL, SAME, HAPS]):
        ax.bar(x, vals, bottom=bottom, width=0.74, color=color,
               edgecolor="white", linewidth=0.45)
        bottom += vals

    ax.axvline(1.5, color="#777777", linewidth=0.65, linestyle=(0, (1.5, 2.4)))
    ax.set_xlim(-0.55, 3.55); ax.set_ylim(0, 102)
    ax.set_yticks([0, 25, 50, 75, 100])
    ax.set_ylabel("Proportion of genes (%)", labelpad=2.0)
    ax.set_xticks(x, labels)
    ax.tick_params(axis="x", length=0, pad=2)
    ax.tick_params(axis="y", pad=1.5)
    ax.spines[["top", "right"]].set_visible(False)
    save_pdf(fig, path)


def intervals(values: pd.Series, gap: float) -> dict[str, tuple[float, float]]:
    usable = 1.0 - gap * (len(ORDER) - 1)
    out = {}
    top = 1.0
    for category in ORDER:
        height = usable * float(values[category]) / float(values.sum())
        out[category] = (top-height, top)
        top -= height + gap
    return out


def render_d(path: Path) -> None:
    data = pd.read_csv(D_SOURCE, sep="\t")
    total = int(data.orthogroups.sum())
    assert total == 7059
    assert set(zip(data.overall_class_FL, data.overall_class_TN)) == {(a,b) for a in ORDER for b in ORDER}

    left = data.groupby("overall_class_FL").orthogroups.sum().reindex(ORDER)
    right = data.groupby("overall_class_TN").orthogroups.sum().reindex(ORDER)
    gap = 0.012
    usable = 1.0 - gap * 3
    li, ri = intervals(left, gap), intervals(right, gap)
    lc = {k: li[k][0] for k in ORDER}
    rc = {k: ri[k][0] for k in ORDER}

    fig = plt.figure(figsize=(38 * MM, 42 * MM), facecolor="white")
    panel_letter(fig, "d")
    handles = [Patch(fc=D_COLORS[k], ec="#333333", lw=.35, label=k) for k in ORDER]
    fig.legend(handles=handles, frameon=False, ncol=2, loc="upper center",
               bbox_to_anchor=(0.57, 0.975), columnspacing=0.8,
               handlelength=1.15, handletextpad=0.35, labelspacing=0.35)

    ax = fig.add_axes([0.11, 0.18, 0.80, 0.62])
    ax.set_xlim(0, 1); ax.set_ylim(-0.075, 1.02); ax.axis("off")
    xL0, xL1, xR0, xR1 = 0.08, 0.20, 0.80, 0.92

    # Stable source-first ordering reproduces the accepted alluvial structure.
    order_rank = {v:i for i,v in enumerate(ORDER)}
    rows = data.assign(_l=data.overall_class_FL.map(order_rank), _r=data.overall_class_TN.map(order_rank)).sort_values(["_l", "_r"])
    for row in rows.itertuples():
        a, b, n = row.overall_class_FL, row.overall_class_TN, int(row.orthogroups)
        h = usable * n / total
        y0a, y1a = lc[a], lc[a] + h
        y0b, y1b = rc[b], rc[b] + h
        c = 0.23
        verts = [(xL1,y0a),(xL1+c,y0a),(xR0-c,y0b),(xR0,y0b),
                 (xR0,y1b),(xR0-c,y1b),(xL1+c,y1a),(xL1,y1a),(xL1,y0a)]
        codes = [MplPath.MOVETO,MplPath.CURVE4,MplPath.CURVE4,MplPath.CURVE4,
                 MplPath.LINETO,MplPath.CURVE4,MplPath.CURVE4,MplPath.CURVE4,MplPath.CLOSEPOLY]
        ax.add_patch(PathPatch(MplPath(verts,codes), fc=D_COLORS[a], ec="none", alpha=.40))
        lc[a] += h; rc[b] += h

    for category in ORDER:
        lo, hi = li[category]
        ax.add_patch(Rectangle((xL0,lo),xL1-xL0,hi-lo,fc=D_COLORS[category],ec="#333333",lw=.5,zorder=4))
        lo, hi = ri[category]
        ax.add_patch(Rectangle((xR0,lo),xR1-xR0,hi-lo,fc=D_COLORS[category],ec="#333333",lw=.5,zorder=4))
    ax.text((xL0+xL1)/2, -0.04, "FL", ha="center", va="top")
    ax.text((xR0+xR1)/2, -0.04, "TN", ha="center", va="top")
    fig.text(0.5, 0.045, "Shared ASE-eligible OGs, n = 7,059", ha="center", va="bottom")
    save_pdf(fig, path)


def combine_ab(a_path: Path, b_path: Path, out_path: Path) -> None:
    a = fitz.open(a_path); b = fitz.open(b_path)
    w = 58 * 72 / 25.4
    ha = 22 * 72 / 25.4
    hb = 30 * 72 / 25.4
    doc = fitz.open()
    page = doc.new_page(width=w, height=ha+hb)
    page.show_pdf_page(fitz.Rect(0, 0, w, ha), a, 0)
    page.show_pdf_page(fitz.Rect(0, ha, w, ha+hb), b, 0)
    doc.set_metadata({"title":"Fig3ab — faithful redraw", "author":"layout-only redraw", "subject":"Reviewed Figure 3 sources"})
    doc.save(out_path, garbage=4, deflate=True)
    doc.close(); a.close(); b.close()


def png_from_pdf(pdf: Path, png: Path) -> None:
    doc = fitz.open(pdf)
    pix = doc[0].get_pixmap(matrix=fitz.Matrix(600/72, 600/72), alpha=False)
    pix.save(png)
    doc.close()


def sha256(p: Path) -> str:
    h=hashlib.sha256()
    with p.open("rb") as f:
        for block in iter(lambda:f.read(1024*1024), b""): h.update(block)
    return h.hexdigest()


def main() -> None:
    set_style()
    for p in [D_SOURCE, B_SOURCE, REFERENCE]:
        if not p.is_file(): raise FileNotFoundError(p)
    with tempfile.TemporaryDirectory(prefix="fig3abd_redraw_", dir=OUT) as td:
        td=Path(td)
        render_a(td/"Fig3a.pdf")
        render_b(td/"Fig3b.pdf")
        combine_ab(td/"Fig3a.pdf", td/"Fig3b.pdf", td/"Fig3ab.pdf")
        render_d(td/"Fig3d.pdf")
        for stem in ["Fig3a","Fig3b","Fig3ab","Fig3d"]:
            png_from_pdf(td/f"{stem}.pdf", td/f"{stem}_600dpi.png")
        # Replace only the explicitly requested revised panels atomically.
        for p in sorted(td.glob("Fig3*")):
            os.replace(p, OUT/p.name)

    outputs=[OUT/f"{s}.{ext}" for s in ["Fig3a","Fig3b","Fig3ab","Fig3d"] for ext in ["pdf","png"]]
    # The PNG names include resolution suffix.
    outputs=[]
    for s in ["Fig3a","Fig3b","Fig3ab","Fig3d"]:
        outputs += [OUT/f"{s}.pdf", OUT/f"{s}_600dpi.png"]
    (OUT/"FIGURE3_ABD_REVISED_CHECKSUMS.sha256").write_text("\n".join(f"{sha256(p)}  {p.name}" for p in outputs)+"\n")
    (OUT/"FIGURE3_ABD_REVISED_SOURCES.tsv").write_text(
        "panel\trole\tsha256\tpath\n"
        f"ab\tdata_and_plot_source\t{sha256(B_SOURCE)}\t{B_SOURCE}\n"
        f"d\ttrue_4x4_transition_data\t{sha256(D_SOURCE)}\t{D_SOURCE}\n"
        f"ab,d\tvisual_reference\t{sha256(REFERENCE)}\t{REFERENCE}\n"
    )

if __name__ == "__main__":
    main()
