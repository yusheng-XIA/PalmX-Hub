#!/usr/bin/env python3
"""pan-BGC overview redraw — reproduces the design of the reference figure
`BGC.png` (Nature-review working copy) from the real pan-BGC data.

Data (read-only):
  ${ANALYSIS_DIR}/22_answer_reviews/00_ms/03_V3/04_figure4/
      panBGC_visualization_revised_diversity_20260708/tables_v20260912/
          st20_matrix.tsv                 52 families x 33 materials, copy number 0/1/2/3
          v20260912_family_order.tsv      family class / major type / prevalence / BGC count
          v20260912_material_order.tsv    column order + per-material BGC count
          st20_summary.tsv                per-material cluster counts

Design (from BGC.png):
  heatmap 52x33 (grey / blue / orange) + class strip + type strip (6 categories)
  + prevalence dot panel (n = 33) + gold BGC-copy bar panel with dashed 48 reference
  + dashed class separators, F-family row labels, 4 bottom legends.

Outputs: figures/Fig_panBGC_overview_v5.{pdf,svg,png}, tables/*.tsv, qc/*.

Note on fidelity: the reference figure is a design mock-up whose internal geometry is
NOT self-consistent (see qc/QA_REPORT.md). This script keeps the design but drives every
element from the real matrix and enforces strict row alignment.
"""
from __future__ import annotations

import hashlib
import json
import platform
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# journal typography: Arial (falls back to the metric-compatible Liberation Sans)
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Arial", "Liberation Sans", "Helvetica", "DejaVu Sans"]
from matplotlib.colors import ListedColormap, to_rgba
from matplotlib.patches import Rectangle

# --------------------------------------------------------------------------------------
# paths / constants
# --------------------------------------------------------------------------------------
ROOT = Path(__file__).resolve().parent.parent
DATA = Path(
    "${ANALYSIS_DIR}/22_answer_reviews/00_ms/03_V3/"
    "04_figure4/panBGC_visualization_revised_diversity_20260708/tables_v20260912"
)
FIG = ROOT / "figures"
TAB = ROOT / "tables"
QC = ROOT / "qc"

MM = 1.0 / 25.4
FIG_W_MM, FIG_H_MM = 103.0, 62.0           # double-column width per projects/PLOTTING.md
FIG_W, FIG_H = FIG_W_MM * MM, FIG_H_MM * MM

# V4 low-saturation light palette; materials use uniform neutral typography.
C_ABSENT, C_PRESENT, C_MULTI = "#F1F2EE", "#78A99B", "#CF907D"
BAR_GOLD, DOT_INK, GUIDE = "#829CB3", "#445B68", "#EBEEEE"
CLASS_COLORS = {"Core": "#91A9BE", "Soft-core": "#D9BE82",
                "Shell": "#94B5A4", "Unique": "#CE9C9E"}
TYPE_COLORS = {"Saccharide": "#A8B7CC", "Cyclopeptide": "#D6B08A", "Putative": "#AFBF9B",
               "Fatty acid": "#BBAACB", "Polyketide": "#BDAA93", "Others": "#C4C9C8"}
DASH_SEP = dict(color="#3A3A3A", linewidth=0.55, dashes=(2.4, 2.0))
SEP_FAMILIES = {"fatty_acid", "fatty_acid-polyketide"}   # -> "Fatty acid" (count 3)

# layout, in figure fractions (measured from the reference figure, 1520x1012 px)
L = dict(
    strip_class=(1.5/103, 1.3/103), strip_type=(8.2/103, 1.3/103),
    heat=(10/103, 71/103), prev=(83/103, 6.5/103), copies=(93/103, 8/103),
    y_bottom=23/62, y_height=33/62,
    lab_x=7.8/103, cls_lab_x=1/103, axis_label_y=8/62,
    title_y=1-0.7/62, header_y=1-3/62,
)
LEG = dict(sw=(1.3/103, 1.3/62), col_gap=1.2/103)
FS = dict(row=6.0, col=6.0, tick=6.0, group=6.0, header=6.0, title=6.0,
          axis=6.0, legend=6.0, legend_title=6.0)

CLASS_ORDER = {"core": 0, "soft_core": 1, "shell": 2, "unique": 3}
CLASS_LABEL = {"core": "Core", "soft_core": "Soft-core",
               "shell": "Shell", "unique": "Unique"}


# --------------------------------------------------------------------------------------
# data
# --------------------------------------------------------------------------------------
def load():
    meta, order = {}, []
    for line in (DATA / "v20260912_family_order.tsv").read_text().splitlines()[1:]:
        p = line.split("\t")
        order.append(p[0])
        meta[p[0]] = dict(cls=p[1], major_type=p[2],
                          prevalence=int(p[3]), copies=int(p[4]))
    materials = [l.split("\t")[0] for l in
                 (DATA / "v20260912_material_order.tsv").read_text().splitlines()[1:]]
    total = {l.split("\t")[0]: int(l.split("\t")[1]) for l in
             (DATA / "v20260912_material_order.tsv").read_text().splitlines()[1:]}
    matrix = {}
    lines = (DATA / "st20_matrix.tsv").read_text().splitlines()
    for line in lines[1:]:
        p = line.split("\t")
        by_material = dict(zip(lines[0].split("\t")[1:], map(int, p[1:])))
        matrix[p[0]] = [by_material[m] for m in materials]
    M = np.array([matrix[f] for f in order])
    assert M.shape == (52, 33) and np.all((M >= 0) & (M <= 3))
    assert M.sum() == 814 and (M > 0).sum() == 795 and (M >= 2).sum() == 18
    assert len(set(order)) == 52 and len(set(materials)) == 33
    for j, m in enumerate(materials):
        assert M[:, j].sum() == total[m], (m, M[:, j].sum(), total[m])
    for i, f in enumerate(order):
        assert M[i].sum() == meta[f]['copies']
        assert (M[i] > 0).sum() == meta[f]['prevalence']
    for c, n, copies in [('core',9,314),('soft_core',2,66),('shell',29,420),('unique',12,14)]:
        fs = [f for f in order if meta[f]['cls'] == c]
        assert len(fs) == n and sum(meta[f]['copies'] for f in fs) == copies
    return meta, order, materials, total, matrix


def type6(major_type: str) -> str:
    if major_type == "saccharide":
        return "Saccharide"
    if major_type == "cyclopeptide":
        return "Cyclopeptide"
    if major_type == "putative":
        return "Putative"
    if major_type == "polyketide":
        return "Polyketide"
    if major_type in SEP_FAMILIES:
        return "Fatty acid"
    return "Others"


def family_id(fam: str) -> str:
    return "F%d" % int(fam.split("_")[1])


def display_order(meta, order) -> list[str]:
    """class (unique->core), and inside a class the reverse of the source sorting rule
    (major type, prevalence desc, copies desc, table order)."""
    types, seen = meta, []
    for f in order:
        if meta[f]["major_type"] not in seen:
            seen.append(meta[f]["major_type"])
    rank = {t: i for i, t in enumerate(seen)}
    idx = {f: i for i, f in enumerate(order)}
    ordered = sorted(order, key=lambda f: (CLASS_ORDER[meta[f]["cls"]],
                                           rank[meta[f]["major_type"]],
                                           -meta[f]["prevalence"], -meta[f]["copies"], idx[f]))
    return list(reversed(ordered))


# --------------------------------------------------------------------------------------
# figure
# --------------------------------------------------------------------------------------
def build(rows, meta, materials, total, matrix, out_dir):
    y0, hh = L["y_bottom"], L["y_height"]
    n_row, n_col = len(rows), len(materials)
    fig = plt.figure(figsize=(FIG_W, FIG_H))
    fig.patch.set_facecolor("white")

    ax_c = fig.add_axes([L["strip_class"][0], y0, L["strip_class"][1], hh])
    ax_t = fig.add_axes([L["strip_type"][0], y0, L["strip_type"][1], hh])
    ax_h = fig.add_axes([L["heat"][0], y0, L["heat"][1], hh])
    ax_p = fig.add_axes([L["prev"][0], y0, L["prev"][1], hh])
    ax_b = fig.add_axes([L["copies"][0], y0, L["copies"][1], hh])
    ov = fig.add_axes([0, 0, 1, 1], zorder=5)          # overlay for lines / manual legends
    ov.set_axis_off()
    ov.set_xlim(0, 1)
    ov.set_ylim(0, 1)

    # ---- matrix ---------------------------------------------------------------
    M = np.array([[min(matrix[f][materials.index(m)], 2) for m in materials] for f in rows])
    ax_h.pcolormesh(np.arange(n_col + 1), np.arange(n_row + 1), M,
                    cmap=ListedColormap([C_ABSENT, C_PRESENT, C_MULTI]),
                    vmin=0, vmax=2, edgecolors="#E1E7E7", linewidth=0.12)
    ax_h.set_xlim(0, n_col)
    ax_h.set_ylim(0, n_row)
    ax_h.invert_yaxis()                                  # row 0 at the top
    ax_h.set_xticks(np.arange(n_col) + 0.5)
    ax_h.set_xticklabels(materials, rotation=90, fontsize=FS["col"],
                         ha="right", va="center", rotation_mode="anchor")
    for i, tick in enumerate(ax_h.get_xticklabels()):
        tick.set_color("#444444")
        tick.set_fontweight("normal")
    ax_h.set_yticks([])
    ax_h.tick_params(axis="x", length=0, pad=1.5)
    for sp in ax_h.spines.values():
        sp.set_linewidth(0.5)
        sp.set_color("#3A3A3A")

    # ---- strips ---------------------------------------------------------------
    for ax, colors, w in ((ax_c, [CLASS_COLORS[CLASS_LABEL[meta[f]["cls"]]] for f in rows],
                           L["strip_class"][1]),
                          (ax_t, [TYPE_COLORS[type6(meta[f]["major_type"])] for f in rows],
                           L["strip_type"][1])):
        for i, col in enumerate(colors):
            ax.add_patch(Rectangle((0, i), 1, 1, facecolor=col, edgecolor="none"))
        ax.set_xlim(0, 1)
        ax.set_ylim(n_row, 0)                          # same orientation as the heatmap
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_axis_off()

    # ---- class group labels, row labels, dashed separators ---------------------
    blocks, start = [], 0
    for i in range(1, n_row + 1):
        if i == n_row or meta[rows[i]]["cls"] != meta[rows[start]]["cls"]:
            blocks.append((start, i - 1))
            start = i
    lab = lambda r: y0 + hh * (n_row - r - 0.5) / n_row
    group_rows = []
    for a, b in blocks:
        cls = CLASS_LABEL[meta[rows[a]]["cls"]]
        centre = (a + b) / 2
        group_rows.append(round(centre))
    # row labels: block start / middle / end, de-duplicated and kept clear of the
    # class-group labels so nothing collides vertically
    want = []
    for a, b in blocks:
        want += [a, (a + b) // 2, b]
    keep = []
    for r in sorted(set(want)):
        if any(abs(r - k) <= 3 for k in keep):
            continue                     # too close to an accepted label
        keep.append(r)
    for r in keep:
        ov.text(L["lab_x"], lab(r), family_id(rows[r]), ha="right", va="center",
                fontsize=FS["row"], color="#222222")
    x_lo = L["strip_type"][0]
    x_hi = L["heat"][0] + L["heat"][1] + 0.002
    for a, b in blocks[:-1]:
        yf = y0 + hh * (n_row - b - 1) / n_row
        ov.plot([x_lo, x_hi], [yf, yf], **DASH_SEP, solid_capstyle="butt")

    ov.text(3/103, L["header_y"], "Class",
            ha="center", va="center", fontsize=FS["header"], fontweight="bold")
    ov.text(9/103, L["header_y"], "Type",
            ha="center", va="center", fontsize=FS["header"], fontweight="bold")

    # ---- prevalence panel ------------------------------------------------------
    prev = [meta[f]["prevalence"] for f in rows]
    ax_p.hlines(np.arange(n_row) + 0.5, 0, 33, color=GUIDE, linewidth=0.5, zorder=1)
    ax_p.scatter(prev, np.arange(n_row) + 0.5, s=3.0, color=DOT_INK, linewidths=0, zorder=3, clip_on=False)
    ax_p.set_xlim(0, 33)
    ax_p.set_ylim(n_row, 0)
    ax_p.set_xticks([0, 33])
    ax_p.set_yticks([])
    ax_p.tick_params(axis="x", labelsize=FS["tick"], length=2.2, width=0.5,
                     color="#3A3A3A", pad=2)
    for sp in ("top", "right", "left"):
        ax_p.spines[sp].set_visible(False)
    ax_p.spines["bottom"].set_linewidth(0.5)
    ov.text(L["prev"][0] + L["prev"][1] / 2, L["title_y"], "Prevalence\n(n = 33)",
            ha="center", va="top", fontsize=FS["title"], fontweight="bold",
            linespacing=1.25)

    # ---- BGC copies panel ------------------------------------------------------
    copies = [meta[f]["copies"] for f in rows]
    ax_b.barh(np.arange(n_row) + 0.5, copies, height=0.72, color=BAR_GOLD,
              edgecolor="none", zorder=2)
    ax_b.axvline(48, color=BAR_GOLD, linewidth=0.7, linestyle=(0, (2.4, 1.8)), zorder=1)
    ax_b.set_xlim(0, 55.5)                      # blank right margin holds the F1 callout
    ax_b.set_ylim(n_row, 0)
    ax_b.set_xticks([0, 50])
    ax_b.set_yticks([])
    ax_b.tick_params(axis="x", labelsize=FS["tick"], length=2.2, width=0.5,
                     color="#3A3A3A", pad=2)
    for sp in ("top", "right", "left"):
        ax_b.spines[sp].set_visible(False)
    ax_b.spines["bottom"].set_linewidth(0.5)
    top = int(np.argmax(copies))
    ov.text(L["copies"][0] + L["copies"][1], y0 - 5.5/FIG_H_MM,
            f"{family_id(rows[top])} ({copies[top]})", ha="right", va="center",
            fontsize=FS["row"], color="#506D86", fontweight="bold")
    ov.text(L["copies"][0] + L["copies"][1] / 2, L["title_y"], "Copies\n(total)",
            ha="center", va="top", fontsize=FS["title"], fontweight="bold",
            linespacing=1.25)

    # ---- bottom legends (manual, exact positions) ------------------------------
    counts_cls = {c: sum(1 for f in rows if CLASS_LABEL[meta[f]["cls"]] == c)
                  for c in CLASS_COLORS}
    counts_typ = {t: sum(1 for f in rows if type6(meta[f]["major_type"]) == t)
                  for t in TYPE_COLORS}

    def swatch(x, y, color):
        ov.add_patch(Rectangle((x, y - LEG["sw"][1] / 2), *LEG["sw"],
                               facecolor=color, edgecolor="none"))

    def entry(x, y, color, label):
        swatch(x, y, color)
        ov.text(x + LEG["sw"][0] + LEG["col_gap"] * 0.35, y, label, ha="left",
                va="center", fontsize=FS["legend"], color="#111111")

    # Compact aligned key: copy/class row and two rows of type entries.
    def title(x_mm, y_mm, label):
        ov.text(x_mm/FIG_W_MM, y_mm/FIG_H_MM, label, fontsize=6,
                fontweight="bold", ha="left", va="center")
    title(1.5, 7.5, "Copy no.")
    for x,c,label in [(12,C_ABSENT,"0"),(18,C_PRESENT,"1"),(24,C_MULTI,"≥2")]:
        entry(x/103,7.5/62,c,label)
    title(34,7.5,"Class")
    for x,k in zip([43,56,73,89],["Core","Soft-core","Shell","Unique"]):
        entry(x/103,7.5/62,CLASS_COLORS[k],f"{k} ({counts_cls[k]})")
    title(1.5,5,"BGC type")
    for x,keys in [(17,["Saccharide","Fatty acid"]),
                   (46,["Cyclopeptide","Polyketide"]),(75,["Putative","Others"])]:
        for y,k in zip([5,2.5],keys):
            entry(x/103,y/62,TYPE_COLORS[k],f"{k} ({counts_typ[k]})")

    plt.rcParams.update({"pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none"})
    FIG.mkdir(parents=True, exist_ok=True)
    stem = FIG / "Fig_panBGC_overview_v5"
    fig.savefig(f"{stem}.pdf")
    fig.savefig(f"{stem}.svg")
    fig.savefig(f"{stem}.png", dpi=600)
    fig.savefig(QC / "v5_preview_150dpi.png", dpi=150)
    return fig, stem


# --------------------------------------------------------------------------------------
def qa(fig):
    """text overlap check + rendered font sizes (final physical size)."""
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    fw_px, fh_px = fig.canvas.get_width_height()
    items = []
    for ax in fig.axes:
        texts = list(ax.texts) + [ax.title, ax.xaxis.label, ax.yaxis.label]
        if ax.axison:
            texts += list(ax.get_xticklabels()) + list(ax.get_yticklabels())
        for t in texts:
            if t.get_text() and t.get_visible():
                bb = t.get_window_extent(renderer=r)
                items.append((t.get_text(), bb, t.get_fontsize(),
                              round(bb.x0 / fw_px, 4), round(bb.y0 / fh_px, 4)))
    bad = []
    for i in range(len(items)):
        for j in range(i + 1, len(items)):
            a, b = items[i][1], items[j][1]
            ix = min(a.x1, b.x1) - max(a.x0, b.x0)
            iy = min(a.y1, b.y1) - max(a.y0, b.y0)
            if ix > 2.0 and iy > 2.0:
                bad.append((items[i][0][:24], items[j][0][:24], round(ix, 1), round(iy, 1),
                            items[i][2], items[i][3], items[j][3]))
    outside = [(t[0][:28], t[3], round(t[1].x1 / fw_px, 4))
               for t in items if t[1].x0 < -0.5 or t[1].x1 > fw_px + 0.5
               or t[1].y0 < -0.5 or t[1].y1 > fh_px + 0.5]
    biggest = max(items, key=lambda t: t[2])
    return bad, outside, biggest


def main():
    for d in (FIG, TAB, QC):
        d.mkdir(parents=True, exist_ok=True)
    from matplotlib import font_manager
    fontdir = Path('${DATA_DIR2}/.local/share/fonts/msttcorefonts')
    for name in ('Arial.TTF', 'Arialbd.TTF', 'Ariali.TTF'):
        font_manager.fontManager.addfont(str(fontdir / name))
    assert 'Arial' in font_manager.FontProperties(fname=font_manager.findfont('Arial', fallback_to_default=False)).get_name()
    meta, order, materials, total, matrix = load()
    rows = display_order(meta, order)
    if not (len(rows) == 52 and len(materials) == 33):
        raise SystemExit("unexpected matrix shape")
    if sum(meta[f]["copies"] for f in rows) != sum(sum(r) for r in matrix.values()):
        raise SystemExit("BGC copy totals are inconsistent")
    fig, stem = build(rows, meta, materials, total, matrix, FIG)

    TAB.mkdir(exist_ok=True)
    (TAB / 'v5_column_order.tsv').write_text('Material\tTotal_Clusters\n' + ''.join(f'{m}\t{total[m]}\n' for m in materials))
    (TAB / 'v5_display_matrix.tsv').write_text('PanBGC_Family\t' + '\t'.join(materials) + '\n' + ''.join(f + '\t' + '\t'.join(map(str, matrix[f])) + '\n' for f in rows))
    with open(TAB / "v5_row_order.tsv", "w") as fh:
        fh.write("Display_Row\tFamily_ID\tPanBGC_Family\tClass\tMajor_Type\tType_Legend\t"
                 "Prevalence\tBGC_Copies\n")
        for i, f in enumerate(rows):
            fh.write(f"{i + 1}\t{family_id(f)}\t{f}\t{CLASS_LABEL[meta[f]['cls']]}\t"
                     f"{meta[f]['major_type']}\t{type6(meta[f]['major_type'])}\t"
                     f"{meta[f]['prevalence']}\t{meta[f]['copies']}\n")

    bad, outside, biggest = qa(fig)
    text_sizes = [t.get_fontsize() for ax in fig.axes
                  for t in (list(ax.texts) + (list(ax.get_xticklabels()) + list(ax.get_yticklabels()) if ax.axison else []))
                  if t.get_visible() and t.get_text()]
    assert text_sizes and all(6 <= z <= 8 for z in text_sizes), text_sizes
    fw, fh = fig.get_size_inches()
    used = sorted(set(FS.values()))
    size_mm = (fw * 25.4, fh * 25.4)
    n_text = 0
    QC.mkdir(exist_ok=True)
    lines = [
        "# Redraw QA — pan-BGC overview",
        "",
        f"- figure: {size_mm[0]:.1f} x {size_mm[1]:.1f} mm (requested width 103 mm)",
        "- nominal font sizes used (pt): 6.0 pt for all labels at the final 103 × 62 mm size; "
        "all within the 5–8 pt range required for this figure size",
        f"- text–text overlaps detected: {len(bad)}",
        f"- text objects extending beyond the canvas: {len(outside)}",
    ]
    if outside:
        lines += [""] + [f"  - {t} (x0={x0}, x1={x1})" for t, x0, x1 in outside[:20]]
    if bad:
        lines += ["", "| text A | text B | overlap x | overlap y | A pt | A x0 | B x0 |",
                  "|---|---|---|---|---|---|---|"]
        lines += [f"| {a} | {b} | {ix} | {iy} | {pa} | {xa} | {xb} |"
                  for a, b, ix, iy, pa, xa, xb in bad[:40]]
    (QC / "QA_v5_graphics.md").write_text("\n".join(lines) + "\n")
    (QC / "v5_run_manifest.json").write_text(json.dumps({
        "script": str(Path(__file__).resolve()),
        "python": platform.python_version(),
        "matplotlib": matplotlib.__version__,
        "inputs": {p.name: hashlib.md5(p.read_bytes()).hexdigest()
                   for p in sorted(DATA.glob("*.tsv"))},
        "outputs": {p.name: hashlib.md5(p.read_bytes()).hexdigest()
                    for p in sorted(FIG.glob("Fig_panBGC_overview_v5.*"))},
        "figure_mm": [round(size_mm[0], 2), round(size_mm[1], 2)],
    }, indent=2) + "\n")
    assert not bad and not outside, (bad, outside)
    print(f"saved: {stem}.pdf/.svg/.png")
    print(f"rows={len(rows)} materials={len(materials)} "
          f"copies={sum(meta[f]['copies'] for f in rows)} "
          f"overlaps={len(bad)} outside_canvas={len(outside)} "
          f"largest_text={biggest[0][:18]!r}@{biggest[2]}pt")


if __name__ == "__main__":
    main()
