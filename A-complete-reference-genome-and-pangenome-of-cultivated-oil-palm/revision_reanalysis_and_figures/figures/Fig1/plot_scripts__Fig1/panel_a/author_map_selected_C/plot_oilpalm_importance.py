#!/usr/bin/env python3
"""Three oil-palm importance concepts, 183 mm wide, PDF + 600 dpi PNG.

Data: existing FAOSTAT World extract (2023), Descals cover grids (2019).
Shares explicitly use the 9 crop / 11 oil categories in the local extract.
The share ratio is not interpreted as agronomic yield or land saving.
Usage: python3 plot_oilpalm_importance.py [A|B|C]
"""
from pathlib import Path
import json
import os
import sys

ROOT = Path(globals().get("__file__", Path.cwd() / "plot_oilpalm_importance.py")).resolve().parent
OUT = ROOT / "importance_designs"
OUT.mkdir(exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(ROOT / "qa" / ".mplconfig"))
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.patches import Rectangle, FancyBboxPatch, ConnectionPatch
from matplotlib.lines import Line2D

# Reuse approved map data, label coordinates and helpers, without constructing
# or exporting the earlier figure. It remains an unchanged companion script.
map_path = ROOT / "plot_fig1a_density.py"
map_ns = {"__file__": str(map_path)}
exec(compile(map_path.read_text().split("# A narrower overview")[0], str(map_path), "exec"), map_ns)
draw_region, scale_bar = map_ns["draw_region"], map_ns["scale_bar"]
MAIN, REGIONS, DENS05 = map_ns["MAIN"], map_ns["REGIONS"], map_ns["DENS05"]
CMAP, OCEAN, INK = map_ns["CMAP"], map_ns["OCEAN"], map_ns["INK"]
ACCENT, MUTED, REST, LINE = "#C65B28", "#67777E", "#E1E6E7", "#C9D2D5"
W = 183 / 25.4
matplotlib.rcParams.update({
    "pdf.fonttype": 42, "ps.fonttype": 42, "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 6.5, "text.color": INK, "axes.labelcolor": INK,
    "xtick.color": MUTED, "ytick.color": INK, "axes.linewidth": .4,
    "savefig.facecolor": "white",
})

YEAR = 2023
fao = pd.read_csv(ROOT / "faostat_world_oilcrops_2019_latest.tsv", sep="\t")
fao = fao[(fao["Area"] == "World") & (fao["Year"] == YEAR)]
area = fao[fao["Element Code"] == 5312].set_index("Item_short")["Value"]
oil = fao[(fao["Element Code"] == 5510) & fao["Item_short"].str.contains(" oil", case=False)].set_index("Item_short")["Value"]
assert area.index.is_unique and oil.index.is_unique
assert len(area) == 9 and len(oil) == 11
assert np.all(area > 0) and np.all(oil > 0)
PA = float(area["Oil palm fruit"] / area.sum() * 100)
PO = float((oil["Palm oil"] + oil["Palm kernel oil"]) / oil.sum() * 100)
AM = float(area["Oil palm fruit"] / 1e6)
OM = float((oil["Palm oil"] + oil["Palm kernel oil"]) / 1e6)
CROP_ROWS = [
    ("Oil palm", PA, PO),
    ("Soybean", float(area["Soya beans"] / area.sum() * 100), float(oil["Soybean oil"] / oil.sum() * 100)),
    ("Rapeseed", float(area["Rapeseed"] / area.sum() * 100), float(oil["Rapeseed oil"] / oil.sum() * 100)),
    ("Sunflower", float(area["Sunflower seed"] / area.sum() * 100), float(oil["Sunflower oil"] / oil.sum() * 100)),
]


def tx(fig, x, y, s, size=6.5, **kwargs):
    """Place text in physical inches, measured from bottom left."""
    defaults = dict(fontsize=size, ha="left", va="center", color=INK)
    defaults.update(kwargs)
    return fig.text(x / W, y / fig.get_figheight(), s, **defaults)


def rect(fig, x, y, w, h, color, **kwargs):
    patch = Rectangle((x, y), w, h, transform=fig.dpi_scale_trans,
                      facecolor=color, edgecolor="none", **kwargs)
    fig.add_artist(patch)
    return patch


def share_bar(fig, x, y, w, h, pct):
    rect(fig, x, y, w, h, REST, zorder=1)
    rect(fig, x, y, w * pct / 100, h, ACCENT, zorder=2)


def heading(fig, title, subtitle=None):
    h = fig.get_figheight()
    tx(fig, .10, h - .12, title, 9, fontweight="bold")
    if subtitle:
        tx(fig, .10, h - .30, subtitle, 6.5, color=MUTED)


def footnote(fig, with_map=False):
    tx(fig, .10, .16, "Shares: 9 major oil crops (area) and 11 vegetable oils (output). Oil palm output = palm oil + palm-kernel oil.", 5.5, color=MUTED)
    source = "FAOSTAT, World, 2023."
    if with_map:
        source += "  Map: Descals et al. (2021), closed-canopy oil palm in 2019."
    tx(fig, .10, .06, source, 5.5, color=MUTED)


def add_colorbar(fig, x, y, w, title="Oil-palm cover (%)", title_y=None):
    h = fig.get_figheight()
    cax = fig.add_axes([x / W, y / h, w / W, .045 / h])
    cb = fig.colorbar(plt.cm.ScalarMappable(cmap=CMAP, norm=LogNorm(1, 100)),
                      cax=cax, orientation="horizontal")
    cb.set_ticks([1, 10, 100]); cb.set_ticklabels(["1", "10", "100"])
    cb.ax.minorticks_off()
    cb.ax.tick_params(labelsize=5.5, length=1.2, width=.35, pad=1.1)
    cb.outline.set_linewidth(.3); cb.outline.set_edgecolor(MUTED)
    tx(fig, x + w / 2, title_y if title_y is not None else y + .10, title,
       5.5, ha="center")
    return cax


def design_a():
    """Message-led: two comparable 100% bars, large exact labels."""
    fig = plt.figure(figsize=(W, 2.38), dpi=100)
    heading(fig, "Oil palm's contribution to vegetable-oil production",
            "Harvested area and oil output among major oil crops and oils  |  World, 2023")
    columns = [(.17, "Harvested area", PA, f"{AM:.1f} million ha of {area.sum()/1e6:.1f} million ha", "Other major oil crops"),
               (3.86, "Vegetable-oil output", PO, f"{OM:.1f} million tonnes of {oil.sum()/1e6:.1f} million tonnes", "Other vegetable oils")]
    for x, label, value, absolute, other in columns:
        tx(fig, x, 1.76, label, 8, fontweight="bold")
        tx(fig, x, 1.37, f"{value:.1f}%", 23, color=ACCENT, fontweight="bold")
        tx(fig, x, 1.06, absolute, 6.5, color=MUTED)
        share_bar(fig, x, .71, 3.10, .20, value)
        tx(fig, x, .56, "Oil palm", 6.5, color=ACCENT, fontweight="bold")
        tx(fig, x + 3.10, .56, other, 6.5, color=MUTED, ha="right")
        tx(fig, x + 3.10, .98, "100%", 5.5, color=MUTED, ha="right")
    footnote(fig)
    return fig


def design_b():
    """Geography plus a crop comparison: each graphic answers one question."""
    fig = plt.figure(figsize=(W, 3.10), dpi=100)
    tx(fig, .08, 2.98, "a", 9, fontweight="bold")
    tx(fig, .24, 2.98, "Where oil palm is grown", 8, fontweight="bold")
    tx(fig, 7.06, 2.98, "Closed-canopy cover, 2019", 6.5, ha="right", color=MUTED)
    x, y, w = .09, 1.55, 7.01
    mh = w * 54 / 286
    ax = fig.add_axes([x / W, y / 3.10, w / W, mh / 3.10])
    draw_region(ax, MAIN)
    for txt, lx, ly in [("Americas", -103, -17), ("Africa", -12, -13), ("Asia–Pacific", 99, -16)]:
        ax.text(lx, ly, txt, fontsize=6.5, fontweight="bold", color=INK, va="center")
    for txt, lx, ly in [("Atlantic\nOcean", -28, -14), ("Indian\nOcean", 80, -3), ("Pacific\nOcean", 155, 22)]:
        ax.text(lx, ly, txt, fontsize=5.5, style="italic", color="#7893A1", ha="center", va="center")
    # Compact legend over Indian Ocean; no production bubbles on the map.
    fx = lambda lon: x + (lon - MAIN[0]) / 286 * w
    fy = lambda lat: y + (lat - MAIN[2]) / 54 * mh
    add_colorbar(fig, fx(58), fy(-22), fx(86) - fx(58), title_y=fy(-17.2))
    tx(fig, .08, 1.35, "b", 9, fontweight="bold")
    tx(fig, .24, 1.35, "Area and oil-output shares among major crops", 8, fontweight="bold")
    a = fig.add_axes([.98 / W, .44 / 3.10, 4.67 / W, .69 / 3.10])
    for i, (name, av, ov) in enumerate(CROP_ROWS):
        yy = 3 - i
        color = ACCENT if i == 0 else "#81939B"
        a.plot([av, ov], [yy, yy], color=color, lw=1.0, zorder=2)
        a.plot(av, yy, marker="o", ms=4.3, mfc="white", mec=color, mew=.9, zorder=3)
        a.plot(ov, yy, marker="o", ms=4.3, mfc=color, mec=color, mew=.6, zorder=3)
        # Text in aligned columns avoids ambiguity for short connecting lines.
        yy_fig = .44 + (yy + .5) / 4 * .69
        tx(fig, 6.10, yy_fig, f"{av:.1f}%", 6.5, ha="right", color=color, fontweight="bold" if i == 0 else "normal")
        tx(fig, 6.94, yy_fig, f"{ov:.1f}%", 6.5, ha="right", color=color, fontweight="bold" if i == 0 else "normal")
    a.set_xlim(0, 50); a.set_ylim(-.5, 3.5)
    a.set_yticks([3, 2, 1, 0]); a.set_yticklabels([r[0] for r in CROP_ROWS], fontsize=6.5)
    a.get_yticklabels()[0].set_color(ACCENT); a.get_yticklabels()[0].set_fontweight("bold")
    a.set_xticks([0, 10, 20, 30, 40, 50]); a.tick_params(axis="x", labelsize=5.5, length=2, pad=1)
    a.tick_params(axis="y", length=0, pad=4)
    a.grid(axis="x", color="#E9ECEC", lw=.4, zorder=0)
    for side in ["top", "right", "left"]: a.spines[side].set_visible(False)
    a.spines["bottom"].set_color(LINE)
    a.set_xlabel("Share within the respective area or oil-output total (%)", fontsize=5.5, labelpad=2)
    tx(fig, 6.10, 1.18, "Area", 6.0, ha="right", color=MUTED)
    tx(fig, 6.94, 1.18, "Oil output", 6.0, ha="right", color=MUTED)
    # Explicit shape semantics, separate from the percentage columns.
    fig.legend(handles=[Line2D([], [], marker="o", color="none", mfc="white", mec=MUTED, ms=4, label="Area"),
                        Line2D([], [], marker="o", color="none", mfc=MUTED, mec=MUTED, ms=4, label="Oil output")],
               loc="center right", bbox_to_anchor=(.98, 1.35 / 3.10), ncol=2,
               frameon=False, fontsize=5.5, handletextpad=.4, columnspacing=.8)
    footnote(fig, with_map=True)
    return fig


def design_c():
    """Retain the three regional zooms; add a compact contribution summary."""
    H = 3.10
    fig = plt.figure(figsize=(W, H), dpi=100)
    heading(fig, "Oil palm: growing regions and contribution to vegetable-oil output")
    tx(fig, .10, 2.73, "Global distribution, 2019", 8, fontweight="bold")
    tx(fig, 5.07, 2.73, "Contribution, 2023", 8, fontweight="bold")
    mx, my, mw = .10, 1.66, 4.72
    mh = mw * 54 / 286
    ax = fig.add_axes([mx / W, my / H, mw / W, mh / H], zorder=1)
    draw_region(ax, MAIN)
    for n, name, bb, col, _ in REGIONS:
        lo0, lo1, la0, la1 = bb
        ax.add_patch(Rectangle((lo0, la0), lo1-lo0, la1-la0, facecolor="none",
                               edgecolor=col, lw=.45, ls=(0, (3.5, 2.5)), alpha=.85, zorder=7))
        ax.text(lo0 + 4, la0 - 4.7, str(n), fontsize=5.5, fontweight="bold",
                color=col, ha="center", va="center", zorder=8,
                bbox=dict(boxstyle="circle,pad=.14", facecolor=OCEAN, edgecolor=col, linewidth=.4))
    add_colorbar(fig, 3.24, 2.69, 1.47, title_y=2.83)
    sx, sw = 5.07, 1.98
    for label, pct, label_y, bar_y in [("Harvested area", PA, 2.39, 2.18),
                                     ("Vegetable-oil output", PO, 1.98, 1.77)]:
        tx(fig, sx, label_y, label, 6.5)
        tx(fig, sx + sw, label_y, f"{pct:.1f}%", 9, fontweight="bold", color=ACCENT, ha="right")
        share_bar(fig, sx, bar_y, sw, .105, pct)
    tx(fig, sx, 1.60, "Share of the selected major-crop / oil totals", 5.5, color=MUTED)

    ih, header, bottom, x = 1.10, .155, .23, .055
    widths = [ih * (bb[1]-bb[0])/(bb[3]-bb[2]) for _,_,bb,_,_ in REGIONS]
    gap = (W - 2*x - sum(widths))/2
    for (n, name, bb, col, pts), w in zip(REGIONS, widths):
        card = FancyBboxPatch((x, bottom), w, ih+header,
                             boxstyle="round,pad=0,rounding_size=.030",
                             facecolor=OCEAN, edgecolor="none", transform=fig.dpi_scale_trans, zorder=.4)
        fig.add_artist(card)
        axi = fig.add_axes([x/W, bottom/H, w/W, ih/H], zorder=1)
        draw_region(axi, bb, grid=DENS05); axi.patch.set_visible(False)
        for coll in axi.collections: coll.set_clip_path(card)
        frame = FancyBboxPatch((x, bottom), w, ih+header,
                              boxstyle="round,pad=0,rounding_size=.030", facecolor="none",
                              edgecolor=LINE, linewidth=.4, transform=fig.dpi_scale_trans, zorder=12)
        fig.add_artist(frame)
        hy = bottom+ih+header*.52
        tx(fig, x+.085, hy, str(n), 5.5, fontweight="bold", color=col, ha="center",
           bbox=dict(boxstyle="circle,pad=.18", facecolor=OCEAN, edgecolor=col, linewidth=.45))
        tx(fig, x+.18, hy, name, 6.5, fontweight="bold", color=col)
        zoom = (w/(bb[1]-bb[0]))/(mw/286)
        tx(fig, x+w-.065, hy, f"{zoom:.1f}×", 5.5, color=MUTED, ha="right")
        for lab, px, py, lx, ly, ha in pts:
            axi.plot(px, py, "o", ms=1.15, color=MUTED, zorder=8)
            axi.annotate(lab, (px, py), xytext=(lx, ly), fontsize=5.5, color=INK,
                         ha=ha, va="center", zorder=9,
                         arrowprops=dict(arrowstyle="-", color=MUTED, lw=.3, shrinkA=2, shrinkB=1.5,
                                         connectionstyle="angle3,angleA=0,angleB=90" if lab=="Gabon" else "arc3"))
        scale_bar(axi, bb[0]+1.9, bb[2]+3.6, [0,500,1000], lat_ref=(bb[2]+bb[3])/2,
                  label_dx=-1.6 if n==3 else 0)
        source = (bb[1], bb[2])
        dest_x = x+w-.06 if n==1 else (x+w*.57 if n==2 else x+.05)
        fig.add_artist(ConnectionPatch(source, (dest_x/W, (bottom+ih+header)/H),
                                       coordsA=ax.transData, coordsB=fig.transFigure, axesA=ax,
                                       arrowstyle="-", connectionstyle="arc3,rad=.08",
                                       lw=.4, color=col, alpha=.35, zorder=6, clip_on=False))
        x += w+gap
    footnote(fig, with_map=True)
    return fig


BUILDERS = {"A": design_a, "B": design_b, "C": design_c}
choice = sys.argv[1].upper() if len(sys.argv)>1 and sys.argv[1].upper() in BUILDERS else None
figures = {}
for key in ([choice] if choice else BUILDERS):
    fig = BUILDERS[key]()
    figures[key] = fig
    fig.savefig(OUT / f"Fig_oilpalm_importance_{key}.pdf", dpi=600)
    fig.savefig(OUT / f"Fig_oilpalm_importance_{key}.png", dpi=600)
    print(f"{key}: {W*25.4:.3f} × {fig.get_figheight()*25.4:.3f} mm, 600 dpi")

summary = {
    "year": YEAR, "area_categories": area.to_dict(), "oil_categories": oil.to_dict(),
    "oil_palm_area_mha": AM, "oil_palm_oil_mt": OM,
    "area_total_mha": float(area.sum()/1e6), "oil_total_mt": float(oil.sum()/1e6),
    "oil_palm_area_share_pct": PA, "oil_palm_oil_share_pct": PO,
    "denominator_scope": "9 selected major oil crops; 11 selected vegetable oils, not complete FAOSTAT commodity coverage",
    "map_year": 2019, "oil_palm_output_definition": "Palm oil + palm kernel oil",
}
(OUT / "data_summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2))
pd.DataFrame(CROP_ROWS, columns=["Crop", "Area_share_pct", "Oil_output_share_pct"]).to_csv(OUT / "crop_share_comparison.csv", index=False)
