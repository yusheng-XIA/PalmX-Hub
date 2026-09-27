#!/usr/bin/env python3
"""Fig. 1a, 2026-09-18: compact overview and softly framed regional zooms.

Based on plot_fig1a_density_v6.py; same grids, geographic extents, equal-aspect
longitude/latitude projection and density normalization. Nature-style draft,
183 mm wide, 600 dpi PNG; PDF retains vector geography, lines and TrueType text.
Run from this directory, including with the supplied check_text_overlap.py.
"""
from pathlib import Path
import os

ROOT = Path(globals().get("__file__", Path.cwd() / "plot_fig1a_density.py")).resolve().parent
QA = ROOT / "qa"
QA.mkdir(exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(QA / ".mplconfig"))
import numpy as np
import geopandas as gpd
from shapely.geometry import box
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, LinearSegmentedColormap
from matplotlib.patches import Rectangle, FancyBboxPatch, ConnectionPatch

matplotlib.rcParams.update({
    "pdf.fonttype": 42, "ps.fonttype": 42,
    "font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 5.5, "axes.linewidth": 0.4,
    "savefig.facecolor": "white",
})
# Exactly 183 mm; 7.2 inches in the brief is its rounded equivalent.
W = 183 / 25.4
# Stronger land/ocean separation, with restrained coastlines and soft frames.
INK, LANDC, OCEAN, BORDER, COAST = "#303B40", "#CDD3CE", "#F3F7FA", "#FFFFFF", "#ADB9B1"
# Original v6 density ramp. Contrast refinement applies to the basemap.
CMAP = LinearSegmentedColormap.from_list("palm", ["#FFF0BF", "#FDCB6E", "#F59440", "#DA452C", "#8C0A22"])
LEAD = "#66747A"
FRAME = {1: "#A65B50", 2: "#467965", 3: "#386E94"}
VMIN, VMAX = 1.0, 100.0


def load(fn):
    with np.load(ROOT / fn) as d:
        dens = (d["ind"] + d["sm"]) / np.maximum(d["tot"], 1) * 100.0
        return np.ma.masked_less(dens, VMIN), d["lon_edges"], d["lat_edges"]


DENS10 = load("descals_density_0p1deg.npz")
DENS05 = load("descals_density_0p05deg.npz")
land_all = gpd.read_file(ROOT / "ne_50m_land.geojson")
bord_all = gpd.read_file(ROOT / "ne_50m_admin_0_boundary_lines_land.geojson")
MAIN = (-118, 168, -27, 27)
# Original anchor coordinates and region boundaries; label positions alone move.
REGIONS = [
    (1, "Latin America", (-104, -44, -13, 21), FRAME[1],
     [("Mexico", -93.2, 17.3, -86.0, 18.0, "left"),
      ("Guatemala", -90.6, 15.6, -102.2, 11.9, "left"),
      ("Honduras", -87.3, 15.4, -77.5, 16.4, "left"),
      ("Costa Rica", -84.0, 9.6, -102.2, 7.0, "left"),
      ("Colombia", -73.6, 8.0, -76.0, 13.9, "left"),
      ("Ecuador", -79.5, 0.6, -102.2, 1.8, "left"),
      ("Peru", -75.6, -7.6, -102.2, -3.5, "left"),
      ("Brazil", -47.9, -2.4, -52.5, 8.2, "center")]),
    (2, "Africa", (-16, 32, -9, 15), FRAME[2],
     [("Côte d'Ivoire", -5.4, 5.9, -14.4, 2.0, "left"),
      ("Ghana", -1.4, 6.1, -1.0, 2.0, "center"),
      ("Nigeria", 6.4, 6.4, 4.7, 2.9, "center"),
      ("Cameroon", 10.0, 4.3, 4.8, -1.3, "center"),
      ("Gabon", 10.4, -0.6, 5.7, -4.4, "center"),
      ("DRC", 22.0, -3.0, 5.7, -7.3, "center")]),
    (3, "Southeast Asia", (92, 156, -12, 16), FRAME[3],
     [("Thailand", 99.5, 9.0, 110.5, 12.9, "left"),
      ("Malaysia", 102.0, 4.0, 106.5, 7.1, "left"),
      ("Sumatra", 101.5, 0.3, 93.5, -4.7, "left"),
      ("Kalimantan", 114.5, -1.0, 112.5, -10.5, "center"),
      ("Philippines", 124.0, 7.6, 130.1, 12.4, "left"),
      ("Papua New Guinea", 147.0, -6.0, 136.0, 4.4, "center")]),
]


def draw_region(ax, bbox, grid=None):
    dens, lon_e, lat_e = DENS10 if grid is None else grid
    lo0, lo1, la0, la1 = bbox
    ax.set_facecolor(OCEAN)
    land = gpd.clip(land_all, box(lo0 - 2, la0 - 2, lo1 + 2, la1 + 2))
    land.plot(ax=ax, facecolor=LANDC, edgecolor=COAST, lw=0.22, zorder=1)
    bord = gpd.clip(bord_all, box(lo0 - 2, la0 - 2, lo1 + 2, la1 + 2))
    bord.plot(ax=ax, color=BORDER, lw=0.3, zorder=2)
    # Crop only the render window: no regridding, interpolation or filtering.
    j0, j1 = max(0, np.searchsorted(lon_e, lo0) - 1), min(len(lon_e) - 1, np.searchsorted(lon_e, lo1) + 1)
    i0, i1 = max(0, np.searchsorted(lat_e, la0) - 1), min(len(lat_e) - 1, np.searchsorted(lat_e, la1) + 1)
    ax.pcolormesh(lon_e[j0:j1 + 1], lat_e[i0:i1 + 1], dens[i0:i1, j0:j1],
                  cmap=CMAP, norm=LogNorm(VMIN, VMAX), shading="flat", zorder=3, rasterized=True)
    ax.set_xlim(lo0, lo1)
    ax.set_ylim(la0, la1)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)


def scale_bar(ax, lon, lat, km_ticks, lat_ref=0.0, end_ha="center", label_dx=0.0):
    """Original longitude scale convention, at each panel's reference latitude."""
    kmdeg = 111.32 * np.cos(np.radians(lat_ref))
    for i, k in enumerate(km_ticks[1:], 1):
        xp = lon + km_ticks[i - 1] / kmdeg
        x1 = lon + k / kmdeg
        ax.add_patch(Rectangle((xp, lat), x1 - xp, 0.27, facecolor=LEAD if i % 2 else "white",
                               edgecolor=LEAD, lw=0.3, zorder=8))
    ax.text(lon, lat - 0.6, "0", ha="center", va="top", fontsize=5.5, color=INK, zorder=8)
    ax.text(lon + km_ticks[-1] / kmdeg + label_dx, lat - 0.6, f"{km_ticks[-1]:,} km",
            ha=end_ha, va="top", fontsize=5.5, color=INK, zorder=8)


def badge(ax, x, y, n, colour):
    return ax.text(x, y, str(n), fontsize=5.5, fontweight="bold", color=colour,
                   ha="center", va="center", zorder=9,
                   bbox=dict(boxstyle="circle,pad=0.18", facecolor=OCEAN,
                             edgecolor=colour, linewidth=0.45))


# A narrower overview makes its three original geographic boxes physically
# smaller; insets stay near the v6 size, making the enlargement explicit.
main_w_frac = 0.86
main_x = (1 - main_w_frac) / 2
main_h_in = W * main_w_frac * (MAIN[3] - MAIN[2]) / (MAIN[1] - MAIN[0])
ins_h = 1.10
header_in, gap_in, top_in, bot_in = 0.155, 0.105, 0.115, 0.075
H = top_in + main_h_in + gap_in + header_in + ins_h + bot_in
fig = plt.figure(figsize=(W, H), dpi=100)
y_main = (bot_in + ins_h + header_in + gap_in) / H
ax = fig.add_axes([main_x, y_main, main_w_frac, main_h_in / H], zorder=1)
draw_region(ax, MAIN)

names = {1: "Americas", 2: "Africa", 3: "Asia–Pacific"}
for n, name, bb, col, _ in REGIONS:
    lo0, lo1, la0, la1 = bb
    ax.add_patch(Rectangle((lo0, la0), lo1 - lo0, la1 - la0, facecolor="none",
                           edgecolor=col, lw=0.4, alpha=0.70,
                           ls=(0, (3.5, 2.5)), joinstyle="round", zorder=7))
    badge(ax, lo0 + 3.5, la0 - 4.0, n, col)
    ax.text(lo0 + 8.1, la0 - 4.0, names[n], fontsize=5.5, fontweight="bold",
            color=col, ha="left", va="center", zorder=8)
for t, x, y in [("Atlantic\nOcean", -28, -14), ("Indian\nOcean", 80, -4),
                ("Pacific\nOcean", 155, 22.3)]:
    ax.text(x, y, t, fontsize=5.5, color="#829BAA", style="italic",
            ha="center", va="center", zorder=8)
scale_bar(ax, -116, -22.8, [0, 2500, 5000], end_ha="right")
ax.annotate("", xy=(-113, -12.9), xytext=(-113, -18.5),
            arrowprops=dict(arrowstyle="-|>", color=LEAD, lw=0.5), zorder=8)
ax.text(-113, -12.0, "N", ha="center", va="bottom", fontsize=5.5, color=INK, zorder=8)
fig.text(0.009, 1 - 0.018 / H, "a", fontsize=9, fontweight="bold", va="top", ha="left", color=INK)


def fx(lon):
    return main_x + (lon - MAIN[0]) / (MAIN[1] - MAIN[0]) * main_w_frac


def fy(lat):
    return y_main + (lat - MAIN[2]) / (MAIN[3] - MAIN[2]) * main_h_in / H


cax = fig.add_axes([fx(58), fy(-21.9), fx(86) - fx(58), fy(-19.9) - fy(-21.9)], zorder=2)
sm = plt.cm.ScalarMappable(cmap=CMAP, norm=LogNorm(VMIN, VMAX))
cb = fig.colorbar(sm, cax=cax, orientation="horizontal")
cb.set_ticks([1, 10, 100])
cb.set_ticklabels(["1", "10", "100"])
cb.ax.tick_params(labelsize=5.5, length=1.4, width=0.35, pad=1.4, colors=INK)
cb.ax.minorticks_off()
cb.outline.set_linewidth(0.3)
cb.outline.set_edgecolor(LEAD)
fig.text(fx(72), fy(-18.5), "Closed-canopy oil palm\n(km² per 100 km²)",
         fontsize=5.5, ha="center", va="bottom", color=INK, linespacing=1.10)

# Equal map heights and identical card/header treatment. Titles occupy an ocean-
# coloured margin above the unchanged geographic extents, never the mapped land.
widths = [ins_h * (bb[1] - bb[0]) / (bb[3] - bb[2]) for _, _, bb, _, _ in REGIONS]
margin_in = 0.055
gap = (W - 2 * margin_in - sum(widths)) / 2
x = margin_in
insets, panel_cards, zoom_factors = [], [], []
for (n, name, bb, col, pts), w in zip(REGIONS, widths):
    card = FancyBboxPatch((x, bot_in), w, ins_h + header_in,
                         boxstyle="round,pad=0,rounding_size=0.030",
                         facecolor=OCEAN, edgecolor="none", transform=fig.dpi_scale_trans, zorder=0.4)
    fig.add_artist(card)
    axi = fig.add_axes([x / W, bot_in / H, w / W, ins_h / H], zorder=1)
    draw_region(axi, bb, grid=DENS05)
    axi.patch.set_visible(False)
    for coll in axi.collections:
        coll.set_clip_path(card)
    # Subtle neutral outline replaces the heavy coloured rectangle.
    frame = FancyBboxPatch((x, bot_in), w, ins_h + header_in,
                          boxstyle="round,pad=0,rounding_size=0.030",
                          facecolor="none", edgecolor="#C9D2D5", linewidth=0.4,
                          transform=fig.dpi_scale_trans, zorder=12)
    fig.add_artist(frame)
    head_y = (bot_in + ins_h + header_in * 0.52) / H
    fig.text((x + 0.087) / W, head_y, str(n), fontsize=5.5, fontweight="bold",
             color=col, ha="center", va="center",
             bbox=dict(boxstyle="circle,pad=0.18", facecolor=OCEAN,
                       edgecolor=col, linewidth=0.45))
    fig.text((x + 0.18) / W, head_y, name, fontsize=6.5, fontweight="bold",
             color=col, ha="left", va="center")
    zoom = (w / (bb[1] - bb[0])) / (W * main_w_frac / (MAIN[1] - MAIN[0]))
    fig.text((x + w - 0.065) / W, head_y, f"{zoom:.1f}×", fontsize=5.5,
             color="#8A969B", ha="right", va="center")
    zoom_factors.append(zoom)
    lo0, lo1, la0, la1 = bb
    for lab, px, py, tx, ty, ha in pts:
        axi.plot(px, py, "o", ms=1.15, color=LEAD, zorder=8)
        axi.annotate(lab, (px, py), xytext=(tx, ty), textcoords="data",
                     fontsize=5.5, color=INK, ha=ha, va="center",
                     arrowprops=dict(arrowstyle="-", color=LEAD, lw=0.3,
                                     shrinkA=2.0, shrinkB=1.5,
                                     connectionstyle="angle3,angleA=0,angleB=90" if lab == "Gabon" else "arc3"), zorder=9)
    # Small text-only offset in Asia keeps the endpoint label clear of both
    # the zero label and Christmas Island; bar geometry is unchanged.
    scale_bar(axi, lo0 + 1.9, la0 + 3.6, [0, 500, 1000],
              lat_ref=(la0 + la1) / 2, label_dx=-1.6 if n == 3 else 0.0)
    # Unambiguous source-to-zoom line. Each starts at an original box corner;
    # low contrast and a shallow curve keep the density pattern dominant.
    source = (lo0, la0) if n == 2 else (lo1, la0)
    target = ((x + 0.055) / W if n == 2 else (x + w - 0.055) / W,
              (bot_in + ins_h + header_in) / H)
    connector = ConnectionPatch(source, target, coordsA=ax.transData,
                                coordsB=fig.transFigure, axesA=ax,
                                arrowstyle="-", connectionstyle=f"arc3,rad={-0.14 if n == 2 else 0.10}",
                                lw=0.4, color=col, alpha=0.35, zorder=6,
                                clip_on=False)
    fig.add_artist(connector)
    insets.append(axi)
    panel_cards.append(card)
    x += w + gap

fig.savefig(ROOT / "Fig1a_density.pdf", dpi=600)
fig.savefig(ROOT / "Fig1a_density.png", dpi=600)
print(f"saved Fig1a_density.pdf/.png  {W*25.4:.3f} x {H*25.4:.3f} mm; "
      f"main {main_h_in:.3f} in, insets {ins_h:.3f} in; "
      f"zoom factors {', '.join(f'{z:.3f}' for z in zoom_factors)}")
