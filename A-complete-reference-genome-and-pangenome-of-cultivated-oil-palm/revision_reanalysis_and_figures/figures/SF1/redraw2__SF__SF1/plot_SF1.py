#!/usr/bin/env python3
"""Supplementary Fig. 1 redraw (redraw2): Hi-C contact maps of the ten haplotype assemblies.

Data: the KR-normalized 500-kb contact matrices exactly as plotted by HapHiC 1.0.6 (HapHiC_plot.py
--bin_size 500 --normalization KR --vmax_coef 4.0), recovered bin-by-bin from the lossless
white->red rasters in results/hic/<sample>/contact_map.pdf (pixel value = min(KR/vmax, 1)).
vmax = the value written in each HapHiC_plot.log (trace/SF/sf1/vmax.tsv), bins/chromosomes from
<sample>.identity.agp (ceil(len / 500 kb) bins per chromosome; 3,463 bins for TN-Hap1 = HapHiC matrix).
Cross-check (${COMPUTE_HOST}, redraw2_sf1/kr_probe3.py): an independent KR balancing of the raw
contact_matrix.pkl correlates r = 0.983 with the recovered bin values, implied vmax ~0.0068 vs 0.00674 logged.
Presentation changes only: bin-resolution map (2 x 2 bin means = 1 Mb, ~1,700 px per map),
sequential white-to-dark-red scale (linear, 0..vmax unchanged), chromosome grid restored (HapHiC
--border_style grid), all 16 chromosome labels (01-16, staggered), colour bars with 0 / vmax/2 / vmax.
"""
import sys
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from PIL import Image

HERE = Path(__file__).resolve().parent
SP = HERE.parents[3]
SRC = SP / "redraw/SF1/hic_extract"
VMAX = {l.split("\t")[0]: float(l.split("\t")[1]) for l in (SP / "trace/SF/sf1/vmax.tsv").read_text().splitlines() if l.strip()}
sys.path.insert(0, str(SP / "fix/beautify/common"))
import palA  # noqa: E402

OUT = HERE / "out"
MM = 1 / 25.4
Image.MAX_IMAGE_PIXELS = None
DARK = "#242A30"
CMAP = mpl.colors.LinearSegmentedColormap.from_list(
    "hic", ["#FFFFFF", "#FAD2C4", "#F29078", "#DC4B3E", "#B01F2A", "#6E0F1A"])
GRID = "#CDD1D5"

mpl.rcParams.update({
    "font.family": "Arial", "font.size": 6, "axes.titlesize": 6.5, "xtick.labelsize": 5.5,
    "ytick.labelsize": 5, "axes.linewidth": 0.5, "xtick.major.width": 0.5, "ytick.major.width": 0.4,
    "xtick.major.size": 2, "ytick.major.size": 1.5, "axes.unicode_minus": True, "pdf.fonttype": 42,
    "mathtext.fontset": "custom", "mathtext.rm": "Arial", "mathtext.it": "Arial:italic",
})
MATERIALS = [("BK", "TN"), ("dura", "TK"), ("pisifera", "NS"), ("nrly", "Nigerian"), ("MZ4", "E. oleifera")]


def chrom_bins(sample):
    rows = [l.split("\t") for l in (SRC / sample / f"{sample}.identity.agp").read_text().splitlines() if l.strip()]
    nb = np.array([-(-int(r[2]) // 500_000) for r in rows])
    return [r[0] for r in rows], nb


def bin_matrix(sample, n):
    """HapHiC raster -> n x n matrix of KR/vmax in [0,1] (row 0 = first bin; HapHiC origin bottom_left)."""
    a = np.asarray(Image.open(SRC / sample / "img-000.png").convert("RGB")).astype(np.float32)
    assert (a[..., 0] >= 254).all()
    v = 1 - (a[..., 1] + a[..., 2]) / 2 / 255
    P = v.shape[0]
    c = ((np.arange(n) + 0.5) * P / n).astype(int)
    return v[::-1][np.ix_(c, c)]


def block2(m):
    n = m.shape[0] // 2 * 2
    b = m[:n, :n].reshape(n // 2, 2, n // 2, 2).mean(axis=(1, 3))
    return b


fig = plt.figure(figsize=(180 * MM, 80 * MM))
L, R, T, B = 0.064, 0.995, 0.975, 0.075
ncol = 5
cw = (R - L) / ncol
rh = (T - B) / 2
ax_w_mm = 28.0
ax_w = ax_w_mm / 180
ax_h = ax_w_mm / 80
geom = []
for ci, (code, name) in enumerate(MATERIALS):
    for ri, hap in enumerate(("1", "2")):
        sample = f"{code}_hap{hap}"
        names, nb = chrom_bins(sample)
        n = int(nb.sum())
        M = bin_matrix(sample, n)
        vmax = VMAX[sample]
        ext = n * 0.5
        x0 = L + ci * cw + 0.004
        y0 = T - (ri + 1) * rh + (rh - ax_h) * 0.45
        ax = fig.add_axes([x0, y0, ax_w, ax_h])
        ax.imshow(block2(M) * vmax, cmap=CMAP, vmin=0, vmax=vmax, origin="lower",
                  extent=[0, (n // 2 * 2) * 0.5, 0, (n // 2 * 2) * 0.5], interpolation="none")
        edges = np.cumsum(nb) * 0.5
        for e in edges[:-1]:
            ax.axhline(e, color=GRID, lw=0.25, zorder=3)
            ax.axvline(e, color=GRID, lw=0.25, zorder=3)
        ax.set_xlim(0, ext)
        ax.set_ylim(0, ext)
        ax.set_xticks(np.arange(0, ext, 500))
        ax.tick_params(axis="x", pad=1)
        mids = (np.concatenate([[0], edges[:-1]]) + edges) / 2
        if ci == 0:
            ax.set_yticks(mids, [s.replace("chr", "") for s in names])
            ax.tick_params(axis="y", pad=0.8, length=1.2)
            # staggered labels: even-numbered chromosomes sit one label-width further out on a longer tick
            for k, (tk, tl) in enumerate(zip(ax.yaxis.get_major_ticks(), ax.get_yticklabels())):
                if k % 2 == 1:
                    tk.tick1line.set_markersize(7.0)
                    tk.set_pad(7.0 + 0.8)
        else:
            ax.set_yticks([])
        for s in ax.spines.values():
            s.set_linewidth(0.5)
        lab = (r"$\mathit{E.\ oleifera}$" if name == "E. oleifera" else name) + f"-Hap{hap}"
        ax.set_title(lab, pad=4.5)
        if ri == 1:
            ax.set_xlabel("Position (Mb)", labelpad=1)
        pos = ax.get_position()
        cax = fig.add_axes([pos.x1 + 0.005, pos.y0, 0.0055, pos.height])
        cb = mpl.colorbar.ColorbarBase(cax, cmap=CMAP, norm=mpl.colors.Normalize(0, vmax))
        cb.set_ticks([0, vmax / 2, vmax])
        cb.set_ticklabels(["0", "", ""])
        cax.text(0.5, 1.012, f"{vmax:.4f}", transform=cax.transAxes, ha="center", va="bottom", fontsize=5)
        cb.outline.set_linewidth(0.4)
        cax.tick_params(labelsize=5, length=1.5, width=0.4, pad=1)
        geom.append((f"{name}-Hap{hap}", sample, n, round(ext, 1), vmax, f"{vmax:.4f}"))

fig.text(0.008, (T + B) / 2, "Chromosome", rotation=90, va="center", ha="center", fontsize=6)
for ext in ("pdf", "png"):
    fig.savefig(OUT / f"Supplementary_Fig_01.{ext}", dpi=600, facecolor="white")
with open(OUT / "SF1_panel_geometry.tsv", "w") as fh:
    fh.write("display_label\tsample\tbins_500kb\tmap_extent_Mb\tcolour_bar_max_HapHiC_log\tshown_as\n")
    for r in geom:
        fh.write("\t".join(map(str, r)) + "\n")
print("\n".join(map(str, geom)))
