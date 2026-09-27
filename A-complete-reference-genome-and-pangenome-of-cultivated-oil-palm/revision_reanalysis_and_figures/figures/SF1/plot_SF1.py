#!/usr/bin/env python3
"""Supplementary Fig. 1 re-layout: Hi-C contact maps of ten haplotype assemblies.

The heat maps are the unmodified rasters rendered by HapHiC (500-kb bins, KR
normalization, --vmax_coef 4.0, white-red colour map, chromosome grid) and
extracted from each results/hic/<sample>/contact_map.pdf with pdfimages. Only the
layout, axes, labels and colour bars are redrawn. Colour-bar maxima are recovered
from the tick-label positions in the original PDFs.
"""
import re
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from PIL import Image

HERE = Path(__file__).resolve().parent
SRC = HERE / "SF1/hic_extract"
OUT = HERE / "out"
MM = 1 / 25.4
Image.MAX_IMAGE_PIXELS = None

mpl.rcParams.update({
    "font.family": "Arial", "font.size": 6, "axes.titlesize": 7, "xtick.labelsize": 5,
    "ytick.labelsize": 4.6, "axes.linewidth": 0.5, "xtick.major.width": 0.5,
    "ytick.major.width": 0.5, "xtick.major.size": 2, "ytick.major.size": 0,
    "axes.unicode_minus": True, "pdf.fonttype": 42,
    "mathtext.fontset": "custom", "mathtext.rm": "Arial", "mathtext.it": "Arial:italic",
})

# columns = materials, rows = haplotype 1 / haplotype 2 (names as used in the main text)
MATERIALS = [("BK", "TN"), ("dura", "TK"), ("pisifera", "NS"), ("nrly", "Nigerian"), ("MZ4", "E. oleifera")]
PAGE_H = 425.197  # pt, identical for all ten HapHiC PDFs


def words(sample):
    html = (SRC / sample / "text_bbox.html").read_text()
    return [(float(a), float(b), float(c), float(d), w) for a, b, c, d, w in
            re.findall(r'<word xMin="([\d.]+)" yMin="([\d.]+)" xMax="([\d.]+)" yMax="([\d.]+)">([^<]+)</word>', html)]


def placement(sample):
    out = {}
    for line in (SRC / sample / "placement.txt").read_text().split("\n"):
        if line.strip():
            k, a, b, c, d, e, f = line.split()
            out[k] = tuple(map(float, (a, b, c, d, e, f)))
    return out


def geometry(sample):
    pl = placement(sample)
    w = words(sample)
    # heat map extent in Mb from the x-axis tick labels (bottom row of integers)
    xt = [((x0 + x1) / 2, float(t)) for x0, y0, x1, y1, t in w if re.fullmatch(r"\d+", t) and y0 > 370]
    (xa, va), (xb, vb) = xt[0], xt[-1]
    pt_per_mb = (xb - xa) / (vb - va)
    img_w = pl["I1"][0]
    extent_mb = img_w / pt_per_mb
    # colour-bar maximum from the tick labels beside the colour bar
    ct = [(PAGE_H - (y0 + y1) / 2, float(t)) for x0, y0, x1, y1, t in w if re.fullmatch(r"0\.\d+", t)]
    ct.sort()
    (ya, ca), (yb, cb) = ct[0], ct[-1]
    val_per_pt = (cb - ca) / (yb - ya)
    cb_y0, cb_h = pl["I2"][5], pl["I2"][3]
    vmax = ca + (cb_y0 + cb_h - ya) * val_per_pt
    return extent_mb, vmax


def chroms(sample):
    rows = [l.split("\t") for l in (SRC / sample / f"{sample}.identity.agp").read_text().split("\n") if l.strip()]
    lens = np.array([int(r[2]) for r in rows]) / 1e6
    return [r[0] for r in rows], lens


fig = plt.figure(figsize=(180 * MM, 80 * MM))
left, right, top, bottom = 0.035, 0.99, 0.95, 0.085
ncol = len(MATERIALS)
cell_w = (right - left) / ncol
cell_h = (top - bottom) / 2
summary = []
for ci, (code, name) in enumerate(MATERIALS):
    for ri, hap in enumerate(("1", "2")):
        sample = f"{code}_hap{hap}"
        extent, vmax = geometry(sample)
        names, lens = chroms(sample)
        img = Image.open(SRC / sample / "img-000.png").convert("RGB")
        img = img.resize((1600, 1600), Image.LANCZOS)
        x0 = left + ci * cell_w
        y0 = top - (ri + 1) * cell_h
        ax_w = cell_w * 0.76
        ax_h = ax_w * fig.get_figwidth() / fig.get_figheight()
        ax = fig.add_axes([x0 + cell_w * 0.07, y0 + (cell_h - ax_h) * 0.6, ax_w, ax_h])
        ax.imshow(np.asarray(img), extent=[0, extent, 0, extent], origin="upper", interpolation="none")
        scale = extent / lens.sum()
        mids = (np.cumsum(lens) - lens / 2) * scale
        if ci == 0:
            ax.set_yticks(mids, [n.replace("chr", "") if (k < 13 or k == 15) else "" for k, n in enumerate(names)])
        else:
            ax.set_yticks([])
        ax.set_xticks(np.arange(0, extent, 500))
        ax.tick_params(axis="x", pad=1)
        ax.tick_params(axis="y", pad=0.8)
        for s in ax.spines.values():
            s.set_linewidth(0.5)
        label = f"{name}-Hap{hap}"
        ax.set_title(label, pad=2, fontstyle="normal")
        if name == "E. oleifera":
            ax.set_title("")
            ax.text(0.5, 1.02, r"$\mathit{E.\,oleifera}$" + f"-Hap{hap}", transform=ax.transAxes,
                    ha="center", va="bottom", fontsize=7)
        if ri == 1:
            ax.set_xlabel("Position (Mb)", labelpad=1)
        # slim colour bar
        pos = ax.get_position()
        cax = fig.add_axes([pos.x1 + 0.006, pos.y0, 0.0055, pos.height])
        cmap = mpl.colors.LinearSegmentedColormap.from_list("wr", ["white", "red"])
        cb = mpl.colorbar.ColorbarBase(cax, cmap=cmap, norm=mpl.colors.Normalize(0, vmax))
        cb.set_ticks([0, vmax])
        cb.set_ticklabels(["0", f"{vmax:.4f}"])
        cb.outline.set_linewidth(0.4)
        cax.tick_params(labelsize=4.4, length=1.5, width=0.4, pad=1)
        summary.append((label, sample, round(extent, 1), round(vmax, 5)))

fig.text(0.009, 0.52, "Chromosome", rotation=90, va="center", ha="center", fontsize=6)
for ext in ("pdf", "png"):
    fig.savefig(OUT / f"Supplementary_Fig_01.{ext}", dpi=600, facecolor="white")
with open(HERE / "SF1/SF1_panel_geometry.tsv", "w") as fh:
    fh.write("display_label\tsample\tmap_extent_Mb\tcolour_bar_max\n")
    for r in summary:
        fh.write("\t".join(map(str, r)) + "\n")
print("\n".join(map(str, summary)))
