#!/usr/bin/env python3
"""Extended Data Fig. 1 re-layout (180 x 170 mm).

a  study design: author vector PDF 油棕种质.pdf (22 Sep) with text corrections (fix_flowchart_pdf.py),
   title removed; rendered at 1200 dpi
b  SubPhaser ancestry circos: cropped unchanged from the current ED Fig. 1 raster (HB-06 image1.png)
c  HiFi/ONT coverage: tracks rasterized at 600 dpi from gaps_readscov.pdf (vector); chromosome labels,
   legend ('Fillled' typo fixed) and per-column scale bars redrawn. The source scales each column to
   its own longest chromosome (chr01A, 192.51 Mb; chr01B, 176.98 Mb), so each column gets its own bar.
d  IGV gap-closure evidence: 2 of the 4 loci in gapsfilling.pdf (embedded screenshots extracted unchanged)
e  Pore-C contact map of FL: heat-map pixels cropped unchanged from the current ED Fig. 1 raster
   (no Pore-C source file found under 00_ms); axes, labels and colour bar redrawn from measured geometry
"""
from pathlib import Path

import fitz
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from PIL import Image

HERE = Path(__file__).resolve().parent
S = HERE / "src"
MM = 1 / 25.4
Image.MAX_IMAGE_PIXELS = None
mpl.rcParams.update({"font.family": "Arial", "font.size": 6, "axes.linewidth": 0.5, "xtick.labelsize": 5.5,
                     "ytick.labelsize": 5.5, "xtick.major.width": 0.5, "ytick.major.width": 0.5,
                     "xtick.major.size": 2, "ytick.major.size": 2, "axes.unicode_minus": True, "pdf.fonttype": 42})
W, H = 180, 170
fig = plt.figure(figsize=(W * MM, H * MM))


def ax_mm(x, y_top, w, h, **kw):
    return fig.add_axes([x / W, 1 - (y_top + h) / H, w / W, h / H], **kw)


def letter(x, y_top, s):
    fig.text(x / W, 1 - y_top / H, s, fontsize=9, fontweight="bold", va="top", ha="left")


def show(ax, arr):
    ax.imshow(arr, interpolation="lanczos")
    ax.axis("off")


cur = Image.open(HERE.parent.parent / "x_HB-06/word/media/image1.png").convert("RGB")
CS = cur.size[0] / 1270  # preview->full scale used when the crops were located

# ---------------- a ----------------
_pa = fitz.open(HERE / "ED1a_study_design.pdf")[0].get_pixmap(dpi=1200)
a_img = Image.frombytes("RGB", (_pa.width, _pa.height), _pa.samples)
aw = 92; ah = aw * a_img.size[1] / a_img.size[0]
ax = ax_mm(3, 4, aw, ah); show(ax, np.asarray(a_img)); letter(0, 1.5, "a")

# ---------------- b ----------------
reg = np.asarray(cur.crop((int(750 * CS), int(15 * CS), int(1268 * CS), int(550 * CS)))).astype(int)
nz = np.where((255 - reg).sum(2) > 40)
b_img = Image.fromarray(reg[nz[0].min():nz[0].max() + 1, nz[1].min():nz[1].max() + 1].astype(np.uint8))
bh = ah; bw = bh * b_img.size[0] / b_img.size[1]
bx = 3 + aw + 4
ax = ax_mm(bx + (W - bx - bw) / 2, 4, bw, bh); show(ax, np.asarray(b_img)); letter(bx - 1, 1.5, "b")
b_img.save(HERE / "ED1b_circos_crop.png")

# ---------------- c ----------------
doc = fitz.open(S / "gaps_readscov.pdf"); pg = doc[0]
lens = pd.read_csv(S / "FL_chr_lengths.tsv", sep="\t").set_index("Formal_chromosome_ID").Length_bp / 1e6
row_top = 4 + ah + 6
c_w, c_h = 106, 48
fig_c = []
for k, (side, x0pt, x1pt, ref) in enumerate([("A", 50, 606, "chr01A"), ("B", 652, 1210, "chr01B")]):
    pix = pg.get_pixmap(dpi=600, clip=fitz.Rect(x0pt, 0, x1pt, 680))
    arr = np.frombuffer(pix.samples, np.uint8).reshape(pix.height, pix.width, pix.n)[:, :, :3].copy()
    if side == "B":  # white-out the legend box (does not overlap any track)
        sx = 600 / 72
        arr[int(590 * sx):int(668 * sx), int((930 - x0pt) * sx):int((1102 - x0pt) * sx)] = 255
    # grey chromosome bars -> row centres and bar lengths
    grey = (np.abs(arr.astype(int) - 200).sum(2) < 40)
    rows = np.where(grey.sum(1) > 0.10 * grey.shape[1])[0]
    groups = np.split(rows, np.where(np.diff(rows) > 5)[0] + 1)
    centres = [g.mean() for g in groups if len(g) > 3]
    bar_len = grey[int(centres[0])].nonzero()[0]
    px_per_mb = (bar_len.max() - bar_len.min()) / lens[ref]
    # orange filled-gap marks
    orange = (arr[:, :, 0] > 170) & (arr[:, :, 1] > 100) & (arr[:, :, 1] < 150) & (arr[:, :, 2] < 110)
    gap_rows = [i for i, c in enumerate(centres) if orange[max(0, int(c) - 25):int(c) + 25].any()]
    fig_c.append((side, len(centres), round(px_per_mb, 3), [f"chr{r + 1:02d}{side}" for r in gap_rows]))
    colw = (c_w - 10) / 2
    ax = ax_mm(3 + 8 + k * (colw + 6), row_top, colw, c_h - 8)
    ax.imshow(arr, interpolation="lanczos", aspect="auto")
    ax.set_xlim(0, arr.shape[1]); ax.set_ylim(arr.shape[0], 0)
    ax.set_yticks(centres, [f"{i + 1:02d}{side}" for i in range(len(centres))])
    ax.set_xticks([]); ax.tick_params(axis="y", length=0, pad=1, labelsize=5.2)
    for s in ax.spines.values():
        s.set_visible(False)
    # scale bar 50 Mb
    L = 50 * px_per_mb
    ax.plot([arr.shape[1] - L - 20, arr.shape[1] - 20], [arr.shape[0] + 60] * 2, color="black", lw=0.8, clip_on=False)
    ax.text(arr.shape[1] - L / 2 - 20, arr.shape[0] + 110, "50 Mb", ha="center", va="top", fontsize=5.2)
    ax.set_title(f"FL-Hap{'1' if side == 'A' else '2'} (chr01{side}–chr16{side})", fontsize=6, pad=2)
letter(0, row_top - 3, "c")
hand = [mpl.patches.Patch(color="#F28080", label="HiFi read coverage"),
        mpl.patches.Patch(color="#4DBFAD", label="ONT read coverage"),
        mpl.patches.Patch(color="#C07A45", label="Filled gap")]
fig.legend(handles=hand, loc="lower left", bbox_to_anchor=(11 / W, 1 - (row_top + c_h + 1.5) / H), ncol=3,
           frameon=False, fontsize=5.5, handlelength=1.0, handleheight=0.7, columnspacing=1.0)

# ---------------- e ----------------
ex = 3 + c_w + 5
e_side = 41
e_crop = cur.crop((390 + 148, 6118 + 115, 390 + 2801, 6118 + 2768))
extent = (2805 - 144) / 370.71 * 500  # Mb, from x-axis tick spacing (500 Mb = 370.7 px)
vmax = 0.008 + (786 - 737) / 170 * 0.001  # colour-bar top from tick spacing (0.001 = 170 px)
ax = ax_mm(ex + 8, row_top + 1, e_side, e_side)
ax.imshow(np.asarray(e_crop), extent=[0, extent, 0, extent], interpolation="lanczos")
# chromosome label centres measured on the source y axis (32 ticks, top chr16B -> bottom chr01A)
ymid_px = [129, 165, 209, 261, 314, 369, 427, 489, 553, 615, 689, 777, 853, 916, 983, 1061, 1155, 1260, 1361,
           1461, 1560, 1662, 1762, 1861, 1961, 2063, 2159, 2253, 2348, 2446, 2563, 2700][::-1]
y_mb = [(2772 - y) / 370.71 * 500 for y in ymid_px]  # bottom frame at px 2772
pair_mid = [(y_mb[2 * i] + y_mb[2 * i + 1]) / 2 for i in range(16)]
ax.set_yticks(pair_mid, [f"{i + 1:02d}" for i in range(16)]); ax.tick_params(axis="y", labelsize=4.8, pad=1, length=1.5)
ax.set_xticks(np.arange(0, extent, 1000)); ax.set_xlabel("Position (Mb)", labelpad=1)
ax.set_ylabel("Chromosome (A and B homologues)", labelpad=1)
ax.set_title("FL Pore-C contacts (500-kb bins)", fontsize=6, pad=2)
cax = ax_mm(ex + 8 + e_side + 1.5, row_top + 9, 1.6, e_side - 16)
cmap = mpl.colors.LinearSegmentedColormap.from_list("wr", ["white", "#E3120B"])
cb = mpl.colorbar.ColorbarBase(cax, cmap=cmap, norm=mpl.colors.Normalize(0, vmax))
cb.set_ticks([0, 0.004, 0.008]); cb.outline.set_linewidth(0.4)
cax.tick_params(labelsize=4.8, length=1.5, width=0.4, pad=1)
cax.set_ylabel("KR-normalized contacts", fontsize=5, labelpad=2)
letter(ex, row_top - 3, "e")

# ---------------- d ----------------
d_top = row_top + c_h + 5
d_h = H - d_top - 1.5
# Both IGV loci are terminal extensions (not internal gaps). Final coordinates from minimap2 mapping of
# the draft-window flanks to the final FL assembly (see Source Data ED1cd_gap_closures).
loci = [("hap1_Chr10", (55, 56), "FL-Hap1 chr12A, left-end extension (~60 kb)"),
        ("hap2_Chr29", (60, 59), "FL-Hap2 chr07B, right-end extension (~1.10 Mb)")]
dw = (W - 6 - 4) / 2
for k, (name, xrefs, title) in enumerate(loci):
    ims = [Image.open(__import__("io").BytesIO(doc_d.extract_image(x)["image"])).convert("RGB")
           for doc_d in [fitz.open(S / "gapsfilling.pdf")] for x in xrefs]
    both = Image.new("RGB", (ims[0].size[0] + ims[1].size[0] + 40, max(i.size[1] for i in ims)), "white")
    both.paste(ims[0], (0, 0)); both.paste(ims[1], (ims[0].size[0] + 40, 0))
    both.save(HERE / f"ED1d_{name}.png")
    hh = min(d_h - 3, dw * both.size[1] / both.size[0]); ww = hh * both.size[0] / both.size[1]
    ax = ax_mm(3 + k * (dw + 4) + (dw - ww) / 2, d_top + 3, ww, hh); show(ax, np.asarray(both))
    fig.text((3 + k * (dw + 4)) / W, 1 - (d_top + 2.6) / H, title, fontsize=6, va="bottom")
letter(0, d_top - 1, "d")

for ext in ("pdf", "png"):
    fig.savefig(HERE / f"Extended_Data_Fig_01.{ext}", dpi=600, facecolor="white")
pd.DataFrame(fig_c, columns=["column", "chromosomes", "px_per_Mb_at_600dpi", "rows_with_filled_gap"]).to_csv(
    HERE / "ED1c_layout_check.tsv", sep="\t", index=False)
print(fig_c); print("a h mm", round(ah, 1), "extent Mb", round(extent, 1), "vmax", round(vmax, 5))
