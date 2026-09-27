#!/usr/bin/env python3
"""Figure 5h,i redrawn for the coverage-masked W = 4, P = 15 path and re-embedded into Figure5.pdf.

Base (read-only): scratchpad/deliver/Main_Figures_revised/Figure5.pdf (current, unmasked W = 4 version).
1. redact_hi.RECTS: the h and i content (text and vector art) is removed from a copy of the page; panel letters
   'h', 'i' and every other panel stay byte-for-byte in place (render outside the rectangles is pixel-identical,
   checked below).
2. The two panels are redrawn with matplotlib in page coordinates (1 unit = 1 pt; Arial, Type 42), using the
   geometry measured on the current page (row pitch, bar height, 0.754 pt/Mb in h; 0.5946-pt cells, 5.678-pt rows
   in i; font sizes 5.43-5.90 pt; triangle size; legend grid) and the author's visual rules
   (85_plot_african35_legacy_v5_gwas_vector.py: MOSAIC_COLORS shading by descending window count, HEATMAP_COLORS
   on ln(1 + load), rows of i sorted by ascending chr01B load, red boxes = selected segments with dSV/dSNP labels).
   Colours are kept as in the current figure (a later pass unifies the palette).
   New in i: cells in which the donor's alignment covers < 50% of the window (not eligible for selection) are
   drawn in grey (#C9C9C9), with a key next to the colour bar.
3. The overlay PDF is placed on the redacted page with show_pdf_page (vector).
Checks: the legacy style reproduces the current legend colours from the current W = 4 donor ranking.
Inputs (this directory): segments_mask.tsv, by_chrom_mask.tsv, per_chrom_capture_mask.tsv, position_marks_mask.tsv,
donor_windows_mask.tsv, zoom_matrix_chr01B.tsv (gen_path_mask.py).
Outputs: Figure5_candidate_mask.pdf, figure_edit_log_mask.tsv
"""
import sys
from pathlib import Path
import fitz
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, to_rgb
from matplotlib.patches import Rectangle, Polygon
import numpy as np
import pandas as pd

H = Path(__file__).resolve().parents[1]
S = Path(__file__).resolve().parents[3]
BASE = S / "deliver/Main_Figures_revised/Figure5.pdf"
OUT = Path(__file__).resolve().parent / "Figure5_candidate_hlabels.pdf"
sys.path.insert(0, str(Path(__file__).resolve().parent))
import redact_hi

MOSAIC_COLORS = ["#E9F3FA", "#D8E8F8", "#C8D8E8", "#AFCFE3", "#8FBAD6", "#68A8D8", "#5A94BF", "#7898C8",
                 "#627FB0", "#4D6F9A", "#3D638E", "#2D536F", "#213F58"]
HEATMAP_COLORS = ["#FBFCFE", "#EFF6FB", "#E2EFF7", "#CDE2F1", "#AFCFE6", "#88BBDC", "#68A8D8", "#3F83B0", "#245E7A"]
CHROM_BG = (0.953, 0.969, 0.98)
GOLD, GREY_T = (0.851, 0.643, 0.255), (0.722, 0.773, 0.8)
BOX_RED, LAB_RED = (0.784, 0.282, 0.282), "#B83F3E"
MASK_GREY = "#C9C9C9"
TXT = "#222222"
CHROMS = [f"chr{i:02d}B" for i in range(1, 17)]
ZC = "chr01B"
mcmap = LinearSegmentedColormap.from_list("mosaic", MOSAIC_COLORS)
hcmap = LinearSegmentedColormap.from_list("heat", HEATMAP_COLORS)


def disp(d):
    import re
    m = re.fullmatch(r"(dura|pisifera|nrly|bk)_hap([12])", d)
    return f"{dict(dura='TK', pisifera='NS', nrly='Nigerian', bk='TN')[m.group(1)]}-Hap{m.group(2)}" if m else d


def shades(order):
    n = len(order)
    return {d: mcmap(0.06 + 0.88 * (i / max(1, n - 1)))[:3] for i, d in enumerate(order)}


# ---- style check: the current legend colours are the legacy shading of the current W = 4 ranking
cur_rank = pd.read_csv(S / "fix/fig5hi_mask/donor_windows_mask.tsv", sep="\t").Donor_ID.tolist()   # current figure = masked 284-locus path
page0 = fitz.open(BASE)[0]
leg = sorted([x for x in page0.get_drawings() if x["seqno"] > 30000 and 630 < x["rect"].y0 < 682 and x["rect"].x1 < 210
              and abs(x["rect"].width - 4.33) < 0.02], key=lambda x: (round(x["rect"].y0, 1), x["rect"].x0))
sh0 = shades(cur_rank)
assert len(leg) == 35 and all(np.abs(np.array(x["fill"]) - np.array(sh0[d])).max() < 0.006 for x, d in zip(leg, cur_rank))

# ---- data
seg = pd.read_csv(H / "segments_mask.tsv", sep="\t")
bc = pd.read_csv(H / "by_chrom_mask.tsv", sep="\t").set_index("Chrom")
pc = pd.read_csv(H / "per_chrom_capture_mask.tsv", sep="\t").set_index("Chrom")
mk = pd.read_csv(H / "position_marks_mask.tsv", sep="\t")
rank = pd.read_csv(H / "donor_windows_mask.tsv", sep="\t").Donor_ID.tolist()
zm = pd.read_csv(H / "zoom_matrix_chr01B.tsv", sep="\t")
assert len(rank) == 35
sh = shades(rank)
chromlen_mb = {c: float(seg[seg.Chrom == c].Segment_End_0based.max()) / 1e6 for c in CHROMS}

PW, PH = page0.rect.width, page0.rect.height
plt.rcParams.update({"font.family": "Arial", "pdf.fonttype": 42, "ps.fonttype": 42})
fig = plt.figure(figsize=(PW / 72, PH / 72))
fig.patch.set_alpha(0)
ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(0, PW); ax.set_ylim(PH, 0); ax.axis("off"); ax.patch.set_alpha(0)


def text(x, y, s, size, color="k", ha="left", bold=False):
    ax.text(x, y, s, fontsize=size, color=color, ha=ha, va="baseline", fontweight="bold" if bold else "normal")


def tri(xc, ytop, w, h, color, edge, lw):
    ax.add_patch(Polygon([(xc - w / 2, ytop), (xc + w / 2, ytop), (xc, ytop + h)], closed=True, facecolor=color,
                         edgecolor=edge, linewidth=lw, joinstyle="miter"))


log = []
# ================= h =================
X0, KMB, BAR_Y0, BAR_H, PITCH = 50.76, 0.75420, 469.26, 5.05, 9.00653
for i, c in enumerate(CHROMS):
    y0 = BAR_Y0 + i * PITCH
    ax.add_patch(Rectangle((X0, y0), KMB * chromlen_mb[c], BAR_H, facecolor=CHROM_BG, edgecolor="none", linewidth=0))
    for r in seg[seg.Chrom == c].itertuples():
        xa = X0 + KMB * r.Segment_Start_0based / 1e6; xb = X0 + KMB * r.Segment_End_0based / 1e6
        ax.add_patch(Rectangle((xa, y0), xb - xa, BAR_H, facecolor=sh[r.Donor_ID], edgecolor="white", linewidth=0.068))
    for m in mk[mk.Chrom == c].itertuples():
        tri(X0 + KMB * m.Pos / 1e6, 467.18 + i * PITCH, 2.755, 2.2, GOLD if m.Mark else GREY_T, "white", 0.13)
    text(48.19, 473.548 + i * PITCH, c, 5.659, "k", ha="right")
    b = bc.loc[c]; p = pc.loc[c]
    lab = f"{int(b.Residual_DSV)} / {int(b.Residual_DSNP)}" + (f"   * {int(p.Captured)}/{int(p.Fav_Total)}" if int(p.Fav_Total) > 0 else "")
    text(190.285, 473.548 + i * PITCH, lab, 5.426)
    log.append(("5h label", c, lab))
text(185.29, 465.329, "dSV/dSNP", 5.709, "#727171", ha="right")
text(190.285, 465.329, "* captured/total", 5.426, "#727171")
ax.plot([25.87, 236.93], [612.3, 612.3], color="k", linewidth=0.456, solid_capstyle="projecting")
for k in range(8):
    x = X0 + 18.855 * k
    ax.plot([x, x], [612.3, 613.21], color="k", linewidth=0.292)
    text(x, 618.218, str(25 * k), 5.899, TXT, ha="center")
text(135.35, 624.972, "Chromosomal position (Mb)", 5.472, TXT, ha="center")
text(26.35, 632.789, "Donor haplotype (n = 35)", 5.472, "k", bold=True)
LX = [32.353, 75.568, 118.939, 162.406, 206.141]
LY = [639.716, 646.644, 653.571, 660.499, 667.427, 674.354, 681.692]
for n, d in enumerate(rank):
    r_, c_ = divmod(n, 5)
    ax.add_patch(Rectangle((LX[c_] - 5.99, LY[r_] - 3.0), 4.33, 2.49, facecolor=sh[d], edgecolor="none", linewidth=0))
    text(LX[c_], LY[r_], disp(d), 5.472)
for xc, col, lab, xt in ((28.17, GOLD, "fav. GWAS allele captured", 32.2), (115.23, GREY_T, "GWAS locus missed", 119.251)):
    tri(xc, 685.0 if col == GOLD else 685.66, 3.88, 3.36, col, col, 0.353)
    text(xt, 689.027, lab, 5.472)

# ================= i =================
zm["Donor"] = zm.Donor_ID.map(disp)
donors = sorted(zm.Donor_ID.unique())
nw = int(zm.Window_Index.max()) + 1
load = zm.pivot(index="Donor_ID", columns="Window_Index", values="Total_Load").reindex(index=donors).to_numpy()
elig = zm.pivot(index="Donor_ID", columns="Window_Index", values="Eligible").reindex(index=donors).to_numpy().astype(bool)
dsv = zm.pivot(index="Donor_ID", columns="Window_Index", values="DSV_Count").reindex(index=donors).to_numpy()
dsnp = zm.pivot(index="Donor_ID", columns="Window_Index", values="DSNP_Count").reindex(index=donors).to_numpy()
order = np.argsort(load.sum(1), kind="stable")                       # legacy: ascending chr01B load
rows = [donors[i] for i in order]; row_of = {d: r for r, d in enumerate(rows)}
HX0, CW, HY0, RH = 281.03, (491.524 - 281.03) / 354, 464.299, (663.043 - 464.299) / 35
assert nw == 354
vmax = float(np.log1p(load).max())
for r, d in enumerate(rows):
    j = donors.index(d)
    for w in range(nw):
        col = MASK_GREY if not elig[j, w] else hcmap(np.log1p(load[j, w]) / vmax)[:3]
        ax.add_patch(Rectangle((HX0 + w * CW, HY0 + r * RH), CW, RH, facecolor=col, edgecolor="none", linewidth=0,
                               antialiased=False))
    text(279.05, HY0 + r * RH + 4.689, disp(d), 5.598, TXT, ha="right")
zs = seg[seg.Chrom == ZC]
obst = []                                                               # (x0, x1, y0, y1): boxes, then labels
for s_ in zs.itertuples():
    r = row_of[s_.Donor_ID]
    obst.append((HX0 + s_.Start_Window_Index * CW, HX0 + s_.End_Window_Index_Exclusive * CW, HY0 + r * RH, HY0 + (r + 1) * RH))
for s_ in zs.itertuples():
    r = row_of[s_.Donor_ID]; j = donors.index(s_.Donor_ID)
    xa = HX0 + s_.Start_Window_Index * CW; xb = HX0 + s_.End_Window_Index_Exclusive * CW
    ax.add_patch(Rectangle((xa, HY0 + r * RH), xb - xa, RH, fill=False, edgecolor=BOX_RED, linewidth=0.456))
    lab = f"{int(dsv[j, s_.Start_Window_Index:s_.End_Window_Index_Exclusive].sum())}/{int(dsnp[j, s_.Start_Window_Index:s_.End_Window_Index_Exclusive].sum())}"
    wlab = 0.556 * 5.598 * len(lab)                                   # Arial digits/slash ~0.556 em
    base = HY0 + r * RH + 4.689
    cands = []
    for dy in (0, -RH, RH, -2 * RH, 2 * RH):
        cands += [(xb + 0.6, base + dy), (xa - 0.6 - wlab, base + dy)]
    def free(x, y):
        if x < HX0 + 0.3 or x + wlab > 512 or y - 4.2 < HY0 - 0.5 or y > HY0 + 35 * RH + 0.5:
            return False
        return all(x + wlab < o0 - 0.2 or x > o1 + 0.2 or y + 0.9 < p0 + 0.05 or y - 4.1 > p1 - 0.05
                   for o0, o1, p0, p1 in obst)
    x, base = next((c for c in cands if free(*c)), cands[0])
    obst.append((x, x + wlab, base - 4.1, base + 0.9))
    # labels removed (fig35_minor, 2026-09-26): per-window counts are in Source Data
    log.append(("5i box", s_.Donor_ID, s_.Start_Window_Index, s_.End_Window_Index_Exclusive, lab, round(x, 2), round(base, 2)))
for m in mk[mk.Chrom == ZC].itertuples():
    tri(HX0 + CW * 2 * m.Pos / 1e6, 460.98, 2.785, 2.30, GOLD if m.Mark else GREY_T, "white", 0.1295)
for k, v in enumerate((0, 50, 100, 150)):
    x = HX0 + CW * 2 * v
    ax.plot([x, x], [663.03, 664.56], color="k", linewidth=0.292)
    text(x, 669.95, str(v), 5.598, TXT, ha="center")
text(326.0, 675.443, "chr01B position (Mb)", 5.472, TXT, ha="center")
CB0, CBK = 277.05, 16.3115                                              # colour bar: 16.31 pt per ln unit
nb = 300
for q in range(nb):
    v0 = vmax * q / nb
    ax.add_patch(Rectangle((CB0 + CBK * v0, 678.26), CBK * vmax / nb + 0.01, 5.17, facecolor=hcmap(v0 / vmax)[:3],
                           edgecolor="none", linewidth=0, antialiased=False))
for k in range(int(vmax) + 1):
    x = CB0 + CBK * k
    ax.plot([x, x], [678.26, 679.29], color="white", linewidth=0.342)
    ax.plot([x, x], [682.39, 683.43], color="white", linewidth=0.342)
    text(x, 688.581, str(k), 5.472, TXT, ha="center")
cb_end = CB0 + CBK * vmax
text(cb_end + 2.6, 682.87, "ln(1 + dSV + dSNP)", 5.472)
kx = cb_end + 2.6 + 52.5
ax.add_patch(Rectangle((kx, 678.26), 4.0, 5.17, facecolor=MASK_GREY, edgecolor="none", linewidth=0))
text(kx + 5.5, 682.87, "< 50% aligned (not eligible)", 5.472)

tmp = H / "_hi_overlay.pdf"
fig.savefig(tmp, transparent=True)
plt.close(fig)

# ---- compose
doc = fitz.open(BASE); page = doc[0]
redact_hi.redact(page)
ov = fitz.open(tmp)
page.show_pdf_page(page.rect, ov, 0, overlay=True)
doc.save(OUT, garbage=3, deflate=True)
tmp.unlink()

# ---- check: outside the redacted rectangles the render is unchanged
a = fitz.open(BASE)[0].get_pixmap(dpi=200); b = fitz.open(OUT)[0].get_pixmap(dpi=200)
A_ = np.frombuffer(a.samples, np.uint8).reshape(a.h, a.w, a.n); B_ = np.frombuffer(b.samples, np.uint8).reshape(b.h, b.w, b.n)
diff = (A_ != B_).any(2); s = 200 / 72; msk = np.ones_like(diff)
for rr in redact_hi.RECTS:
    msk[int(rr.y0 * s):int(rr.y1 * s) + 1, int(rr.x0 * s):int(rr.x1 * s) + 1] = False
print("pixels changed outside h/i rectangles:", int((diff & msk).sum()))
pd.DataFrame(log).to_csv(Path(__file__).resolve().parent / "figure_edit_log_hlabels.tsv", sep="\t", index=False, header=False)
print("5i vmax", round(vmax, 3), "rows", [disp(d) for d in rows][:6], "...")
