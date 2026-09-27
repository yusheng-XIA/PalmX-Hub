#!/usr/bin/env python3
"""Extended Data Fig. 1 re-layout (180 x 170 mm).

a  study design: author vector PDF 油棕种质.pdf (22 Sep) with text corrections (fix_flowchart_pdf.py),
   title removed; rendered at 1200 dpi
b  (redraw2) SubPhaser circos redrawn from its own Circos input tables
   (00_ms/05_MS/0918_revision/completed_genome_v2.sg.k15_q200_f2.circos_renamed/data/*.txt, copy in ../../src/ED1b)
c  (redraw2) HiFi/ONT coverage redrawn as vector from the per-bin rectangles of the author's vector plot
   gaps_readscov.pdf (100-kb bins, relative depth in the author's plot units); both haplotypes on one scale
d  read-alignment panel from fix/ed1d_redraw/ed1d_panel.py (draw_d), unchanged
e  (redraw2) FL Pore-C contact map redrawn from the HapHiC raw contact matrix (contact_matrix.pkl, 500-kb bins,
   final assembly completed_genome_v2) after KR normalisation with HapHiC's own code (../../src/ED1e/ed1e_kr.py);
   vmax = 4 x median intra-chromosomal value = 0.0082895 (identical to the HapHiC log of the published map)
"""
from pathlib import Path

import fitz
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from PIL import Image

WORK = Path(__file__).resolve().parent               # redraw2/ED/work/ED1 (this script)
FIX = WORK.parents[3]                                # scratchpad/fix
OUT = WORK.parents[1] / "out"                        # redraw2/ED/out
HERE = FIX / "beautify/work/ED1"                     # panel a / c sources (read only)
S = HERE / "src"
SRC2 = WORK.parents[1] / "src"                       # ED1b SubPhaser circos tables, ED1e KR matrix
D1D = FIX / "ed1d_redraw"                            # panel d is drawn by the ED1d module (unchanged)
import sys; sys.path.insert(0, str(D1D)); from ed1d_panel import draw_d, write_source_data  # noqa: E402
sys.path.insert(0, str(WORK.parents[1] / "common")); import palA  # noqa: E402
from matplotlib.path import Path as MPath  # noqa: E402
from matplotlib.patches import PathPatch, Polygon, Rectangle  # noqa: E402
(WORK / "data").mkdir(exist_ok=True)
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
    fig.text(x / W, 1 - y_top / H, s, fontsize=8, fontweight="bold", va="top", ha="left")


def show(ax, arr):
    ax.imshow(arr, interpolation="lanczos")
    ax.axis("off")



# ---------------- a ----------------
# beautify: panel a is the Arial-converted vector flowchart (ED1a_to_arial.py), placed as vector after saving
A_PDF = HERE / "ED1a_study_design_arial.pdf"
_ra = fitz.open(A_PDF)[0].rect
aw = 92; ah = aw * _ra.height / _ra.width
A_RECT_MM = (3, 4, aw, ah)
letter(0, 1.5, "a")

# ---------------- b (redraw2: SubPhaser circos from its Circos input tables) ----------------
SG1, SG2, LTR_OTHER = "#E0B04A", "#4F8FC6", "#B9BEC4"   # SG1 gold / SG2 blue as in the SubPhaser original, softened
fl_len = pd.read_csv(S / "FL_chr_lengths.tsv", sep="\t").set_index("Formal_chromosome_ID").Length_bp
kar = pd.read_csv(SRC2 / "ED1b/data/genome_karyotype.txt", sep=r"\s+", header=None,
                  names=["t", "p", "id", "lab", "s", "e", "col"])
kar["formal"] = kar.id.str.replace(r"chr(\d\d)_hap1", r"chr\1A", regex=True).str.replace(r"chr(\d\d)_hap2", r"chr\1B",
                                                                                         regex=True)
assert (kar.e.values == fl_len.loc[kar.formal].values).all()   # hap1 = A (FL-Hap1), hap2 = B (FL-Hap2), exact lengths


def rd_track(f):
    d = pd.read_csv(SRC2 / "ED1b/data" / f, sep=r"\s+", header=None)
    v = d[3].astype(str).str.split(",", expand=True).astype(float)
    d = pd.concat([d.iloc[:, :3], v], axis=1)
    d.columns = ["chr", "s", "e"] + [f"v{i}" for i in range(v.shape[1])]
    return d


enr, rat = rd_track("sg_enrich.txt"), rd_track("sg_ratio.txt")
sg1, sg2, ltr = rd_track("subgenome.SG1.txt"), rd_track("subgenome.SG2.txt"), rd_track("ltr_density.txt")
links = pd.read_csv(SRC2 / "ED1b/data/block_link.txt", sep=r"\s+", header=None,
                    names=["c1", "s1", "e1", "c2", "s2", "e2", "col"])
GAP_DEG, TOP_GAP = 0.55, 7.0                       # between chromosomes; wider gap at 12 o'clock for ring numbers
tot = kar.e.sum()
deg_per_bp = (360 - TOP_GAP - GAP_DEG * (len(kar) - 1)) / tot
a_start, cum = {}, 90 - TOP_GAP / 2               # clockwise from the top
for _, r in kar.iterrows():
    a_start[r.id] = cum
    cum -= r.e * deg_per_bp + GAP_DEG


def ang(c, pos):
    return np.deg2rad(a_start[c] - pos * deg_per_bp)


def xy(t, r):
    return r * np.cos(t), r * np.sin(t)


def arc(c, s, e, r, n=None):
    n = n or max(2, int((e - s) / 1e6) + 2)
    t = np.linspace(ang(c, s), ang(c, e), n)
    return np.c_[r * np.cos(t), r * np.sin(t)]


def band(ax, c, s, e, r0, r1, col, **kw):
    ax.add_patch(Polygon(np.r_[arc(c, s, e, r1), arc(c, s, e, r0)[::-1]], closed=True, facecolor=col, lw=0, **kw))


def hist(ax, d, c, r0, r1, lo, hi, col):
    """filled step histogram between value levels lo..hi (fractions of the track height) for one chromosome"""
    q = d[d.chr == c].sort_values("s")
    if q.empty:
        return
    t = np.ravel(np.c_[ang(c, q.s.values), ang(c, q.e.values)])
    rlo = r0 + np.repeat(lo, 2) * (r1 - r0)
    rhi = r0 + np.repeat(hi, 2) * (r1 - r0)
    pts = np.r_[np.c_[rhi * np.cos(t), rhi * np.sin(t)], np.c_[rlo * np.cos(t), rlo * np.sin(t)][::-1]]
    ax.add_patch(Polygon(pts, closed=True, facecolor=col, lw=0))


bw = W - (3 + aw + 4) - 1
bside = min(bw, ah)
bx = 3 + aw + 4 + (bw - bside) / 2
ax = ax_mm(bx, 3.5, bside, bside)
ax.set_xlim(-1.13, 1.13); ax.set_ylim(-1.13, 1.13); ax.set_aspect("equal"); ax.axis("off")
RINGS = [(0.88, 0.99), (0.76, 0.87), (0.64, 0.75), (0.52, 0.63), (0.39, 0.50)]   # Circos r0/r1 of the source
for c in kar.id:
    # ring 1: significant enrichment of SG-specific k-mers (1-Mb windows)
    q = enr[enr.chr == c].sort_values("s").reset_index(drop=True)
    lab = np.where(q.v0 > 0, 1, np.where(q.v1 > 0, 2, 0))
    run = np.r_[0, np.cumsum(lab[1:] != lab[:-1])]           # merge consecutive windows with the same call
    for _, g in q.groupby(run):
        L_ = lab[g.index[0]]
        if L_:
            band(ax, c, g.s.min(), g.e.max(), *RINGS[0], SG1 if L_ == 1 else SG2)
    # ring 2: normalised proportion of SG1/SG2-specific k-mers (stacked, SG1 inside)
    q = rat[rat.chr == c].sort_values("s")
    hist(ax, rat, c, *RINGS[1], np.zeros(len(q)), q.v0.values, SG1)
    hist(ax, rat, c, *RINGS[1], q.v0.values, (q.v0 + q.v1).values, SG2)
    # rings 3-4: counts of SG1- and SG2-specific k-mers (each scaled to its own maximum, as in Circos)
    for d_, rr, colr in ((sg1, RINGS[2], SG1), (sg2, RINGS[3], SG2)):
        q = d_[d_.chr == c].sort_values("s")
        hist(ax, d_, c, *rr, np.zeros(len(q)), (q.v0 / d_.v0.max()).values, colr)
        ax.plot(*arc(c, 0, kar.set_index("id").e[c], rr[0]).T, color="#D5D8DC", lw=0.25)
    # ring 5: LTR-RT density, stacked SG1-specific / SG2-specific / other (scaled to the maximum stacked total)
    q = ltr[ltr.chr == c].sort_values("s")
    m = (ltr.v0 + ltr.v1 + ltr.v2).max()
    c0, c1, c2 = q.v0.values / m, (q.v0 + q.v1).values / m, (q.v0 + q.v1 + q.v2).values / m
    hist(ax, ltr, c, *RINGS[4], np.zeros(len(q)), c0, SG1)
    hist(ax, ltr, c, *RINGS[4], c0, c1, SG2)
    hist(ax, ltr, c, *RINGS[4], c1, c2, LTR_OTHER)
# links: homologous blocks (ribbons through the centre, Circos bezier_radius = 0), coloured by chromosome
LINK_R = 0.37
import colorsys  # noqa: E402
LCOL = {f"color=chr{i + 1}": colorsys.hls_to_rgb(((i * 7) % 16) / 16, 0.62, 0.42) for i in range(16)}
for _, l in links.drop_duplicates().iterrows():
    p1, p2 = arc(l.c1, l.s1, l.e1, LINK_R, 6), arc(l.c2, l.s2, l.e2, LINK_R, 6)
    verts = np.r_[p1, [[0, 0]], p2[::-1][:1], p2[::-1][1:], [[0, 0]], p1[:1]]
    codes = [MPath.MOVETO] + [MPath.LINETO] * (len(p1) - 1) + [MPath.CURVE3, MPath.CURVE3] + \
            [MPath.LINETO] * (len(p2) - 1) + [MPath.CURVE3, MPath.CURVE3]
    ax.add_patch(PathPatch(MPath(verts, codes), facecolor=LCOL[l.col], edgecolor="none", alpha=0.55, lw=0))
for _, r in kar.iterrows():   # chromosome labels
    t = ang(r.id, r.e / 2)
    deg = np.rad2deg(t)
    rot = deg - 90 if np.sin(t) >= 0 else deg + 90
    ax.text(*xy(t, 1.055), r.formal.replace("chr", ""), rotation=rot, rotation_mode="anchor", ha="center",
            va="center", fontsize=5.2)
for i, (r0, r1) in enumerate(RINGS):   # ring numbers in the top gap
    ax.text(0, (r0 + r1) / 2, str(i + 1), ha="center", va="center", fontsize=5, color="#4A4A4A")
ax.legend(handles=[mpl.patches.Patch(color=SG1, label="SG1"), mpl.patches.Patch(color=SG2, label="SG2"),
                   mpl.patches.Patch(color=LTR_OTHER, label="Other LTR-RTs (ring 5)")],
          loc="lower left", bbox_to_anchor=(-0.03, -0.02), frameon=False, fontsize=5.5, handlelength=0.9,
          handleheight=0.8, handletextpad=0.4, labelspacing=0.3, borderaxespad=0)
letter(3 + aw + 3, 1.5, "b")
pd.DataFrame({"circos_chromosome": kar.id, "formal_chromosome": kar.formal, "length_bp": kar.e,
              "SubPhaser_karyotype_subgenome": kar.col}).to_csv(WORK / "data/ED1b_chromosome_map.tsv", sep="\t",
                                                                 index=False)

# ---------------- c (redraw2: vector coverage from the author's per-bin rectangles) ----------------
doc = fitz.open(S / "gaps_readscov.pdf"); pg = doc[0]
HIFI, ONT, GREYF, GAPF = (0.94, 0.46, 0.48), (0.29, 0.75, 0.67), (0.78, 0.78, 0.78), (0.75, 0.49, 0.29)
rows = []
for dr in pg.get_drawings():
    f = dr.get("fill")
    if f is None:
        continue
    r = dr["rect"]
    rows.append((tuple(round(x, 2) for x in f), r.x0, r.y0, r.x1, r.y1))
rc = pd.DataFrame(rows, columns=["fill", "x0", "y0", "x1", "y1"])
bars = rc[rc.fill == GREYF].copy()
bars["side"] = np.where(bars.x0 < 605, "A", "B")
bars = bars.sort_values(["side", "y0"])
bars["chrom"] = [f"chr{i % 16 + 1:02d}{s}" for i, s in enumerate(bars.side)]
bars["yc"] = (bars.y0 + bars.y1) / 2
bars["len_bp"] = fl_len.loc[bars.chrom].values
bars["bp_per_pt"] = bars.len_bp / (bars.x1 - bars.x0)
assert bars.groupby("side").bp_per_pt.agg(lambda v: v.max() / v.min() - 1).max() < 0.004   # one scale per column
cov = []
for key, trk in ((HIFI, "HiFi"), (ONT, "ONT")):
    t = rc[(rc.fill == key) & ((rc.x1 - rc.x0) < 2)].copy()           # drop the legend swatch
    for side in "AB":
        bs = bars[bars.side == side]
        tt = t[(t.x0 < 605) == (side == "A")].copy()
        idx = np.abs(((tt.y0 + tt.y1) / 2).values[:, None] - bs.yc.values[None, :]).argmin(1)
        b = bs.iloc[idx].reset_index(drop=True)
        cov.append(pd.DataFrame({"track": trk, "chrom": b.chrom, "start_bp": ((tt.x0.values - b.x0) * b.bp_per_pt),
                                 "end_bp": ((tt.x1.values - b.x0) * b.bp_per_pt),
                                 "depth_plot_units": (tt.y1 - tt.y0).values}))
cov = pd.concat(cov, ignore_index=True)
cov["start_bp"] = (cov.start_bp / 1e5).round().astype(int) * 100000     # 100-kb bins
cov["end_bp"] = np.minimum((cov.end_bp / 1e5).round().astype(int) * 100000, fl_len.loc[cov.chrom].values)
gaps = rc[(rc.fill == GAPF) & ((rc.y1 - rc.y0) < 5)].copy()
gb = [bars.iloc[np.abs(bars.yc.values - (g.y0 + g.y1) / 2).argmin()] for _, g in gaps.iterrows()]
gaps["chrom"] = [b.chrom for b in gb]
gaps["start_bp"] = [(g.x0 - b.x0) * b.bp_per_pt for (_, g), b in zip(gaps.iterrows(), gb)]
gaps["end_bp"] = [(g.x1 - b.x0) * b.bp_per_pt for (_, g), b in zip(gaps.iterrows(), gb)]
gaps = gaps.sort_values("chrom")[["chrom", "start_bp", "end_bp"]]
cov.to_csv(WORK / "data/ED1c_HiFi_ONT_coverage_100kb.tsv", sep="\t", index=False)
gaps.round(0).astype({"start_bp": int, "end_bp": int}).to_csv(WORK / "data/ED1c_filled_gap_marks.tsv", sep="\t",
                                                              index=False)
row_top = 4 + ah + 6
c_w, c_h = 106, 48
HIFI_C, ONT_C, GAP_C = "#D9676B", "#3AA58F", "#8A5A34"
YMAX = 16.0                                    # plot-unit cap of the source drawing (max HiFi bar = 16.0)
colw = (c_w - 10) / 2
xmax = fl_len.max() / 1e6 + 1
for k, side in enumerate("AB"):
    ax = ax_mm(3 + 8 + k * (colw + 6), row_top, colw, c_h - 8)
    for i in range(16):
        cname = f"chr{i + 1:02d}{side}"
        y = 15 - i
        L = fl_len[cname] / 1e6
        ax.add_patch(Rectangle((0, y - 0.07), L, 0.14, facecolor="#D0D3D6", edgecolor="none", zorder=2))
        for trk, col, sgn in (("HiFi", HIFI_C, 1), ("ONT", ONT_C, -1)):
            q = cov[(cov.track == trk) & (cov.chrom == cname)].sort_values("start_bp")
            xs = np.ravel(np.c_[q.start_bp, q.end_bp]) / 1e6
            hs = np.repeat(np.minimum(q.depth_plot_units.values, YMAX) / YMAX * 0.36, 2)
            ax.fill_between(xs, y + sgn * 0.08, y + sgn * (0.08 + hs), color=col, lw=0, zorder=1)
        for _, g in gaps[gaps.chrom == cname].iterrows():
            w_ = max((g.end_bp - g.start_bp) / 1e6, 1.2)          # drawn at least 1.2 Mb wide to be visible
            mid = (g.start_bp + g.end_bp) / 2e6
            ax.add_patch(Rectangle((mid - w_ / 2, y - 0.34), w_, 0.68, facecolor=GAP_C, edgecolor="none", zorder=3))
    ax.set_xlim(0, xmax); ax.set_ylim(-0.6, 15.75)
    ax.set_yticks(range(16), [f"{i + 1:02d}{side}" for i in range(16)][::-1])
    ax.tick_params(axis="y", length=0, pad=1, labelsize=5.2)
    ax.set_xticks([0, 50, 100, 150], ["0", "50", "100", "150"])
    ax.tick_params(axis="x", labelsize=5.2, pad=1)
    ax.set_xlabel("Position (Mb)", labelpad=1, fontsize=5.8)
    for s_ in ("top", "right", "left"):
        ax.spines[s_].set_visible(False)
    ax.set_title(f"FL-Hap{'1' if side == 'A' else '2'} (chr01{side}–chr16{side})", fontsize=6, pad=2)
letter(0, row_top - 3, "c")
hand = [mpl.patches.Patch(color=HIFI_C, label="HiFi read coverage"),
        mpl.patches.Patch(color=ONT_C, label="ONT read coverage"),
        mpl.patches.Patch(color=GAP_C, label="Filled gap")]
fig.legend(handles=hand, loc="lower left", bbox_to_anchor=(11 / W, 1 - (row_top + c_h + 1.5) / H), ncol=3,
           frameon=False, fontsize=5.5, handlelength=1.0, handleheight=0.7, columnspacing=1.2)
fig_c = gaps.assign(start_Mb=(gaps.start_bp / 1e6).round(2)).values.tolist()

# ---------------- e (redraw2: KR-normalised Pore-C matrix) ----------------
ex = 3 + c_w + 5
e_side = 41
npz = np.load(SRC2 / "ED1e/ED1e_KR_500kb.npz")
kr, vmax, bs = npz["kr"], float(npz["vmax"]), int(npz["bin_size"])
chroms, sizes = [str(c) for c in npz["chroms"]], npz["sizes"].astype(int)
assert abs(vmax - 0.008289517522942717) < 1e-12        # same matrix and normalisation as the published map
nb = np.ceil(sizes / bs).astype(int)                    # HapHiC draw_heatmap uses ceil(size / bin) per scaffold
edges = np.r_[0, np.cumsum(nb)]
assert edges[-1] <= kr.shape[0]
kr = kr[:edges[-1], :edges[-1]]
ax = ax_mm(ex + 8, row_top + 1, e_side, e_side)
cmap = mpl.colors.LinearSegmentedColormap.from_list("wr", ["#FFFFFF", "#F2A38E", "#C8504A", "#8E2A26"])
ax.imshow(np.minimum(kr, vmax), cmap=cmap, vmin=0, vmax=vmax, origin="lower", interpolation="antialiased",
          extent=(0, edges[-1] * bs / 1e6, 0, edges[-1] * bs / 1e6))
for e_ in edges[1:-1]:
    ax.axvline(e_ * bs / 1e6, color="#C9CDD2", lw=0.25, zorder=3)
    ax.axhline(e_ * bs / 1e6, color="#C9CDD2", lw=0.25, zorder=3)
pair_mid = [(edges[2 * i] + edges[2 * i + 2]) / 2 * bs / 1e6 for i in range(16)]
ax.set_yticks(pair_mid, [f"{i + 1:02d}" for i in range(16)]); ax.tick_params(axis="y", labelsize=5.0, pad=1, length=1.5)
ax.set_xticks(np.arange(0, edges[-1] * bs / 1e6, 1000), ["0", "1,000", "2,000", "3,000"])
ax.tick_params(axis="x", labelsize=5.2, pad=1)
ax.set_xlabel("Position (Mb)", labelpad=1, fontsize=5.8)
ax.set_ylabel("Chromosome (A and B homologues)", labelpad=1, fontsize=5.8)
ax.set_title("FL Pore-C contacts (500-kb bins)", fontsize=6, pad=2)
for s_ in ax.spines.values():
    s_.set_linewidth(0.5)
cax = ax_mm(ex + 8 + e_side + 1.5, row_top + 9, 1.6, e_side - 16)
cb = mpl.colorbar.ColorbarBase(cax, cmap=cmap, norm=mpl.colors.Normalize(0, vmax))
cb.set_ticks([0, 0.004, 0.008]); cb.outline.set_linewidth(0.4)
cax.tick_params(labelsize=5.0, length=1.5, width=0.4, pad=1)
cax.set_ylabel("KR-normalized contacts", fontsize=5.2, labelpad=2)
letter(ex, row_top - 3, "e")
extent = edges[-1] * bs / 1e6
pd.DataFrame({"chromosome": chroms, "length_bp": sizes, "bins_500kb": nb, "first_bin": edges[:-1]}).to_csv(
    WORK / "data/ED1e_bin_table.tsv", sep="\t", index=False)

# ---------------- d ----------------
# 2026-09-24: IGV screenshots replaced by read alignments redrawn from the HiFi/ONT remapping BAMs (ed1d_panel.py)
d_top = row_top + c_h + 5
d_log = draw_d(fig, W, H, d_top)
letter(0, d_top - 1, "d")

fig.savefig(WORK / "_ED1_base.pdf", dpi=600, facecolor="white")
# overlay the vector panel a, then render the 600-dpi PNG from the final PDF
out = fitz.open(WORK / "_ED1_base.pdf"); pg = out[0]
k = 72 / 25.4
x, y, w, h = A_RECT_MM
pg.show_pdf_page(fitz.Rect(x * k, y * k, (x + w) * k, (y + h) * k), fitz.open(A_PDF), 0)
out.subset_fonts(); out.save(OUT / "Extended_Data_Fig_01.pdf", garbage=4, deflate=True)
pix = fitz.open(OUT / "Extended_Data_Fig_01.pdf")[0].get_pixmap(dpi=600, alpha=False)
pix.set_dpi(600, 600); pix.save(OUT / "Extended_Data_Fig_01.png")
write_source_data(d_log, WORK / "data")
print(fig_c); print("a h mm", round(ah, 1), "extent Mb", round(extent, 1), "vmax", round(vmax, 5))
