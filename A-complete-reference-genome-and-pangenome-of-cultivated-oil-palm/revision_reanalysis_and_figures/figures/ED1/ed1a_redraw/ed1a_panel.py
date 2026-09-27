#!/usr/bin/env python3
"""Extended Data Fig. 1a (study design) redrawn in matplotlib as vector art, 2026-09-25.

draw_a(fig, W, H, x0, y0, w, h) draws the panel into the 92 x 80.44 mm cell of ED1 (x0, y0 = top-left, mm).
Text: Arial, body 5-5.5 pt, layer titles 6.5 pt bold.  Colours: scheme A (fix/beautify/common/palA.py).
Photographs and icons are the original bitmaps of the author's flowchart (ED1a_study_design_arial.pdf),
extracted at native resolution (work/extract_imgs.py -> work/img/x<xref>.png) and embedded unresampled
(interpolation='none').  Fruit photographs stored as JPEG on white get an alpha channel from a flood fill of the
near-white border region (RGB values unchanged).  The vector overlays of the icons (chromosome bands, monitor trace)
are replayed from the source PDF.  Text is taken verbatim from the current Extended_Data_Fig_01.pdf (after the
polish_fix hybrid-order edit); long explanatory sentences and the footnotes move to the legend (legend_replacement_ED1a.tsv).
"""
from pathlib import Path
import sys

import fitz
import matplotlib as mpl
import numpy as np
from matplotlib import font_manager as fm
from matplotlib.ft2font import LOAD_NO_HINTING
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, PathPatch
from matplotlib.path import Path as MPath
from PIL import Image
from scipy import ndimage

HERE = Path(__file__).resolve().parent
FIX = HERE.parent
sys.path.insert(0, str(FIX / "beautify/common"))
import palA  # noqa: E402

SRC_PDF = FIX / "beautify/work/ED1/ED1a_study_design_arial.pdf"
IMG = HERE / "work/img"
PT = 25.4 / 72                      # mm per point
INK, INK2, FRAME = "#222222", "#5A5A5A", "#8E99A6"
LW_FRAME, LW_BOX, LW_ARROW = 0.5, 0.5, 0.6

MAT = {"TK": palA.TK, "NS": palA.NS, "TN": palA.TN, "NIG": palA.NIG, "EOL": palA.EOL, "FL": palA.FL,
       "EG": palA.EG11}


def dark(hexc, f=0.78):
    h = hexc.lstrip("#")
    return "#%02X%02X%02X" % tuple(round(int(h[i:i + 2], 16) * f) for i in (0, 2, 4))


# ---------------------------------------------------------------- rich text (runs of regular / italic)
_FONTS = {}


def _font(bold, italic):
    k = (bold, italic)
    if k not in _FONTS:
        prop = fm.FontProperties(family="Arial", weight="bold" if bold else "normal",
                                 style="italic" if italic else "normal")
        _FONTS[k] = (prop, mpl.ft2font.FT2Font(fm.findfont(prop)))
    return _FONTS[k]


def adv(s, size, bold=False, italic=False):
    """advance width of s in mm (sum of glyph advances, as placed by the PDF backend)"""
    _, f = _font(bold, italic)
    f.set_size(size, 72)
    return sum(f.load_char(ord(c), flags=LOAD_NO_HINTING).linearHoriAdvance / 65536 for c in s) * PT


def parse(s):
    """'Deli *dura*' -> [('Deli ', False), ('dura', True)]"""
    out, it = [], False
    for i, part in enumerate(s.split("*")):
        if part:
            out.append((part, it))
        it = not it
    return out


LOG = []                            # (text, size, bold, italic) of every run drawn


class Panel:
    def __init__(self, fig, W, H, x0, y0, w, h):
        self.ax = fig.add_axes([x0 / W, 1 - (y0 + h) / H, w / W, h / H])
        self.ax.set_xlim(0, w); self.ax.set_ylim(h, 0); self.ax.axis("off")
        self.w, self.h = w, h

    def width(self, s, size, bold=False):
        return sum(adv(t, size, bold, it) for t, it in parse(s))

    def text(self, x, y, s, size, bold=False, ha="center", color=INK, z=5, bg=None):
        """y = baseline (mm); s uses *...* for italic"""
        wtot = self.width(s, size, bold)
        xs = x - {"left": 0, "center": wtot / 2, "right": wtot}[ha]
        for t, it in parse(s):
            prop, _ = _font(bold, it)
            tt = self.ax.text(xs, y, t, fontproperties=prop, fontsize=size, color=color, ha="left",
                              va="baseline", zorder=z)
            if bg:
                tt.set_bbox(dict(facecolor=bg, edgecolor="none", pad=0.6))
            LOG.append((t, size, bold, it))
            xs += adv(t, size, bold, it)
        return wtot

    def lines(self, x, y, rows, size, lead=1.18, **kw):
        """rows: list of strings or (string, size, bold, color); returns the baseline after the last row"""
        for r in rows:
            s, sz, b, c = (r, size, False, kw.get("color", INK)) if isinstance(r, str) else r
            self.text(x, y, s, sz, bold=b, color=c, ha=kw.get("ha", "center"))
            y += sz * lead * PT
        return y

    def box(self, x, y, w, h, fc, ec, lw=LW_BOX, r=0.7, ls="-", z=1):
        self.ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle=f"round,pad=0,rounding_size={r}", fc=fc, ec=ec,
                                         lw=lw, ls=ls, zorder=z))

    def arrow(self, pts, dashed=False, head=True, lw=LW_ARROW, color="#333333"):
        xy = np.asarray(pts, float)
        ls = (0, (2.2, 1.6)) if dashed else "-"
        if len(xy) > 2:
            self.ax.plot(xy[:-1, 0], xy[:-1, 1], color=color, lw=lw, ls=ls, zorder=4, solid_capstyle="butt",
                         dash_capstyle="butt")
            xy = xy[-2:]
        self.ax.add_patch(FancyArrowPatch(xy[0], xy[1], arrowstyle="-|>,head_length=2.2,head_width=1.2" if head else "-",
                                          mutation_scale=1, lw=lw, ls=ls, color=color, shrinkA=0, shrinkB=0,
                                          zorder=4))

    def cross(self, x, y, size=7):
        self.text(x, y + size * 0.36 * PT, "×", size, bold=True, color="#333333")

    # ------------------------------------------------------------ bitmaps / vector icons from the source PDF
    def group(self, src, x, y, h=None, w=None, ha="center", va="top", vec=False):
        """draw every source image (and vector icon path) inside PDF rect src, scaled into a box at (x, y)"""
        src = fitz.Rect(src)
        s = h / src.height if h else w / src.width
        dw, dh = src.width * s, src.height * s
        x0 = x - {"left": 0, "center": dw / 2, "right": dw}[ha]
        y0 = y - {"top": 0, "center": dh / 2, "bottom": dh}[va]
        m = lambda px, py: (x0 + (px - src.x0) * s, y0 + (py - src.y0) * s)  # noqa: E731
        for xref, bb in _SRC_IMAGES:
            if src.contains(bb):
                (l, t), (r, b) = m(bb.x0, bb.y0), m(bb.x1, bb.y1)
                self.ax.imshow(_img(xref), extent=(l, r, b, t), interpolation="none", zorder=3, aspect="auto")
        for dr in _SRC_DRAW:
            if vec and src.contains(dr["rect"]) and dr["rect"].width < src.width * 0.95:
                verts, codes = [], []
                for it in dr["items"]:
                    if it[0] == "l":
                        pts, cd = [it[1], it[2]], [MPath.MOVETO, MPath.LINETO]
                    elif it[0] == "c":
                        pts, cd = list(it[1:5]), [MPath.MOVETO] + [MPath.CURVE4] * 3
                    elif it[0] == "re":
                        q = it[1]
                        pts = [q.tl, q.tr, q.br, q.bl, q.tl]; cd = [MPath.MOVETO] + [MPath.LINETO] * 4
                    elif it[0] == "qu":
                        q = it[1]
                        pts = [q.ul, q.ur, q.lr, q.ll, q.ul]; cd = [MPath.MOVETO] + [MPath.LINETO] * 4
                    else:
                        continue
                    if verts and codes and np.allclose(verts[-1], m(*pts[0])):
                        pts, cd = pts[1:], cd[1:]
                    verts += [m(p.x, p.y) for p in pts]; codes += cd
                fc = dr.get("fill"); ec = dr.get("color")
                self.ax.add_patch(PathPatch(MPath(verts, codes), fc=fc if fc else "none",
                                            ec=ec if ec else "none", lw=(dr.get("width") or 0) * s / PT,
                                            zorder=3.5))
        return x0, y0, dw, dh


_doc = fitz.open(SRC_PDF)
_pg = _doc[0]
_SRC_IMAGES = [(i["xref"], fitz.Rect(i["bbox"])) for i in _pg.get_image_info(xrefs=True)]
_SRC_DRAW = _pg.get_drawings()
_WHITE_BG = {12, 15, 20, 21, 22, 23, 24, 25}   # fruit photographs on a white ground (JPEG, or masked with the white square kept)


def _img(xref):
    im = np.asarray(Image.open(IMG / f"x{xref}.png").convert("RGBA")).copy()
    if xref in _WHITE_BG:
        near = im[..., :3].min(2) > 232
        lab, _ = ndimage.label(near)
        border = np.unique(np.r_[lab[0], lab[-1], lab[:, 0], lab[:, -1]])
        bg = np.isin(lab, border[border > 0])
        a = np.where(bg, 0, 255).astype(float)
        a = ndimage.gaussian_filter(a, 0.6)              # soften the cut-out edge by ~1 px
        im[..., 3] = np.minimum(im[..., 3], np.clip(a, 0, 255).astype(np.uint8))
    return im


# source rectangles (PDF points of ED1a_study_design_arial.pdf)
R = dict(eg_tree=(389, 21, 497, 102), eo_tree=(891, 22, 993, 99),
         dura=(201.5, 202.5, 343, 312), pisi=(437.5, 205.5, 567, 343), ten=(644, 187.5, 781, 324.5),
         TK=(292.5, 599.5, 345, 652), NS=(421, 592.5, 487, 655.5), TN=(562, 601.5, 622, 658),
         NIG=(699.5, 590, 762, 656), EOL=(840.5, 595.5, 890, 651), FL=(969.5, 586.5, 1035.5, 670.5),
         chrom=(248.5, 711, 278, 798), bar=(440.5, 708, 452.5, 805), seq=(642.5, 717.5, 702.5, 759.5),
         mon=(846.5, 723.5, 917.5, 781.5))


def draw_a(fig, W, H, x0, y0, w, h):
    LOG.clear()
    P = Panel(fig, W, H, x0, y0, w, h)
    TSZ, BSZ, NSZ = 6.2, 5.0, 5.5              # layer title, body, box names
    LX0, LX1 = 0.25, w - 0.25                  # layer frames
    IN0, IN1 = 1.4, w - 1.4                    # content inside a frame
    GAP = 2.2
    Hs = [9.6, 15.4, 13.2, 17.8, 14.0]
    tops = []
    y = 1.35
    for hh in Hs:
        tops.append(y); y += hh + GAP
    assert tops[-1] + Hs[-1] <= h + 1e-6, tops[-1] + Hs[-1]
    titles = ["1. Species", "2. SHELL-associated fruit forms", "3. Breeding backgrounds/diversity sources",
              "4. Six haplotype-resolved accessions used in this study", "5. Genomic resources and datasets"]
    subt = [None, "(defined by the SHELL locus)", None, "(12 haplotypes in total)", None]
    for t, hh, s, sub in zip(tops, Hs, titles, subt):
        P.box(LX0, t, LX1 - LX0, hh, "#FBFBFC", FRAME, lw=LW_FRAME, r=1.0, z=0)
        tw = P.width(s, TSZ, True)
        bw = tw + (P.width(" " + sub, BSZ) if sub else 0)
        P.box(1.2, t - 1.35, bw + 1.4, 2.7, "white", "none", lw=0, r=0.4, z=0.5)      # tab behind the title
        P.text(1.8, t + 0.78, s, TSZ, bold=True, ha="left")
        if sub:
            P.text(1.8 + tw, t + 0.78, " " + sub, BSZ, ha="left", color=INK2)

    # ---- 1. Species
    t, hh = tops[0], Hs[0]
    by, bh = t + 2.1, hh - 3.0
    half = (IN1 - IN0 - 1.2) / 2
    sp = [("African oil palm", "*Elaeis guineensis*", "eg_tree", MAT["EG"]),
          ("American oil palm", "*Elaeis oleifera*", "eo_tree", MAT["EOL"])]
    sp_box = []
    for k, (n1, n2, img, col) in enumerate(sp):
        bx = IN0 + k * (half + 1.2)
        P.box(bx, by, half, bh, palA.tint(col, 0.10), col)
        _, _, pw, _ = P.group(R[img], bx + half - 1.0, by + 0.8, h=bh - 1.6, ha="right")
        cx = bx + (half - pw - 1.0) / 2
        P.text(cx, by + bh / 2 - 0.35, n1, 6.0, bold=True)
        P.text(cx, by + bh / 2 + 2.25, n2, 6.0, bold=True, color=dark(col, 0.9) if k else INK)
        sp_box.append((bx, by, half, bh))

    # ---- 2. SHELL-associated fruit forms
    t, hh = tops[1], Hs[1]
    sc = 8.0 / 137.0                             # one scale for the three fruit bunches (mm per source pt)
    lab_y = t + 3.6
    ph_top = t + 6.2
    forms = [("*Dura* (thick-shelled)", "*Sh+/Sh+*", "dura", 11.0, MAT["TK"]),
             ("*Pisifera* (shell-less)", "*sh−/sh−*", "pisi", 33.5, MAT["NS"]),
             ("*Tenera* (thin-shelled)", "*Sh+/sh−*", "ten", 56.0, MAT["TN"])]
    mid = {}
    for n, g, img, cx, col in forms:
        src = fitz.Rect(R[img])
        P.text(cx, lab_y, n, NSZ, bold=True)
        P.text(cx, lab_y + 2.1, g, NSZ, bold=True, color=dark(col))
        x_, y_, dw, dh = P.group(R[img], cx, ph_top, h=src.height * sc)
        mid[img] = (x_, y_, dw, dh)
    ym = ph_top + 3.8
    P.cross(22.3, ym, 8)
    P.arrow([(mid["pisi"][0] + mid["pisi"][2] + 0.6, ym), (mid["ten"][0] - 0.6, ym)], lw=0.9)
    # interspecific hybrid (tenera x E. oleifera)
    hx = IN0 + 5.5 * (IN1 - IN0 + 1.0) / 6 - 0.5   # centre of the FL box in row 4
    ten_r = mid["ten"][0] + mid["ten"][2] + 0.8
    P.ax.plot([ten_r, hx], [ym, ym], color="#333333", lw=LW_ARROW, ls=(0, (2.2, 1.6)), zorder=4)
    eo_b = sp_box[1][1] + sp_box[1][3]
    P.arrow([(hx, eo_b), (hx, t + hh - 5.0)], dashed=True)
    P.cross(hx - 3.2, ym - 1.3, 6.5)
    P.text(hx - 4.8, t + hh - 2.95, "Interspecific hybrid", NSZ, bold=True)
    P.text(hx - 4.8, t + hh - 0.95, "*E. guineensis* × *E. oleifera*", NSZ, bold=True, color=dark(MAT["FL"], 0.9))

    # ---- 3. Breeding backgrounds / diversity sources
    t, hh = tops[2], Hs[2]
    bg = [("Deli", "(breeding background)", ["Southeast Asian;", "*dura* maternal background"],
           MAT["TK"]),
          ("AVROS", "(breeding background)", ["African germplasm;", "*pisifera* paternal background"], MAT["NS"]),
          ("Nigerian", "(geographical origin/diversity resource)",
           ["West African *E. guineensis* diversity;", "*dura*, *pisifera* and *tenera*"], MAT["NIG"])]
    widths = [23.4, 23.6, 31.6]
    bx = IN0
    by, bh = t + 2.1, hh - 3.0
    for (n, sub, body, col), bw in zip(bg, widths):
        P.box(bx, by, bw, bh, palA.tint(col, 0.10), col)
        cx = bx + bw / 2
        P.lines(cx, by + 2.35, [(n, 6.0, True, dark(col)), (sub, BSZ, False, INK2)] + body, BSZ, lead=1.2)
        bx += bw + 1.2
    SEL_Y = t + hh / 2

    # ---- 4. Six accessions
    t, hh = tops[3], Hs[3]
    acc = [("TK", ["Deli *dura*", "*Sh+/Sh+*"], MAT["TK"]),
           ("NS", ["AVROS *pisifera*", "*sh−/sh−*"], MAT["NS"]),
           ("TN", ["Deli *dura* ×", "AVROS *pisifera*", "*Sh+/sh−*"], MAT["TN"]),
           ("Nigerian", ["West African", "*E. guineensis*", "(*tenera*) *Sh+/sh−*"], MAT["NIG"]),
           ("*E. oleifera*", ["American", "oil palm"], MAT["EOL"]),
           ("FL (Reyou-2)", ["*E. guineensis* ×", "*E. oleifera*"], MAT["FL"])]
    keys = ["TK", "NS", "TN", "NIG", "EOL", "FL"]
    n6 = len(acc)
    g6 = 1.0
    bw = (IN1 - IN0 - g6 * (n6 - 1)) / n6
    by, bh = t + 2.1, hh - 3.0
    fl_cx = None
    for i, ((n, body, col), k) in enumerate(zip(acc, keys)):
        bx = IN0 + i * (bw + g6)
        cx = bx + bw / 2
        P.box(bx, by, bw, bh, palA.tint(col, 0.10), col)
        P.lines(cx, by + 2.4, [(n, 6.0, True, dark(col))] + body, BSZ, lead=1.2)
        P.group(R[k], cx, by + bh - 0.6, h=5.0, va="bottom")
        if k == "FL":
            fl_cx = cx
    # selection arrow: hybrid label -> FL box
    P.arrow([(fl_cx, tops[1] + Hs[1] + 0.2), (fl_cx, by - 0.05)], dashed=True)
    sw = P.width("selection", NSZ, True)
    prop, _ = _font(True, False)
    P.ax.text(fl_cx + 0.9, SEL_Y, "selection", fontproperties=prop, fontsize=NSZ, rotation=90, ha="left",
              va="center", color=INK, zorder=5)
    LOG.append(("selection", NSZ, True, False))

    # ---- 5. Genomic resources and datasets
    t, hh = tops[4], Hs[4]
    by, bh = t + 2.1, hh - 3.0
    res = [("chrom", 1.9, ["12 haplotype-resolved", "assemblies", "from the 6 accessions", "(2 haplotypes per",
                           "accession)"]),
           ("bar", 1.0, ["27 additional", "chromosome-scale", "assemblies", "(33 biological", "materials in total)"]),
           ("seq", 4.8, ["308 deeply", "resequenced", "accessions"]),
           ("mon", 5.2, ["Gene pangenome", "(64,577 gene families)", "+ graph pangenome", "(structural variation", "landscape)"])]
    conn = ["+", "+", "→"]
    cw = 2.3
    need = []
    for key, iw, rows in res:
        tw = max(P.width(r, BSZ) for r in rows)
        need.append(iw + tw + 2.2)
    extra = (IN1 - IN0 - cw * 3 - sum(need)) / 4
    assert extra > -0.01, extra
    bx = IN0
    for j, ((key, iw, rows), nw) in enumerate(zip(res, need)):
        bw_ = nw + extra
        P.box(bx, by, bw_, bh, "#F3F5F7", FRAME)
        src = fitz.Rect(R[key])
        ih = min(bh - 1.6, iw * src.height / src.width)
        P.group(R[key], bx + 0.8 + iw / 2, by + bh / 2, h=ih, va="center", vec=True)
        txt_cx = bx + 0.8 + iw + (bw_ - 0.8 - iw) / 2
        lead = 1.12
        n_ = len(rows)
        y0_ = by + bh / 2 - (n_ - 1) * BSZ * lead * PT / 2 + BSZ * 0.35 * PT
        P.lines(txt_cx, y0_, rows, BSZ, lead=lead)
        bx += bw_
        if j < 3:
            if conn[j] == "→":
                P.arrow([(bx + 0.35, by + bh / 2), (bx + cw - 0.35, by + bh / 2)], lw=0.7)
            else:
                P.text(bx + cw / 2, by + bh / 2 + 1.1, conn[j], 7, bold=True, color="#333333")
            bx += cw
    return P, list(LOG)


if __name__ == "__main__":                      # stand-alone preview of the cell
    import matplotlib.pyplot as plt
    mpl.rcParams.update({"font.family": "Arial", "pdf.fonttype": 42})
    Wc, Hc = 92, 80.44345
    fig = plt.figure(figsize=(Wc / 25.4, Hc / 25.4))
    _, log = draw_a(fig, Wc, Hc, 0, 0, Wc, Hc)
    out = HERE / "work/ED1a_panel_preview"
    fig.savefig(str(out) + ".pdf", dpi=600, facecolor="white")
    fig.savefig(str(out) + ".png", dpi=600, facecolor="white")
    print("min text size", min(l[1] for l in log))
