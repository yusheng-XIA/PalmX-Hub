#!/usr/bin/env python3
"""Compose the Figure 2 re-layout candidate (vector throughout).

Source: deliver/Main_Figures_revised/Figure2.pdf (current submission version; not modified).
- old 2a, 2b, 2e, 2f copied as vector; old 2c (metabolite heatmap) and old 2d (axes) region removed;
- new c = old 2d redrawn taller from Source Data (panels/c_axes.pdf) in the freed right column;
- old e,f relettered d,e; new OLE16a row f-h (+ optional placeholder i) appended below;
- in-figure text: 'Postharvest (h)' -> 'Post-harvest (h)' (2e header), 'FL oil storage' -> 'FL lipid synthesis',
  'OLE16' row label in the fold-change table gets a second line '(family)'.
usage: python3 compose_fig2.py   -> Figure2_v1.pdf/.png, Figure2_v1_placeholder.pdf/.png
"""
from pathlib import Path

import fitz

HERE = Path(__file__).resolve().parents[1]
OUT = Path(__file__).resolve().parent
S = HERE.parents[1]
SRC = S / "deliver/Main_Figures_revised/Figure2.pdf"
AR = "/System/Library/Fonts/Supplemental/Arial.ttf"
ARB = "/System/Library/Fonts/Supplemental/Arial Bold.ttf"
fR, fB = fitz.Font(fontfile=AR), fitz.Font(fontfile=ARB)
W = 510.2359924316406
H0 = 588.8569946289062
ROW_Y, ROW_H = 590.0, 108.0
REGION_CD = fitz.Rect(262.6, 190.8, W, 383.2)
LAB_X, LAB_Y = 364.0, 516.5


def spans(pg):
    return [s for b in pg.get_text("dict")["blocks"] for l in b.get("lines", []) for s in l["spans"]]


def find(pg, text, x=None, y=None):
    out = [s for s in spans(pg) if s["text"].strip() == text
           and (x is None or abs(s["bbox"][0] - x) < 2) and (y is None or abs(s["bbox"][1] - y) < 2)]
    assert len(out) == 1, (text, x, y, len(out))
    return out[0]


def strike(pg, s):
    ox, oy = s["origin"]
    yc = oy - 0.33 * s["size"]
    pg.add_redact_annot(fitz.Rect(s["bbox"][0] + 0.3, yc - 0.12, s["bbox"][2] - 0.3, yc + 0.12), fill=False)


def edited_source():
    doc = fitz.open(SRC)
    pg = doc[0]
    # 1) remove old c/d region entirely (text, vector paths fully inside, image pixels)
    pg.add_redact_annot(REGION_CD, fill=False)
    pg.apply_redactions(images=fitz.PDF_REDACT_IMAGE_PIXELS, graphics=fitz.PDF_REDACT_LINE_ART_REMOVE_IF_COVERED,
                        text=fitz.PDF_REDACT_TEXT_REMOVE)
    # 2) text edits outside the region
    writes = []
    s = find(pg, "Postharvest (h)", 194.3, 387.3)
    c = (s["bbox"][0] + s["bbox"][2]) / 2
    strike(pg, s); writes.append(("Post-harvest (h)", c, s["origin"][1], s["size"], False, s["color"], "center"))
    s = find(pg, "FL oil storage")
    strike(pg, s)   # v2: label moved next to the C9 cloud it describes (lower-left edge of C9)
    writes.append(("FL lipid synthesis", LAB_X, LAB_Y, s["size"], False, s["color"], "center"))
    s = find(pg, "OLE16", 37.7, 541.5)
    strike(pg, s)
    writes.append(("OLE16", s["bbox"][2], s["origin"][1] - 2.4, s["size"], True, s["color"], "right"))
    writes.append(("(family)", s["bbox"][2], s["origin"][1] + 4.2, 5.0, False, s["color"], "right"))
    for old, new, x, y in (("e", "d", 1.8, 383.8), ("f", "e", 263.6, 383.8)):
        s = find(pg, old, x, y)
        strike(pg, s); writes.append((new, s["bbox"][0], s["origin"][1], s["size"], True, s["color"], "left"))
    pg.apply_redactions(images=fitz.PDF_REDACT_IMAGE_NONE, graphics=fitz.PDF_REDACT_LINE_ART_NONE,
                        text=fitz.PDF_REDACT_TEXT_REMOVE)
    for txt, x, y, sz, bold, col, anch in writes:
        f = fB if bold else fR
        w = f.text_length(txt, fontsize=sz)
        x0 = {"left": x, "right": x - w, "center": x - w / 2}[anch]
        rgb = ((col >> 16) & 255) / 255, ((col >> 8) & 255) / 255, (col & 255) / 255
        pg.insert_text((x0, y), txt, fontname="zArB" if bold else "zArR", fontfile=ARB if bold else AR,
                       fontsize=sz, color=rgb)
    return doc


def compose(row_pdf, out_stem):
    src = edited_source()
    H = ROW_Y + ROW_H
    out = fitz.open()
    pg = out.new_page(width=W, height=H)
    pg.show_pdf_page(fitz.Rect(0, 0, W, H0), src, 0)
    cpan = fitz.open(HERE / "panels/c_axes.pdf")
    r = cpan[0].rect
    pg.show_pdf_page(fitz.Rect(263.0, 191.5, 263.0 + r.width, 191.5 + r.height), cpan, 0)  # contains letter c
    rw = fitz.open(OUT / row_pdf)
    rr = rw[0].rect
    pg.show_pdf_page(fitz.Rect(0, ROW_Y, rr.width, ROW_Y + rr.height), rw, 0)
    out.save(OUT / f"{out_stem}.pdf", garbage=4, deflate=True)
    pg.get_pixmap(dpi=300).save(OUT / f"{out_stem}.png")
    pg.get_pixmap(dpi=150).save(OUT / f"{out_stem}_preview.png")
    print(out_stem, f"{W / 72 * 25.4:.1f} x {H / 72 * 25.4:.1f} mm")


if __name__ == "__main__":
    compose("row_ole16_v2.pdf", "Figure2_v2")
