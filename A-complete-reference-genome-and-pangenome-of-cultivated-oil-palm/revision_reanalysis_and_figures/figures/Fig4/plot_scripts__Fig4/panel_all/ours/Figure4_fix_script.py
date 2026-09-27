#!/usr/bin/env python3
"""Figure 4 label corrections (vector text replacement with PyMuPDF).

Every replaced label is removed with a text-only redaction and rewritten on the
original baseline, size and alignment in Arial (italic for species names).
Graphics and images are untouched.
"""
import sys
from pathlib import Path

import fitz

SRC, DST = Path(sys.argv[1]), Path(sys.argv[2])
FONTS = {
    "reg": "/System/Library/Fonts/Supplemental/Arial.ttf",
    "it": "/System/Library/Fonts/Supplemental/Arial Italic.ttf",
    "bold": "/System/Library/Fonts/Supplemental/Arial Bold.ttf",
}
FOBJ = {k: fitz.Font(fontfile=v) for k, v in FONTS.items()}

# (panel region x0,x1,y0,y1), old text, new runs [(text, style)], alignment
EG = [183, 65, 176, 146, 25, 17, 90, 57, 58, 67, 72, 71, 86, 107, 113, 83, 62, 75, 95, 33, 102, 8, 35, 15, 41, 37]
G = (0, 172, 250, 338)
J = (0, 226, 370, 612)
K = (226, 520, 530, 580)
D = (140, 172, 125, 180)
EDITS = [(G, str(n), [(f"EG_{n:03d}", "reg")], "right") for n in EG] + [
    (G, "Dura", [("TK", "reg")], "right"),
    (G, "Pisifera", [("NS", "reg")], "right"),
    (G, "TN (boke)", [("TN", "reg")], "right"),
    (G, "Houke", [("EG_houke", "reg")], "right"),
    (G, "FL (seedless)", [("FL", "reg")], "right"),
    (G, "Meizhou4", [("E. oleifera", "it")], "right"),
    (G, "NRLY", [("Nigerian", "reg")], "right"),
    (J, "EG_dura", [("TK", "reg")], "right"),
    (J, "EG_pisifera", [("NS", "reg")], "right"),
    (J, "bk_hap1", [("TN-Hap1", "reg")], "right"),
    (J, "bk_hap2", [("TN-Hap2", "reg")], "right"),
    (J, "Africa_hap2", [("FL-Hap2", "reg")], "right"),
    (J, "American_hap1", [("FL-Hap1", "reg")], "right"),
    (J, "EG_niriliya", [("Nigerian", "reg")], "right"),
    (J, "African/American haplotypes", [("FL haplotypes", "reg")], "left"),
    (J, "BK haplotypes", [("TN haplotypes", "reg")], "left"),
    (J, "Reference genomes", [("TK, NS, Nigerian and EG_houke", "reg")], "left"),
    (K, "Pisifera (NS)", [("NS", "reg")], "vtop"),
    (K, "Dura (TK)", [("TK", "reg")], "vtop"),
    (K, "E. oleifera", [("E. oleifera", "it")], "vtop"),
    (D, "Dura", [("TK", "reg")], "left"),
    (D, "Pisifera", [("NS", "reg")], "left"),
    (D, "Oleifera", [("E. oleifera", "it")], "left"),
]


def spans(page):
    for b in page.get_text("dict")["blocks"]:
        for l in b.get("lines", []):
            for s in l["spans"]:
                yield s, l["dir"]


def inside(bbox, reg):
    x0, x1, y0, y1 = reg
    return x0 <= bbox[0] < x1 and y0 <= bbox[1] < y1


doc = fitz.open(SRC)
page = doc[0]
all_spans = list(spans(page))
jobs, log = [], []
for reg, old, runs, align in EDITS:
    hits = [(s, d) for s, d in all_spans if s["text"].strip() == old and inside(s["bbox"], reg)]
    if len(hits) != 1:
        raise SystemExit(f"{old!r}: expected 1 match in {reg}, found {len(hits)}")
    s, d = hits[0]
    jobs.append((s, d, runs, align, old))
    # narrow band inside the glyphs of this line only (rows are tighter than their bboxes)
    x0, y0, x1, y1 = s["bbox"]
    if abs(d[0]) < 0.5:          # rotated text: band across x
        ox = s["origin"][0]
        r = fitz.Rect(ox - 0.55 * s["size"], y0 + 0.2, ox - 0.2 * s["size"], y1 - 0.2)
    else:
        oy = s["origin"][1]
        r = fitz.Rect(x0 + 0.2, oy - 0.55 * s["size"], x1 - 0.2, oy - 0.2 * s["size"])
    page.add_redact_annot(r, fill=None)
page.apply_redactions(images=fitz.PDF_REDACT_IMAGE_NONE, graphics=fitz.PDF_REDACT_LINE_ART_NONE,
                      text=fitz.PDF_REDACT_TEXT_REMOVE)

for k, f in FONTS.items():
    if k != "bold":
        page.insert_font(fontname=f"F4{k}", fontfile=f)

for s, d, runs, align, old in jobs:
    size = s["size"]
    c = s["color"]
    col = ((c >> 16 & 255) / 255, (c >> 8 & 255) / 255, (c & 255) / 255)
    ox, oy = s["origin"]
    width = sum(FOBJ[st].text_length(t, fontsize=size) for t, st in runs)
    vertical = abs(d[0]) < 0.5
    if align == "right":
        x = s["bbox"][2] - width
        pts = [(x, oy)]
    elif align == "left":
        pts = [(ox, oy)]
    elif align == "vtop":            # text rotated 90° (reads bottom-to-top), top end kept fixed
        top = s["bbox"][1]
        pts = [(ox, top + width)]
    x, y = pts[0]
    for t, st in runs:
        if vertical:
            page.insert_text((x, y), t, fontsize=size, fontname=f"F4{st}", rotate=90, color=col)
            y -= FOBJ[st].text_length(t, fontsize=size)
        else:
            page.insert_text((x, y), t, fontsize=size, fontname=f"F4{st}", color=col)
            x += FOBJ[st].text_length(t, fontsize=size)
    log.append((old, "".join(t for t, _ in runs), round(size, 2)))

doc.subset_fonts()
doc.save(DST, garbage=3, deflate=True)
for r in log:
    print("\t".join(map(str, r)))
