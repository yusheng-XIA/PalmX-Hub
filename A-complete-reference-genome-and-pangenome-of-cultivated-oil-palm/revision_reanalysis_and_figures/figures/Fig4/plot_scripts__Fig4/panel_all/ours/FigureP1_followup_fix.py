#!/usr/bin/env python3
"""Follow-up: residual _hap names and ASCII minus signs in revised main figures (vector edit)."""
import re, shutil, sys, collections
from pathlib import Path
import fitz

D = Path("deliver/Main_Figures_revised")
ARIAL = "/System/Library/Fonts/Supplemental/Arial.ttf"


def newtext(t):
    t = re.sub(r"\b(FL|TN)_hap([12])\b", r"\1-Hap\2", t)
    t = re.sub(r"^hap([12])$", r"Hap\1", t)
    t = re.sub(r"(^|[\s(=<>,])-(\d)", "\\1\u2212\\2", t)
    return t


log = []
for f in ["Figure3", "Figure4", "Figure5"]:
    src = D / f"{f}.pdf"
    bak = Path("p1") / f"{f}_before_P1.pdf"
    if not bak.exists():
        shutil.copy2(src, bak)
    doc = fitz.open(bak)
    pg = doc[0]
    before = pg.get_text("words")
    edits = []
    for b in pg.get_text("rawdict")["blocks"]:
        for l in b.get("lines", []):
            for s in l["spans"]:
                t = "".join(c["c"] for c in s["chars"])
                n = newtext(t)
                if n != t:
                    edits.append((s, t, n))
    for s, t, n in edits:
        x0, y0, x1, y1 = s["bbox"]
        h = y1 - y0
        pg.add_redact_annot(fitz.Rect(x0 + 0.1, y0 + 0.25 * h, x1 - 0.1, y1 - 0.25 * h))
    pg.apply_redactions(images=fitz.PDF_REDACT_IMAGE_NONE, graphics=fitz.PDF_REDACT_LINE_ART_NONE,
                        text=fitz.PDF_REDACT_TEXT_REMOVE)
    font = fitz.Font(fontfile=ARIAL)
    for s, t, n in edits:
        c = s["color"]
        rgb = ((c >> 16) & 255) / 255, ((c >> 8) & 255) / 255, (c & 255) / 255
        tw = fitz.TextWriter(pg.rect, color=rgb)
        tw.append(s["origin"], n, font=font, fontsize=s["size"])
        tw.write_text(pg)
        log.append((f, round(s["origin"][0], 1), round(s["origin"][1], 1), t, n))
    out = D / f"{f}.pdf"
    doc.save(out, garbage=3, deflate=True)
    # word-level comparison
    after = fitz.open(out)[0].get_text("words")
    bw = collections.Counter(newtext(w[4]) for w in before)
    aw = collections.Counter(w[4].replace("\u00ad", "-") for w in after)
    miss = bw - aw
    extra = aw - bw
    print(f, "edits", len(edits), "missing", dict(miss), "extra", dict(extra))
for r in log:
    print(r)
