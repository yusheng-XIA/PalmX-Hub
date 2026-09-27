#!/usr/bin/env python3
"""Final vector pass on Figure 5: ASCII hyphen-minus in negative numbers -> U+2212; FL/TN_hap -> -Hap.

Called at the end of fix_figure5.py so that re-running the build never reintroduces '-1000' (5f).
"""
import re
import sys
from pathlib import Path

import fitz

ARIAL = "/System/Library/Fonts/Supplemental/Arial.ttf"


def newtext(t):
    t = re.sub(r"\b(FL|TN)_hap([12])\b", r"\1-Hap\2", t)
    return re.sub(r"(^|[\s(=<>,])-(\d)", "\\1\u2212\\2", t)


def fix(path):
    path = Path(path)
    doc = fitz.open(path)
    pg = doc[0]
    edits = []
    for b in pg.get_text("rawdict")["blocks"]:
        for l in b.get("lines", []):
            for s in l["spans"]:
                t = "".join(c["c"] for c in s["chars"])
                n = newtext(t)
                if n != t:
                    edits.append((s, t, n))
    if not edits:
        print("minus_fix: nothing to change")
        return []
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
    tmp = path.with_suffix(".minusfix.pdf")
    doc.save(tmp, garbage=3, deflate=True)
    doc.close()
    tmp.replace(path)
    log = [(t, n, round(s["origin"][0], 1), round(s["origin"][1], 1)) for s, t, n in edits]
    print("minus_fix:", log)
    return log


if __name__ == "__main__":
    fix(sys.argv[1])
