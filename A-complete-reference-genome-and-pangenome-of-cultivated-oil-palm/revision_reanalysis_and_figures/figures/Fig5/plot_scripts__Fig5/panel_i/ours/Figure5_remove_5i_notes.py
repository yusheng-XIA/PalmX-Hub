#!/usr/bin/env python3
"""Vector pass on Figure 5i: remove the two in-figure notes "Red boxes: selected segments" and
"Labels: dSV/dSNP" (their content is carried by the figure legend). Text-only redaction; graphics untouched.
The colour-bar title ln(1 + dSV + dSNP) is kept. Called after colorbar_title_5i in fix_figure5.py.
"""
import sys
from pathlib import Path

import fitz

NOTES = ("Red boxes: selected segments", "Labels: dSV/dSNP")


def remove(path):
    path = Path(path)
    doc = fitz.open(path)
    pg = doc[0]
    spans = []
    for b in pg.get_text("dict")["blocks"]:
        for l in b.get("lines", []):
            for s in l["spans"]:
                if s["text"].strip() in NOTES and s["bbox"][1] > 660:
                    spans.append(s)
    if not spans:
        print("remove_5i_notes: already removed")
        doc.close()
        return
    assert len(spans) == 2, spans
    for s in spans:
        x0, y0, x1, y1 = s["bbox"]
        h = y1 - y0
        pg.add_redact_annot(fitz.Rect(x0 + 0.1, y0 + 0.25 * h, x1 - 0.1, y1 - 0.25 * h))
    pg.apply_redactions(images=fitz.PDF_REDACT_IMAGE_NONE, graphics=fitz.PDF_REDACT_LINE_ART_NONE,
                        text=fitz.PDF_REDACT_TEXT_REMOVE)
    tmp = path.with_suffix(".no5inotes.pdf")
    doc.save(tmp, garbage=3, deflate=True)
    doc.close()
    tmp.replace(path)
    print("remove_5i_notes: removed", [s["text"] for s in spans])


if __name__ == "__main__":
    remove(sys.argv[1])
