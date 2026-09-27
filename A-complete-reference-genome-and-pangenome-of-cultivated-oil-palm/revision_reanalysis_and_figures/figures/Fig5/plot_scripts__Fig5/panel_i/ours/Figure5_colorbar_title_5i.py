#!/usr/bin/env python3
"""Vector pass on Figure 5i: add the colour-bar title and shift the two legend lines to make room.

The 5i heat map colour is ln(1 + total load), total load = dSV + frequency-inclusive dSNP count
(author Source Data 4 of Source_Data_dSV_dSNP_African35.xlsx). Called after minus_fix in fix_figure5.py.
"""
import sys
from pathlib import Path

import fitz

HERE = Path(__file__).resolve().parent
ARIAL = str(HERE / "ArialFix.ttf")
TITLE = "ln(1 + dSV + dSNP)"
BAR_RIGHT = 372.5          # visible right end of the 5i colour bar (pt)
BAR_Y = (677.8, 684.0)     # colour-bar strip
GAP = 3.0
LEGEND_TEXTS = ("Labels: dSV/dSNP", "Red boxes: selected segments")


def add_title(path):
    path = Path(path)
    doc = fitz.open(path)
    pg = doc[0]
    if TITLE in pg.get_text():
        print("colorbar_title: already present")
        return
    spans = []
    for b in pg.get_text("dict")["blocks"]:
        for l in b.get("lines", []):
            for s in l["spans"]:
                if s["text"].strip() in LEGEND_TEXTS and s["bbox"][1] > 660:
                    spans.append(s)
    assert len(spans) == 2, spans
    font = fitz.Font(fontfile=ARIAL)
    size = spans[0]["size"]
    tw_title = font.text_length(TITLE, fontsize=size)
    x_title = BAR_RIGHT + GAP
    new_legend_x = x_title + tw_title + 7.0
    dx = new_legend_x - min(s["origin"][0] for s in spans)
    right = max(s["bbox"][2] for s in spans) + dx
    assert right < pg.rect.width - 4, right
    for s in spans:
        x0, y0, x1, y1 = s["bbox"]
        h = y1 - y0
        pg.add_redact_annot(fitz.Rect(x0 + 0.1, y0 + 0.25 * h, x1 - 0.1, y1 - 0.25 * h))
    pg.apply_redactions(images=fitz.PDF_REDACT_IMAGE_NONE, graphics=fitz.PDF_REDACT_LINE_ART_NONE,
                        text=fitz.PDF_REDACT_TEXT_REMOVE)
    c = spans[0]["color"]
    rgb = ((c >> 16) & 255) / 255, ((c >> 8) & 255) / 255, (c & 255) / 255
    tw = fitz.TextWriter(pg.rect, color=rgb)
    for s in spans:
        tw.append((s["origin"][0] + dx, s["origin"][1]), s["text"], font=font, fontsize=s["size"])
    base = (BAR_Y[0] + BAR_Y[1]) / 2 + 0.36 * size
    tw.append((x_title, base), TITLE, font=font, fontsize=size)
    tw.write_text(pg)
    tmp = path.with_suffix(".cbtitle.pdf")
    doc.save(tmp, garbage=3, deflate=True)
    doc.close()
    tmp.replace(path)
    print(f"colorbar_title: title at x={x_title:.1f}, legend shifted by {dx:.1f} pt (right edge {right:.1f})")


if __name__ == "__main__":
    add_title(sys.argv[1])
