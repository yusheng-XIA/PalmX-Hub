"""Fig. 1b: restore EG11 / EO12 telomere-positive chromosome ends from the author's tidk re-run.
Input (read-only): deliver/Main_Figures_revised/Figure1.pdf
The two '–' placeholders in the telomere column are redacted (text only) and rewritten at the
same origin, size (5.57 pt) and colour (#222222) in Arial, left-aligned like the other cells."""
import sys, fitz
from pathlib import Path
S = Path(__file__).resolve().parents[2]
SRC = S / "deliver/Main_Figures_revised/Figure1.pdf"
vals = {"EG11": sys.argv[1], "EO12": sys.argv[2]}
out_path = Path(__file__).parent / sys.argv[3]
ARIAL = "/System/Library/Fonts/Supplemental/Arial.ttf"
FONT = fitz.Font(fontfile=ARIAL)
TARGET = {"EG11": (461.72, 301.5), "EO12": (461.72, 396.36)}

doc = fitz.open(SRC); page = doc[0]
spans = [s for b in page.get_text("dict")["blocks"] for l in b.get("lines", []) for s in l["spans"]]
hits = {}
for k, (x, y) in TARGET.items():
    m = [s for s in spans if s["text"].strip() == "–" and abs(s["origin"][0]-x) < .1 and abs(s["origin"][1]-y) < .1]
    assert len(m) == 1, (k, m); hits[k] = m[0]
for s in hits.values():
    x0, y0, x1, y1 = s["bbox"]
    page.add_redact_annot(fitz.Rect(x0, y0, x1, y1))
page.apply_redactions(images=fitz.PDF_REDACT_IMAGE_NONE, graphics=fitz.PDF_REDACT_LINE_ART_NONE)
for k, s in hits.items():
    c = s["color"]; rgb = ((c>>16)&255)/255, ((c>>8)&255)/255, (c&255)/255
    tw = fitz.TextWriter(page.rect, color=rgb)
    tw.append(s["origin"], vals[k], font=FONT, fontsize=s["size"])
    tw.write_text(page)
doc.save(out_path, garbage=4, deflate=True)
print("wrote", out_path)
