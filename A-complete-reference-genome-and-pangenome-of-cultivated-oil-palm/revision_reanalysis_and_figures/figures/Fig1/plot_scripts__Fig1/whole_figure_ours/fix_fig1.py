#!/usr/bin/env python3
"""Figure 1: text/naming fixes and resize to the Nature page area, all in vector form.

1. Text edits: the original span is redacted (text only; drawings and images are kept) and
   rewritten at the same baseline, size and colour in Arial (the figure's own typeface).
2. Panel b material key: dots and labels in the key row are removed and redrawn in the
   same colours with longer labels (E. oleifera / FL / TK (dura) / NS (pisifera) / Nigerian / TN).
3. Resize: crop to the ink bounding box, then scale uniformly to 247 mm height.
"""
from pathlib import Path

import fitz

HERE = Path(__file__).resolve().parent
SRC = HERE / "src.pdf"
OUT = HERE.parent / "deliver" / "Main_Figures_revised"
OUT.mkdir(parents=True, exist_ok=True)

FONTS = {
    "R": "/System/Library/Fonts/Supplemental/Arial.ttf",
    "I": "/System/Library/Fonts/Supplemental/Arial Italic.ttf",
    "B": "/System/Library/Fonts/Supplemental/Arial Bold.ttf",
}
FOBJ = {k: fitz.Font(fontfile=v) for k, v in FONTS.items()}

# (old text, [(style, new text), ...]) ; style R = roman, I = italic
EDITS = {
    # panel b haplotype rows
    "Nigerian_hap1": [("R", "Nigerian-Hap1")], "Nigerian_hap2": [("R", "Nigerian-Hap2")],
    "Dura_hap1": [("R", "TK-Hap1")], "Dura_hap2": [("R", "TK-Hap2")],
    "Pisifera_hap1": [("R", "NS-Hap1")], "Pisifera_hap2": [("R", "NS-Hap2")],
    "TN_hap1": [("R", "TN-Hap1")], "TN_hap2": [("R", "TN-Hap2")],
    "FL_hap1": [("R", "FL-Hap1")], "FL_hap2": [("R", "FL-Hap2")],
    "Oleifera_hap1": [("I", "E. oleifera"), ("R", "-Hap1")],
    "Oleifera_hap2": [("I", "E. oleifera"), ("R", "-Hap2")],
    # panel b haplotype symbol key
    "hap1": [("R", "Hap1")], "hap2": [("R", "Hap2")],
    # panel b telomere column: reference values without a traceable source (removed on author decision, 2026-09-23)
    "8/32": [("R", "–")], "13/32": [("R", "–")],
    # panel c key
    "Genes density": [("R", "Gene density")],
    "FL exclusive gene density": [("R", "FL-specific gene density")],
    "Telemeres": [("R", "Telomeres")],
    # panel e ancestry key (wording as in the figure legend)
    "Dura-like": [("R", "dura-like")],
    "Pisifera-like": [("R", "pisifera-like")],
    "Oleifera-like": [("I", "E. oleifera"), ("R", "-like")],
}
# panel e group labels (only the italic row labels at x ≈ 295)
GROUP_EDITS = {"Dura": [("R", "TK")], "Pisifera": [("R", "NS")], "Nigerian": [("R", "Nigerian")]}

KEY_Y = (189.9, 197.2)            # panel b material key row
KEY_X = (312.0, 480.0)
KEY_LABELS = [[("I", "E. oleifera")], [("R", "FL")], [("R", "TK (dura)")],
              [("R", "NS (pisifera)")], [("R", "Nigerian")], [("R", "TN")]]


def rgb(c):
    return ((c >> 16) & 255) / 255, ((c >> 8) & 255) / 255, (c & 255) / 255


def write(page, origin, parts, size, color):
    x, y = origin
    for style, txt in parts:
        page.insert_text((x, y), txt, fontsize=size, fontname=f"Arial{style}", fontfile=FONTS[style],
                         color=color)
        x += FOBJ[style].text_length(txt, fontsize=size)
    return x


doc = fitz.open(SRC)
page = doc[0]
log = []

# collect spans to rewrite
todo = []
for b in page.get_text("dict")["blocks"]:
    for ln in b.get("lines", []):
        for s in ln["spans"]:
            t = s["text"].strip()
            if t in EDITS:
                todo.append((s, EDITS[t]))
            elif t in GROUP_EDITS and "Italic" in s["font"] and abs(s["origin"][0] - 295.3) < 1:
                todo.append((s, GROUP_EDITS[t]))

# material key row: remember dot colours / geometry, then remove dots and labels
dots = [dr for dr in page.get_drawings()
        if KEY_Y[0] <= dr["rect"].y0 and dr["rect"].y1 <= KEY_Y[1] and KEY_X[0] <= dr["rect"].x0 <= KEY_X[1]
        and dr.get("fill") is not None]
dots.sort(key=lambda d: d["rect"].x0)
assert len(dots) == 6, len(dots)
key_spans = [s for b in page.get_text("dict")["blocks"] for ln in b.get("lines", []) for s in ln["spans"]
             if 188 < s["bbox"][1] < 190 and 315 < s["bbox"][0] < 480 and s["text"].strip()]
key_origin_y = key_spans[0]["origin"][1]
key_color = rgb(key_spans[0]["color"])

for s, _ in todo:
    # stop the box just above the baseline so that text on the next line is never touched
    x0, y0, x1, y1 = s["bbox"]
    page.add_redact_annot(fitz.Rect(x0, y0, x1, min(y1, s["origin"][1] - 0.3)))
page.apply_redactions(images=fitz.PDF_REDACT_IMAGE_NONE, graphics=fitz.PDF_REDACT_LINE_ART_NONE)
page.add_redact_annot(fitz.Rect(KEY_X[0], KEY_Y[0], KEY_X[1], KEY_Y[1]))
page.apply_redactions(images=fitz.PDF_REDACT_IMAGE_NONE,
                      graphics=fitz.PDF_REDACT_LINE_ART_REMOVE_IF_COVERED)

ROW_SIZE = 6.0   # panel b row labels: 6.5 pt -> 6.0 pt so that "E. oleifera-Hap1" stays inside the row box
for s, parts in todo:
    size = ROW_SIZE if s["text"].strip().endswith(("_hap1", "_hap2")) else s["size"]
    write(page, s["origin"], parts, size, rgb(s["color"]))
    log.append((s["text"].strip(), "".join(p[1] for p in parts), round(s["origin"][0], 1), round(s["origin"][1], 1)))

# redraw the key
x = dots[0]["rect"].x0
w, h = dots[0]["rect"].width, dots[0]["rect"].height
y0 = dots[0]["rect"].y0
for dot, label in zip(dots, KEY_LABELS):
    page.draw_oval(fitz.Rect(x, y0, x + w, y0 + h), color=None, fill=dot["fill"])
    end = write(page, (x + 5.62, key_origin_y), label, 6.0, key_color)
    x = end + 4.6
assert x < 518, x
log.append(("material key", "E. oleifera / FL / TK (dura) / NS (pisifera) / Nigerian / TN", 314.0, 192.8))

# resize: crop to ink box, scale to 247 mm height
ink = fitz.Rect(0, 3.8, 518.74, 758.2)
H = 247 / 25.4 * 72
s = H / ink.height
W = ink.width * s
out = fitz.open()
np_ = out.new_page(width=W, height=H)
np_.show_pdf_page(np_.rect, doc, 0, clip=ink)
# Arial maps U+002D and U+00AD to the same glyph and the inserted fonts' ToUnicode picks U+00AD;
# point it back to the ASCII hyphen so that the text extracts correctly.
tmp = HERE / "_tmp.pdf"
out.save(tmp, garbage=4, deflate=True)
out = fitz.open(tmp)
n_fix = 0
for xref in range(1, out.xref_length()):
    if out.xref_is_stream(xref):
        st = out.xref_stream(xref)
        if st and b"begincmap" in st and (b"<00AD>" in st or b"<00ad>" in st):
            out.update_stream(xref, st.replace(b"<00AD>", b"<002D>").replace(b"<00ad>", b"<002d>"))
            n_fix += 1
out.save(OUT / "Figure1.pdf", garbage=4, deflate=True)
tmp.unlink()
print("ToUnicode streams fixed:", n_fix)
print(f"scale {s:.4f}; size {W / 72 * 25.4:.1f} x {H / 72 * 25.4:.1f} mm; min font 5.6 pt -> {5.6 * s:.2f} pt")
with open(HERE / "edit_log.tsv", "w") as fh:
    fh.write("original\tnew\tx\ty\n")
    for r in log:
        fh.write("\t".join(map(str, r)) + "\n")
