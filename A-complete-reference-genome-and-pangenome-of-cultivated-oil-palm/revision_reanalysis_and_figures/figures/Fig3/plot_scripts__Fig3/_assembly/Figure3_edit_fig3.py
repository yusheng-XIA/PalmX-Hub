#!/usr/bin/env python3
"""Figure 3: vector text corrections (PyMuPDF). Removes the original glyphs by redaction
(text only; line art and images untouched) and rewrites at the original baseline, size and colour."""
import fitz

SRC, OUT = "src.pdf", "edited.pdf"
AR = "/System/Library/Fonts/Supplemental/Arial.ttf"
ARI = "/System/Library/Fonts/Supplemental/Arial Italic.ttf"
fR, fI = fitz.Font(fontfile=AR), fitz.Font(fontfile=ARI)

doc = fitz.open(SRC)
page = doc[0]
spans = [s for b in page.get_text("dict")["blocks"] for l in b.get("lines", []) for s in l["spans"]]


def find(text, near=None):
    c = [s for s in spans if s["text"].strip() == text]
    if near:
        c = [s for s in c if abs(s["bbox"][0] - near[0]) < 3 and abs(s["bbox"][1] - near[1]) < 3]
    assert len(c) == 1, (text, len(c))
    return c[0]


def rgb(c):
    return ((c >> 16) & 255) / 255, ((c >> 8) & 255) / 255, (c & 255) / 255


edits = []  # (span, runs[(text, font, raise)], align, size)
log = []


def plan(span, runs, align="left", size=None, why=""):
    edits.append((span, runs, align, size or span["size"]))
    log.append((span["text"].strip(), "".join(r[0] for r in runs), [round(v, 1) for v in span["bbox"]], why))


# 3h header
plan(find("Seed development"), [("Development", "R", 0)], "center", why="mesocarp fruit development, not seed")
# 3l legend
plan(find("Favorable action (293 targets)"), [("Favourable action (293 targets)", "R", 0)], "left", why="British spelling (legend text)")
plan(find("Oleifera"), [("E. oleifera", "I", 0)], "left", why="species name")
plan(find("O = Oleifera"), [("O = ", "R", 0), ("E. oleifera", "I", 0)], "left", why="species name")
for k in range(1, 10):
    plan(find(f"chr{k}"), [(f"chr{k:02d}", "R", 0)], "right", why="chromosome naming as elsewhere (chr01)")
# 3e P values (Source Data Fig.3e_tests: KW P = 2.160e-14 and 1.273e-13)
plan(find("P = 2.2e-14"), [("P = 2.2 × 10", "R", 0), ("−14", "R", 2.4)], "center", why="KW P = 2.16e-14")
plan(find("P = 1.3e-13"), [("P = 1.3 × 10", "R", 0), ("−13", "R", 2.4)], "center", why="KW P = 1.27e-13")
# 3c 5-pt annotations -> 6.15 pt
# 3c: A/B follow the ASE pipeline definition (FL A = FL-Hap2, B = FL-Hap1; TN A = TK-like, B = NS-like)
plan(find("Hap 1 > Hap 2"), [("A > B", "R", 0)], "right", size=6.15, why="A/B per ASE definition; enlarged from 5 pt")
plan(find("Hap 1 < Hap 2"), [("A < B", "R", 0)], "right", size=6.15, why="A/B per ASE definition; enlarged from 5 pt")
for old_t, new_t in [("TN_hap1 (Hap1> Hap2)", "TN A (TK-like)"), ("FL_hap1 (Hap1> Hap2)", "FL A (FL-Hap2)"),
                     ("TN_hap2 (Hap1 < Hap2)", "TN B (NS-like)"), ("FL_hap2 (Hap1 < Hap2)", "FL B (FL-Hap1)")]:
    plan(find(old_t), [(new_t, "R", 0)], "left", why="3c legend: A/B per ASE data definition (author decision)")
# 3l x-axis tick labels (6.1 pt) -> 6.3 pt
for s in spans:
    if 6.0 < s["size"] < 6.15 and 630 < s["bbox"][1] < 634:
        plan(s, [(s["text"].strip(), "R", 0)], "center", size=6.3, why="tick label enlarged from 6.1 pt")

# redact original glyphs (text only)
for s, runs, align, size in edits:
    x0, y0, x1, y1 = s["bbox"]
    page.add_redact_annot(fitz.Rect(x0 + 0.3, y0 + 1.5, x1 - 0.3, y1 - 1.5), fill=None)
page.apply_redactions(images=fitz.PDF_REDACT_IMAGE_NONE, graphics=fitz.PDF_REDACT_LINE_ART_NONE,
                      text=fitz.PDF_REDACT_TEXT_REMOVE)

page.insert_font(fontname="ArR", fontfile=AR)
page.insert_font(fontname="ArI", fontfile=ARI)
for s, runs, align, size in edits:
    x0, y0, x1, y1 = s["bbox"]
    base = s["origin"][1]
    widths = [(fI if f == "I" else fR).text_length(t, fontsize=size) for t, f, _ in runs]
    total = sum(widths)
    if align == "left":
        x = x0
    elif align == "right":
        x = s["origin"][0] + fR.text_length(s["text"].rstrip(), fontsize=s["size"]) - total
    else:
        x = (x0 + x1) / 2 - total / 2
    for (t, f, rise), w in zip(runs, widths):
        page.insert_text((x, base - rise), t, fontname="ArI" if f == "I" else "ArR", fontsize=size,
                         color=rgb(s["color"]))
        x += w

doc.save(OUT, garbage=3, deflate=True)
with open("text_edits.tsv", "w") as fh:
    fh.write("original\tnew\tbbox_pt\treason\n")
    for o, n, b, w in log:
        fh.write(f"{o}\t{n}\t{b}\t{w}\n")
print(len(edits), "edits")
