#!/usr/bin/env python3
"""Text corrections on the author's vector flowchart (油棕种质.pdf, 22 Sep version) for ED Fig. 1a.

Each edit unit = a run of spans on one text line. The run is removed with a text-only redaction
(images and line art untouched) and rewritten at the same baseline and font size with embedded
Times New Roman Bold / Bold Italic (matching the original TimesNRMTPro weight/style), either
left-aligned at the original origin or centred on the original run. The in-figure title is removed
and the page is cropped below it. Output: ED1a_study_design.pdf + ED1a_text_edits.tsv.
"""
import re
from pathlib import Path

import fitz

HERE = Path(__file__).resolve().parent
SUP = "/System/Library/Fonts/Supplemental/"
FF = {"b": SUP + "Times New Roman Bold.ttf", "bi": SUP + "Times New Roman Bold Italic.ttf"}
FONT = {k: fitz.Font(fontfile=v) for k, v in FF.items()}

doc = fitz.open(HERE / "src/油棕种质.pdf")
pg = doc[0]
lines = [l for b in pg.get_text("dict")["blocks"] for l in b.get("lines", [])]


def style(span):
    return "bi" if span["flags"] & 2 else "b"


def fix(t):
    t = t.replace("bacjgrounds", "backgrounds").replace("（", "(").replace("）", ")").replace("＋", "+")
    t = re.sub(r"E\.(Oleifera|oleifera)", "E. oleifera", t)
    t = t.replace("E.guineensis", "E. guineensis")
    t = re.sub(r"\b(sh|Sh)([+-])/(sh|Sh)([+-])", lambda m: f"{m[1]}{'+' if m[2]=='+' else '−'}/{m[3]}{'+' if m[4]=='+' else '−'}", t)
    return t


def find(first_text, y_hint):
    for l in lines:
        for i, s in enumerate(l["spans"]):
            if s["text"] == first_text and abs(s["origin"][1] - y_hint) < 20:
                return l, i
    raise KeyError(first_text)


# (first span text, approx baseline y, number of spans in run (None = to end of line), align, manual overrides)
UNITS = [
    ("sh-/sh-", 245, 1, "left", None), ("Sh+/sh-", 243, 1, "left", None),
    ("E.oleifera × E.guineensis", 366, 1, "center", None),
    ("E.guineensis", 500, None, "center", None),
    ("sh-/sh-", 647, 1, "left", None),
    ("（", 648, None, "center", None),
    ("Sh+/sh-", 662, 1, "left", None),
    ("E.guineensis", 647, None, "center", None),
    ("Sh+/sh-", 663, 1, "left", None),
    ("Oleifera", 618, 1, "center", [("E. oleifera", "bi")]),
    (" E.oleifera", 634, 1, "center", [("American oil palm", "b")]),
    ("E.oleifera × E.guineensis", 637, 1, "center", None),
    ("＋", 819, 1, "center", [("+", "b")]),
    ("＋", 819, 1, "center2", [("+", "b")]),
    ("TK, NS, TN, Nigerian, Oleifera and FL are the six accessions with phased assemblies", 937, 1, "left",
     [("TK, NS, TN, Nigerian, ", "b"), ("E. oleifera", "bi"), (" and FL: accessions with phased assemblies", "b")]),
    (" E.Oleifera", 983, None, "left", None),
    ("Breeding bacjgrounds (Deli, AVROS) and geographical origins are independent of SHELL fruit forms", 979, 1, "left", None),
]

plans, used = [], set()
for first, y, n, align, override in UNITS:
    # the two '＋' spans: take them in order
    cands = [(l, i) for l in lines for i, s in enumerate(l["spans"])
             if s["text"] == first and abs(s["origin"][1] - y) < 20 and (id(l), i) not in used]
    l, i = cands[0]
    used.add((id(l), i))
    spans = l["spans"][i:(i + n) if n else None]
    old = "".join(s["text"] for s in spans)
    new = override or [(fix(s["text"]), style(s) if "Noto" not in s["font"] else "b") for s in spans]
    size = spans[0]["size"]
    x0 = min(s["bbox"][0] for s in spans); x1 = max(s["bbox"][2] for s in spans)
    y0 = min(s["bbox"][1] for s in spans); y1 = max(s["bbox"][3] for s in spans)
    base = spans[-1]["origin"][1] if "Noto" in spans[0]["font"] else spans[0]["origin"][1]
    if "Noto" in spans[0]["font"] and len(spans) > 1:
        base = spans[1]["origin"][1]
    col = spans[0]["color"]
    plans.append(dict(old=old, new=new, size=size, rect=(x0, y0, x1, y1), base=base,
                      origin=spans[0]["origin"][0], align=align, color=col))
    h = y1 - y0
    pg.add_redact_annot(fitz.Rect(x0 + 0.3, y0 + 0.3 * h, x1 - 0.3, y1 - 0.3 * h))

# remove the in-figure title
title = [s for l in lines for s in l["spans"] if "Oil Palm Diversity" in s["text"]][0]
pg.add_redact_annot(fitz.Rect(title["bbox"]))
# hidden leftover text "SHELL genotype" (covered by a white box in the source; would surface after editing)
hid = [s for l in lines for s in l["spans"] if s["text"].strip() == "SHELL genotype"]
for s in hid:
    pg.add_redact_annot(fitz.Rect(s["bbox"]))
print("hidden spans removed:", len(hid))
pg.apply_redactions(images=fitz.PDF_REDACT_IMAGE_NONE, graphics=fitz.PDF_REDACT_LINE_ART_NONE)

for k, v in FF.items():
    pg.insert_font(fontname="TNR" + k, fontfile=v)
rows = []
for p in plans:
    segs = [(t, f) for t, f in p["new"] if t]
    w = sum(FONT[f].text_length(t, p["size"]) for t, f in segs)
    x0, y0, x1, y1 = p["rect"]
    if p["align"] == "left":
        x = p["origin"]
    else:
        x = (x0 + x1) / 2 - w / 2
    for t, f in segs:
        c = p["color"]
        rgb = ((c >> 16) & 255) / 255, ((c >> 8) & 255) / 255, (c & 255) / 255
        pg.insert_text((x, p["base"]), t, fontname="TNR" + f, fontsize=p["size"], color=rgb)
        x += FONT[f].text_length(t, p["size"])
    new_s = "".join(t for t, _ in segs)
    rows.append((round(x0, 1), round(y0, 1), round(x1, 1), round(y1, 1), p["size"], p["old"], new_s))

# crop the title band: top of the first panel frame
pg.set_cropbox(fitz.Rect(0, title["bbox"][3] + 8, pg.rect.width, pg.rect.height))
# drop the now-unused NotoSansCJK font resource (full-width glyphs were all replaced)
cont = b"".join(doc.xref_stream(x) for x in pg.get_contents())
keep = []
for f in pg.get_fonts():
    if "Noto" in f[3]:
        assert f"/{f[4]} ".encode() not in cont, "Noto font still used"
    else:
        keep.append(f"/{f[4]} {f[0]} 0 R")
doc.xref_set_key(pg.xref, "Resources/Font", "<<" + " ".join(keep) + ">>")
doc.subset_fonts()
doc.save(HERE / "ED1a_study_design.pdf", garbage=4, deflate=True)
with open(HERE / "ED1a_text_edits.tsv", "w") as fh:
    fh.write("x0_pt\ty0_pt\tx1_pt\ty1_pt\tfont_pt\told\tnew\n")
    for r in rows:
        fh.write("\t".join(map(str, r)) + "\n")
for r in rows:
    print(r[5], "->", r[6])
