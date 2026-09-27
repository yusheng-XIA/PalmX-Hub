#!/usr/bin/env python3
"""Minimal vector text edits to Figure 2 (PyMuPDF): remove the original text span and
rewrite the new text at the original baseline, size and weight in embedded Arial.

usage: fix_fig2.py SRC.pdf OUT.pdf [--fix-node20]
"""
import re
import sys

import fitz

SRC, OUT = sys.argv[1], sys.argv[2]
FIX_NODE = "--fix-node20" in sys.argv
AR = "/System/Library/Fonts/Supplemental/Arial.ttf"
ARB = "/System/Library/Fonts/Supplemental/Arial Bold.ttf"
MINUS = "−"

doc = fitz.open(SRC)
pg = doc[0]
fR, fB = fitz.Font(fontfile=AR), fitz.Font(fontfile=ARB)

spans = []
for b in pg.get_text("dict")["blocks"]:
    for l in b.get("lines", []):
        for s in l["spans"]:
            spans.append(s)


def find(text, near=None, font=None):
    out = []
    for s in spans:
        if s["text"] == text and (font is None or font in s["font"]):
            if near is None or (abs(s["bbox"][0] - near[0]) < 2 and abs(s["bbox"][1] - near[1]) < 2):
                out.append(s)
    if not out:
        raise SystemExit(f"span not found: {text!r} {near}")
    return out


edits = []   # (span_or_spans, new_text, anchor, bold, size, colour)
log = []


def rewrite(sp, new, anchor="left", bold=None, x=None, size=None):
    sp = sp if isinstance(sp, list) else [sp]
    s0 = sp[0]
    isbold = ("Bold" in s0["font"]) if bold is None else bold
    sz = size or s0["size"]
    font = fB if isbold else fR
    w_new = font.text_length(new, fontsize=sz)
    x0, y = s0["origin"]
    if x is not None:
        x0 = x
    elif anchor == "right":
        x0 = max(s["bbox"][2] for s in sp) - w_new
    elif anchor == "center":
        c = (min(s["bbox"][0] for s in sp) + max(s["bbox"][2] for s in sp)) / 2
        x0 = c - w_new / 2
    col = s0["color"]
    rgb = ((col >> 16) & 255) / 255, ((col >> 8) & 255) / 255, (col & 255) / 255
    for s in sp:
        # thin rectangle through the glyph centres only, so neighbouring labels are untouched
        ox, oy = s["origin"]
        yc = oy - 0.33 * s["size"]
        pg.add_redact_annot(fitz.Rect(s["bbox"][0] + 0.3, yc - 0.12, s["bbox"][2] - 0.3, yc + 0.12), fill=False)
    edits.append((x0, y, new, "ArialB" if isbold else "ArialR", sz, rgb))
    log.append((" + ".join(repr(s["text"]) for s in sp), new))


# ---- 2a: full-width tilde and parentheses (AdobeSongStd) -> ASCII in Arial Bold
for s in [s for s in spans if s["text"] == "～"]:
    rewrite(s, "~", bold=True, anchor="right")
rewrite(find("Ancestral Palmae karyotype"), "Ancestral palm karyotype")
rewrite([find("（", near=(44.0, 18.2))[0], find("APK  n= 5")[0], find("）", near=(84.6, 18.2))[0]],
        "(APK, n = 5)", anchor="center", bold=False, size=7.91)
rewrite(find("Ancestral Core Palmae karyotype"), "Ancestral core palm karyotype")
rewrite([find("（", near=(42.1, 129.1))[0], find("ACPK  n= 10")[0], find(" ）", near=(94.9, 129.1))[0]],
        "(ACPK, n = 10)", anchor="center", bold=False, size=7.91)

# ---- 2b: tip label, true minus signs, optional node-20 correction
rewrite(find("Oleifera", near=(158.8, 229.3)), "FL-Hap1", bold=True)
# internal-node labels that are mis-placed relative to CAFE5 Gamma_asr.tre / Gamma_clade_results.txt
# (only rewritten with --fix-node20; see Figure2_changes.md)
NODE_FIX = {  # (x, y) of '+' span -> correct (+, -) values
    (110.8, 297.6): ("+1412", MINUS + "1680"),   # Calamus + Daemonorops MRCA = node 15
    (38.6, 280.6): ("+56", MINUS + "2302"),      # Calamoideae + (Nypa ...) = node 18
    (10.6, 301.6): ("+13", MINUS + "237"),       # palms + Musa = node 20
}
fix_spans = set()
if FIX_NODE:
    for (x, y), (plus, minus) in NODE_FIX.items():
        sp_plus = find(next(t for t in ("+89", "+1412", "+56") if any(
            abs(s["bbox"][0] - x) < 2 and abs(s["bbox"][1] - y) < 2 and s["text"] == t for s in spans)), near=(x, y))[0]
        sp_minus = [s for s in spans if abs(s["bbox"][0] - x) < 2 and 6 < s["bbox"][1] - y < 11][0]
        rewrite(sp_plus, plus)
        rewrite(sp_minus, minus)
        fix_spans.add(id(sp_minus))
for s in spans:
    t = s["text"]
    x0, y0 = s["bbox"][:2]
    in_2b = 0 < x0 < 265 and 195 < y0 < 360
    if in_2b and re.search(r"(^|\s)-\d", t) and id(s) not in fix_spans:
        rewrite(s, t.strip().replace("-", MINUS), x=(s["origin"][0] + (fR.text_length(" ", fontsize=s["size"]) if t.startswith(" ") else 0)))

# ---- 2d: axis tick labels
for s in spans:
    if s["text"] in ("-0.4", "-0.5") and 270 < s["bbox"][0] < 410:
        rewrite(s, s["text"].replace("-", MINUS), anchor="right")

# ---- 2f
rewrite(find("(1400/6836)"), "(1403/6836)", anchor="center")
rewrite(find("89.1% FL"), "89.1% FL (185 d)", anchor="center")

pg.apply_redactions(images=fitz.PDF_REDACT_IMAGE_NONE, graphics=fitz.PDF_REDACT_LINE_ART_NONE,
                    text=fitz.PDF_REDACT_TEXT_REMOVE)
for x0, y, new, fn, sz, rgb in edits:
    pg.insert_text((x0, y), new, fontname=fn, fontfile=(ARB if fn == "ArialB" else AR), fontsize=sz, color=rgb)

# drop the now-unused CJK font from the page's /Font resource dictionary
for xref, ext, typ, base, name, enc, *_ in pg.get_fonts(full=True):
    if "Song" in base:
        typ_, val = doc.xref_get_key(pg.xref, "Resources/Font")
        if typ_ == "xref":
            fx = int(val.split()[0])
            src = doc.xref_object(fx, compressed=True)
            doc.update_object(fx, re.sub(r"/%s\s+\d+\s+0\s+R" % re.escape(name), "", src))
        elif typ_ == "dict":
            doc.xref_set_key(pg.xref, "Resources/Font", re.sub(r"/%s\s+\d+\s+0\s+R" % re.escape(name), "", val))
        else:
            rt, rv = doc.xref_get_key(pg.xref, "Resources")
            rx = int(rv.split()[0]) if rt == "xref" else pg.xref
            ft, fv = doc.xref_get_key(rx, "Font")
            if ft == "xref":
                fx = int(fv.split()[0])
                doc.update_object(fx, re.sub(r"/%s\s+\d+\s+0\s+R" % re.escape(name), "", doc.xref_object(fx, compressed=True)))
            else:
                doc.xref_set_key(rx, "Font", re.sub(r"/%s\s+\d+\s+0\s+R" % re.escape(name), "", fv))
doc.subset_fonts()
# Arial maps U+002D and U+00AD to the same glyph (GID 16); make the ToUnicode entry a normal hyphen
for fx, *_rest in pg.get_fonts(full=True):
    base = _rest[2]
    if base.startswith("Arial") and "+" not in base:
        tu = doc.xref_get_key(fx, "ToUnicode")
        if tu[0] == "xref":
            tx = int(tu[1].split()[0])
            cmap = doc.xref_stream(tx).decode("latin1").replace("<0010> <00ad>", "<0010> <002d>")
            doc.update_stream(tx, cmap.encode("latin1"))
doc.save(OUT, garbage=4, deflate=True)
for a, b in log:
    print(f"{a} -> {b!r}")
