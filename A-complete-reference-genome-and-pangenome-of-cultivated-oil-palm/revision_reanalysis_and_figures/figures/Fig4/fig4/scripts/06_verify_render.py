#!/usr/bin/env python3
"""Checks on the candidate + 300 dpi before/after images.
- page size unchanged; fonts Arial family only
- words: full-page multiset diff (orig vs candidate) must equal the intended edits
- vector drawings of panels a, e, f, g, h, j identical (rect, colours, width, dashes)
- 600 dpi pixel diff: changed pixels only inside panels b, c, d, i, k
"""
import collections, json
import numpy as np, fitz
from PIL import Image, ImageDraw, ImageFont
from common import *

CAND = FIX / "Figure4_candidate.pdf"
PAN = {"a": (0, 0, 172, 112), "b": (172, 0, 343, 116), "c": (343, 0, 519, 116), "d": (0, 112, 172, 228),
       "e": (172, 112, 343, 228), "f": (343, 112, 519, 233), "g": (0, 228, 172, 348), "h": (172, 228, 343, 348),
       "i": (343, 233, 519, 348), "j": (0, 348, 226, 621), "k": (226, 348, 519, 621)}
o, c = fitz.open(ORIG)[0], fitz.open(CAND)[0]
rep = {}
rep["page_size"] = [list(o.rect), list(c.rect)]
assert o.rect == c.rect
rep["fonts"] = sorted({f[3] for f in c.get_fonts(full=True)})
assert all("Arial" in f for f in rep["fonts"])

wo = collections.Counter(w[4] for w in o.get_text("words"))
wc = collections.Counter(w[4] for w in c.get_text("words"))
rep["words_removed"] = dict(wo - wc); rep["words_added"] = dict(wc - wo)


def sig(page, r):
    R = fitz.Rect(r)
    out = []
    for d in page.get_drawings():
        if R.contains(d["rect"]):
            out.append((d["type"], tuple(round(v, 2) for v in d["rect"]), d.get("fill"), d.get("color"),
                        d.get("width"), d.get("dashes"), len(d["items"])))
    return sorted(out, key=str)


rep["drawings_identical"] = {}
for p in "aefghj":
    rep["drawings_identical"][p] = sig(o, PAN[p]) == sig(c, PAN[p])
for p in "bc":   # b, c: graphics must be identical too (text-only edits)
    rep["drawings_identical"][p] = sig(o, PAN[p]) == sig(c, PAN[p])

po = o.get_pixmap(dpi=600, alpha=False); pc = c.get_pixmap(dpi=600, alpha=False)
A = np.frombuffer(po.samples, np.uint8).reshape(po.height, po.width, 3)
B = np.frombuffer(pc.samples, np.uint8).reshape(pc.height, pc.width, 3)
diff = np.abs(A.astype(int) - B.astype(int)).max(axis=2) > 8
sc = 600 / 72
rep["changed_px_by_panel"] = {}
for p, (x0, y0, x1, y1) in PAN.items():
    sub = diff[int(y0 * sc):int(y1 * sc), int(x0 * sc):int(x1 * sc)]
    n = int(sub.sum())
    if n:
        ys, xs = np.nonzero(sub)
        rep["changed_px_by_panel"][p] = [n, [round(x0 + xs.min() / sc, 1), round(y0 + ys.min() / sc, 1),
                                             round(x0 + xs.max() / sc, 1), round(y0 + ys.max() / sc, 1)]]
    else:
        rep["changed_px_by_panel"][p] = 0
Image.fromarray(B).save(FIX / "compare/Figure4_candidate_600dpi.png")

# side-by-side 300 dpi before/after per edited panel
font = ImageFont.truetype("/System/Library/Fonts/Supplemental/Arial.ttf", 36)
for p in "bcdik":
    r = fitz.Rect(PAN[p])
    im = [Image.frombytes("RGB", (x.width, x.height), x.samples) for x in
          (o.get_pixmap(dpi=300, clip=r, alpha=False), c.get_pixmap(dpi=300, clip=r, alpha=False))]
    W, H = im[0].width, im[0].height
    canvas = Image.new("RGB", (2 * W + 40, H + 60), "white")
    dr = ImageDraw.Draw(canvas)
    for k, (lab, img) in enumerate((("before (submitted)", im[0]), ("after (candidate)", im[1]))):
        canvas.paste(img, (k * (W + 40), 60)); dr.text((k * (W + 40) + 10, 10), f"Fig. 4{p} - {lab}", fill="black", font=font)
    canvas.save(FIX / f"compare/Fig4{p}_before_after_300dpi.png", dpi=(300, 300))
c.parent.get_page_pixmap if False else None
pv = c.get_pixmap(dpi=150, alpha=False); pv.save(FIX / "compare/Figure4_candidate_preview_150dpi.png")
json.dump(rep, open(FIX / "compare/verify_report.json", "w"), indent=1, default=str)
print(json.dumps(rep, indent=1, default=str))
