#!/usr/bin/env python3
"""Targeted text corrections in the ED Fig. 1a study-design raster (流程图.png; no editable source exists).

Word boxes were located by ink segmentation of each text line (not by eye). Each edit paints the
original word box with the local background colour and writes the corrected text in the same
Times New Roman face; the font size is chosen so that the ORIGINAL string reproduces the measured
ink width, and the new text is top-aligned to the original ink. Pixels outside the recorded edit
boxes are verified unchanged. The in-figure title band is cropped.
"""
from pathlib import Path

import numpy as np
from PIL import Image, ImageChops, ImageDraw, ImageFont

HERE = Path(__file__).resolve().parent
Image.MAX_IMAGE_PIXELS = None
SUP = "/System/Library/Fonts/Supplemental/"
F = {"b": SUP + "Times New Roman Bold.ttf", "bi": SUP + "Times New Roman Bold Italic.ttf"}

src = Image.open(HERE / "src/流程图.png").convert("RGB")
img = src.copy()
arr = np.asarray(src).astype(int)
draw = ImageDraw.Draw(img)
edited = []


def ink(x0, x1, y0, y1, thr=120):
    sub = arr[y0:y1, x0:x1]
    bg = np.median(np.concatenate([sub[:4].reshape(-1, 3), sub[-4:].reshape(-1, 3)]), axis=0)
    d = np.abs(sub - bg).sum(2)
    ys, xs = np.where(d > thr)
    box = (x0 + xs.min(), y0 + ys.min(), x0 + xs.max() + 1, y0 + ys.max() + 1)
    col = tuple(int(v) for v in sub.reshape(-1, 3)[d.reshape(-1).argmax()])
    return box, tuple(int(v) for v in bg), col


def ink_width(segs, size):
    x, xs = 0, []
    for t, f in segs:
        ft = ImageFont.truetype(F[f], size)
        bb = ft.getbbox(t)
        xs += [x + bb[0], x + bb[2]]
        x += ft.getlength(t)
    return max(xs) - min(xs)


def fit(segs, w):
    lo, hi = 10, 400
    while hi - lo > 1:
        mid = (lo + hi) // 2
        lo, hi = (mid, hi) if ink_width(segs, mid) <= w else (lo, mid)
    return lo


def replace(x0, x1, y0, y1, old, new, align="left", clip_bottom=None, pad=5):
    ib, bg, col = ink(x0, x1, y0, y1)
    size = fit(old, ib[2] - ib[0])
    bottom = clip_bottom if clip_bottom else ib[3] + pad
    paint = (ib[0] - pad, ib[1] - pad, ib[2] + pad, bottom)
    draw.rectangle((paint[0], paint[1], paint[2] - 1, paint[3] - 1), fill=bg)
    layer = Image.new("RGBA", img.size, (0, 0, 0, 0))
    ld = ImageDraw.Draw(layer)
    top_off = min(ImageFont.truetype(F[f], size).getbbox(t)[1] for t, f in new)
    wn = ink_width(new, size)
    first_l = ImageFont.truetype(F[new[0][1]], size).getbbox(new[0][0])[0]
    x = ib[0] - first_l if align == "left" else (ib[0] + ib[2]) / 2 - wn / 2 - first_l
    y = ib[1] - top_off
    for t, f in new:
        ft = ImageFont.truetype(F[f], size)
        ld.text((x, y), t, font=ft, fill=col + (255,))
        x += ft.getlength(t)
    if clip_bottom:
        a = np.asarray(layer).copy()
        a[clip_bottom:, :, 3] = 0
        layer = Image.fromarray(a)
    bbox = layer.getbbox()
    img.paste(layer, (0, 0), layer)
    box = (min(paint[0], bbox[0]), min(paint[1], bbox[1]), max(paint[2], bbox[2]), max(paint[3], bbox[3]))
    edited.append((box, "".join(t for t, _ in old), "".join(t for t, _ in new), size))


# --- material boxes (row 4) ---
BX, BY = 3252, 4978
replace(BX + 300, BX + 745, BY + 300, BY + 425, [("Sh+/Sh+", "bi")], [("Sh−/Sh−", "bi")], align="center")  # NS
replace(BX + 1515, BX + 1855, 5445, 5513, [("Sh+/Sh+", "bi")], [("Sh+/Sh−", "bi")], align="center",
        clip_bottom=5513)                                                                                    # TN (photo below)
replace(BX + 3740, BX + 4140, BY + 70, BY + 190, [("Oleifera", "b")], [("E. oleifera", "bi")], align="center")
replace(BX + 3710, BX + 4180, BY + 210, BY + 330, [("E.oleifera", "bi")], [("American oil palm", "b")], align="center")
# --- footnotes (word boxes from ink segmentation) ---
replace(685, 1218, 8062, 8213, [("bacjgrounds", "b")], [("backgrounds", "b")])
replace(6371, 9118, 7710, 7830, [("Oleifera and FL are the six accessions with phased assemblies", "b")],
        [("E. oleifera", "bi"), (" and FL: the six accessions with phased assemblies", "b")])
replace(7628, 8774, 8062, 8213, [("E.Oleifera", "bi"), (" × ", "b"), ("E.guineensis", "bi")],
        [("E. oleifera", "bi"), (" × ", "b"), ("E. guineensis", "bi")])

diff = np.asarray(ImageChops.difference(src, img).convert("L")) > 0
mask = np.zeros_like(diff)
for (x0, y0, x1, y1), *_ in edited:
    mask[y0:y1, x0:x1] = True
assert not (diff & ~mask).any(), "pixels changed outside edit boxes"

arr2 = np.asarray(img).astype(int)
frame_rows = [y for y in range(0, 900) if (np.abs(arr2[y, 100:9200] - 255).sum(1) > 60).mean() > 0.6]
top = max(0, frame_rows[0] - 40) if frame_rows else 0
img.crop((0, top, img.size[0], img.size[1])).save(HERE / "ED1a_study_design_corrected.png", dpi=(600, 600))
with open(HERE / "ED1a_text_edits.tsv", "w") as fh:
    fh.write("box_x0\tbox_y0\tbox_x1\tbox_y1\told\tnew\tfont_px\n")
    for (x0, y0, x1, y1), o, n, s in edited:
        fh.write(f"{x0}\t{y0}\t{x1}\t{y1}\t{o}\t{n}\t{s}\n")
print("title crop top =", top)
for e in edited:
    print(e)
