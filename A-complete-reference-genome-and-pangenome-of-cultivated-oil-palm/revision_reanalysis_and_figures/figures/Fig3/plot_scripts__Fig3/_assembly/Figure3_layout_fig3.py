#!/usr/bin/env python3
"""Figure 3 re-layout (vector). 1) remove four blank inter-row bands (15 pt); 2) move the panel-l legend
from below the chromosome plot into the empty right-hand area of the right column (x > 392 pt, rows
chr12-chr16); 3) scale uniformly to <=183 x 247 mm. All placement is show_pdf_page with clips (vector)."""
import fitz

SRC = "edited.pdf"
W = 518.74
KEEP = [(0, 163), (166.5, 299.5), (304, 473.5), (476, 616), (620.5, 870)]   # cut bands lie in all-white rows
H = sum(b - a for a, b in KEEP)


def ymap(y):
    off = 0
    for a, b in KEEP:
        if a <= y <= b:
            return off + (y - a)
        off += b - a
    raise ValueError(y)


src = fitz.open(SRC)
mid = fitz.open()
pg = mid.new_page(width=W, height=H)
off = 0
for a, b in KEEP:
    pg.show_pdf_page(fitz.Rect(0, off, W, off + b - a), src, 0, clip=fitz.Rect(0, a, W, b))
    off += b - a

# legend groups (source rects, measured from rendered pixels, 1-pt margin)
G = {
    "action": fitz.Rect(125.5, 875.8, 234.0, 900.0),
    "chips": fitz.Rect(390.0, 875.8, 477.0, 899.8),
    "ancestry": fitz.Rect(3.0, 875.8, 74.8, 900.0),
    "diploid": fitz.Rect(297.0, 875.8, 350.8, 907.0),
}
x0, y0, gap = 397.0, ymap(731.0), 4.0
box = fitz.Rect(x0 - 4, y0 - 3, W - 0.5, 0)
y = y0
places = []
for k in ("action", "chips", "ancestry", "diploid"):
    r = G[k]
    places.append((k, fitz.Rect(x0, y, x0 + r.width, y + r.height)))
    y += r.height + gap
box.y1 = y - gap + 3
assert box.y1 <= ymap(869.5), box
pg.draw_rect(box, color=(0.85, 0.85, 0.85), fill=(1, 1, 1), width=0.4)
for k, dst in places:
    pg.show_pdf_page(dst, src, 0, clip=G[k])

# uniform scale to the Nature page area
s = min(183 / 25.4 * 72 / W, 247 / 25.4 * 72 / H)
out = fitz.open()
fp = out.new_page(width=W * s, height=H * s)
fp.show_pdf_page(fp.rect, mid, 0)
out.subset_fonts()
out.save("Figure3.pdf", garbage=4, deflate=True)
print(f"content H = {H:.1f} pt; scale = {s:.4f}; final = {W*s/72*25.4:.1f} x {H*s/72*25.4:.1f} mm; "
      f"min font 6.15 pt -> {6.15*s:.2f} pt")
