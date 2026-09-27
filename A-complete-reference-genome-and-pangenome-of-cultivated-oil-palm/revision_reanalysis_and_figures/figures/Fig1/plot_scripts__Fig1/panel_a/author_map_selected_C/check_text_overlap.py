#!/usr/bin/env python3
"""Render a figure-producing script and report overlapping text bounding boxes.
Usage: python3 check_text_overlap.py plot_script.py
The script must leave a `fig` variable in its namespace. savefig calls are neutralised."""
import sys, re
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.text import Text, Annotation

src = open(sys.argv[1]).read()
src = re.sub(r"fig\.savefig\([^\n]*\)", "pass", src)
ns = {}
exec(compile(src, sys.argv[1], "exec"), ns)
fig = eval(sys.argv[2], ns) if len(sys.argv) > 2 else ns["fig"]
fig.canvas.draw()
r = fig.canvas.get_renderer()

texts = []
def add(t, tag):
    if not isinstance(t, Text) or not t.get_visible():
        return
    s = t.get_text().strip()
    if not s:
        return
    bb = Text.get_window_extent(t, renderer=r)  # text box only (Annotation would include its arrow)
    if bb.width == 0 or bb.height == 0:
        return
    texts.append((bb, tag, s.replace("\n", " / ")))

for t in fig.texts:
    add(t, "fig")
for ax in fig.axes:
    for t in ax.texts:
        add(t, "ax")
    if ax.axison:
        for t in ax.get_xticklabels():
            if ax.xaxis.get_visible() and ax.get_xaxis().get_major_ticks() and any(tk.label1.get_visible() for tk in ax.xaxis.get_major_ticks()): add(t, "tick")
        for t in ax.get_yticklabels():
            if ax.yaxis.get_visible() and any(tk.label1.get_visible() for tk in ax.yaxis.get_major_ticks()): add(t, "tick")
    add(ax.title, "title"); add(ax.xaxis.label, "xlabel"); add(ax.yaxis.label, "ylabel")
    if ax.get_legend():
        for t in ax.get_legend().get_texts():
            add(t, "legend")
for lg in fig.legends:
    for t in lg.get_texts():
        add(t, "legend")

def overlap(a, b, pad=0.0):
    return not (a.x1 + pad < b.x0 or b.x1 + pad < a.x0 or a.y1 + pad < b.y0 or b.y1 + pad < a.y0)

pairs = []
for i in range(len(texts)):
    for j in range(i + 1, len(texts)):
        if overlap(texts[i][0], texts[j][0]):
            a, b = texts[i][0], texts[j][0]
            ix = max(0, min(a.x1, b.x1) - max(a.x0, b.x0)); iy = max(0, min(a.y1, b.y1) - max(a.y0, b.y0))
            pairs.append((ix * iy, texts[i], texts[j]))
pairs.sort(key=lambda t: t[0], reverse=True)
W, H = fig.get_size_inches() * fig.dpi
print(f"{len(texts)} text objects; {len(pairs)} overlapping pairs")
for area, (ba, ta, sa), (bb, tb, sb) in pairs:
    print(f"  {area:7.0f}px²  [{ta}] '{sa}'  x  [{tb}] '{sb}'   at ({ba.x0/W:.3f},{ba.y0/H:.3f})")
# text overlapping high-zorder patches (scale bars, badges) in the same axes
from matplotlib.patches import Rectangle
for ax in fig.axes:
    rects = [pt for pt in ax.patches if isinstance(pt, Rectangle) and pt.get_zorder() >= 8]
    for t in ax.texts:
        if not t.get_visible() or not t.get_text().strip():
            continue
        tb = Text.get_window_extent(t, renderer=r)
        for pt in rects:
            pb = pt.get_window_extent(renderer=r)
            if overlap(tb, pb):
                print(f"  text on patch: '{t.get_text().strip()}' overlaps a scale-bar/patch rectangle")
# text crossing the frame of its own axes (would sit on / be cut by the spine)
for ax in fig.axes:
    ab = ax.get_window_extent(renderer=r)
    if not any(sp.get_visible() for sp in ax.spines.values()):
        continue
    for t in ax.texts:
        if not t.get_visible() or not t.get_text().strip():
            continue
        bb = Text.get_window_extent(t, renderer=r)
        inside = bb.x0 >= ab.x0 and bb.x1 <= ab.x1 and bb.y0 >= ab.y0 and bb.y1 <= ab.y1
        outside = bb.x1 < ab.x0 or bb.x0 > ab.x1 or bb.y1 < ab.y0 or bb.y0 > ab.y1
        if not inside and not outside:
            print(f"  text crosses axes frame: '{t.get_text().strip()}'  (text y0={bb.y0:.0f}, axes y0={ab.y0:.0f})")
# also: text outside figure
out = [(t, s) for bb, t, s in texts if bb.x0 < 0 or bb.y0 < 0 or bb.x1 > W or bb.y1 > H]
if out:
    print("text clipped by figure edge:")
    for t, s in out:
        print(f"  [{t}] '{s}'")
