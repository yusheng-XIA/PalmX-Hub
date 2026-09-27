#!/usr/bin/env python3
"""Replace the Fig. 5g dSV-density ring (2,657-dSV set) with the current 1,480-dSV set.

Geometry and encoding follow the original script redraw_circos_60x60mm_20260921.py:
polar axes, theta clockwise from east, 16 chromosomes with a 0.016-rad gap starting at pi/2,
1-Mb bins, bar width 0.96 x bin span, ring base r = 0.45, height 0.070 x count / max(count),
colour #C51B8A. Circle centre and the radius scale are fitted from the existing vector bars.
"""
import math
import sys
from pathlib import Path

import fitz
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
DENS = HERE.parent / "sd_fix/Fig5g_density_1Mb.tsv"
COL = (0xC5 / 255, 0x1B / 255, 0x8A / 255)
BASE, HEIGHT, GAP = 0.45, 0.070, 0.016


def is_dsv(d):
    f = d.get("fill")
    return f is not None and all(abs(a - b) < 0.02 for a, b in zip(f, COL))


def main(pdf_in, pdf_out):
    dens = pd.read_csv(DENS, sep="\t")
    sizes = dens.groupby("Chrom", sort=False).End0.max()
    total = float(sizes.sum())
    usable = 2 * math.pi - GAP * len(sizes)
    off, cur = {}, math.pi / 2
    for c, s in sizes.items():
        span = usable * s / total
        off[c] = (cur, span, s)
        cur += span + GAP

    def theta(c, p):
        a, span, s = off[c]
        return a + span * p / s

    doc = fitz.open(pdf_in)
    page = doc[0]
    bars = [d for d in page.get_drawings() if is_dsv(d) and len(d["items"]) > 1]
    # first point of each bar path lies on the inner arc (r = BASE): fit the circle
    pts = np.array([[d["items"][0][1].x, d["items"][0][1].y] for d in bars])
    A = np.c_[2 * pts, np.ones(len(pts))]
    b = (pts ** 2).sum(1)
    cx, cy, c0 = np.linalg.lstsq(A, b, rcond=None)[0]
    R0 = math.sqrt(c0 + cx ** 2 + cy ** 2)
    k = R0 / BASE                                   # pt per data-radius unit
    res = np.abs(np.hypot(pts[:, 0] - cx, pts[:, 1] - cy) - R0)
    print(f"centre ({cx:.2f},{cy:.2f}) R0={R0:.3f} k={k:.2f} fit max residual {res.max():.3f} pt, n bars {len(bars)}")

    # validation: angular positions/heights of the drawn bars reproduce the 2,657 set
    old = dens.dSV_2657_as_plotted
    om = old.max()
    tops = []
    for d in bars:
        P = np.array([[q.x, q.y] for it in d["items"] for q in it[1:] if isinstance(q, fitz.Point)])
        r = np.hypot(P[:, 0] - cx, P[:, 1] - cy).max() / k
        ang = math.atan2(np.mean(P[:, 1]) - cy, np.mean(P[:, 0]) - cx) % (2 * math.pi)
        tops.append((ang, r))
    exp = sorted(((theta(r.Chrom, (r.Start0 + r.End0) / 2) % (2 * math.pi), BASE + HEIGHT * v / om)
                  for r, v in zip(dens.itertuples(), old) if v > 0))
    obs = sorted(t for t in tops if t[1] > BASE + 0.002)
    ok = len(exp) == len(obs) and max(abs(a[1] - b[1]) for a, b in zip(exp, obs)) < 0.004
    print("validation vs 2,657 set:", len(exp), "expected non-zero bars,", len(obs), "observed; heights match:", ok)
    if not ok:
        sys.exit("geometry validation failed; not replacing")

    # remove the old ring bars from the content stream: each is "q 1 0 0 1 x y cm .773 .106 .541 scn <path> f Q";
    # the legend swatch ("... scn x y w h re f") has no q/cm wrapper and is kept
    import re
    pat = re.compile(rb"q 1 0 0 1 [-\d.]+ [-\d.]+ cm \.773 \.106 \.541 scn [^Qq]*? f Q ?")
    removed = 0
    for xref in page.get_contents():
        st = doc.xref_stream(xref)
        st2, k2 = pat.subn(b"", st)
        if k2:
            doc.update_stream(xref, st2)
            removed += k2
    print("old bar blocks removed from content stream:", removed)
    left = [d for d in page.get_drawings() if is_dsv(d) and len(d["items"]) > 1]
    print("old bars left after removal:", len(left))

    new = dens.dSV_1480
    nm = new.max()
    sh = page.new_shape()
    n = 0
    for r, v in zip(dens.itertuples(), new):
        if v <= 0:
            continue
        t0, t1 = theta(r.Chrom, r.Start0), theta(r.Chrom, r.End0)
        c, w = (t0 + t1) / 2, (t1 - t0) * 0.96
        a0, a1 = c - w / 2, c + w / 2
        ri, ro = BASE * k, (BASE + HEIGHT * v / nm) * k
        angs = np.linspace(a0, a1, 6)
        poly = [fitz.Point(cx + ri * math.cos(a), cy + ri * math.sin(a)) for a in angs]
        poly += [fitz.Point(cx + ro * math.cos(a), cy + ro * math.sin(a)) for a in angs[::-1]]
        sh.draw_polyline(poly + [poly[0]])
        n += 1
    sh.finish(color=None, fill=COL, closePath=True, width=0)
    sh.commit()
    doc.save(pdf_out, garbage=3, deflate=True)
    print("new bars drawn:", n, "max count", nm, "(old max", om, ")")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
