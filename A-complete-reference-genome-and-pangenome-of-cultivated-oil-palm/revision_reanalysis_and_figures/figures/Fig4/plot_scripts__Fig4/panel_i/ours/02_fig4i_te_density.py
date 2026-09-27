#!/usr/bin/env python3
"""Fig. 4i: redraw the four TE-density curves and 95% CI ribbons from the author's audited
9/23 run (RUN-FIG4I-39ASSEMBLY-20260923-001, TE_density_class_position_curve_genome39_material33.tsv)
inside the existing 4i axes. Axes, ticks, labels (Arial 6/6.2/6.5 pt), legend, dashed TSS/TES
lines, colours, line width (0.95 pt) and ribbon opacity (0.16) are the original ones.
Display recipe = author's build_fig4i_genome39.py::plot_final: Savitzky-Golay order 2,
windows 7/11/7 per region, lo=min(lo,y), hi=max(hi,y), x = bin + 0.5 on 0..140.
(The same recipe applied to the old Source Data reproduces the old curves to 0.005 pt.)
in: steps/step1_4d.pdf  out: steps/step2_4i.pdf  +  sd/Fig4i_curve.tsv
"""
import csv, re
import numpy as np
from scipy.signal import savgol_filter
import fitz
from common import *

SRC, DST = FIX / "steps/step1_4d.pdf", FIX / "steps/step2_4i.pdf"
import sys as _s
if len(_s.argv) > 2: SRC, DST = Path(_s.argv[1]), Path(_s.argv[2])
TSV = FIX / "src/TE_density_class_position_curve_genome39_material33.tsv"
AX_X0, AX_W, AX_Y0, AX_H = 371.242, 139.56201, 283.6333, 84.47299      # form coordinates
K = (364.9493 - 283.6333) / 40                                          # pt per % (ticks 0 and 40)
CLASSES = {"Core": (".51 .78 .722", "Fm1293"), "Soft-core": (".902 .365 .427", "Fm1294"),
           "Shell": (".161 .686 .831", "Fm1295"), "Cloud": (".49 .804 .973", "Fm1296")}


def smooth(v):
    r = np.asarray(v, float).copy()
    for a, b, w in ((0, 20, 7), (20, 120, 11), (120, 140, 7)):
        r[a:b] = savgol_filter(r[a:b], window_length=w, polyorder=2, mode="interp")
    return r


rows = list(csv.DictReader(open(TSV), delimiter="\t"))
curves = {}
for c in CLASSES:
    d = sorted((r for r in rows if r["Pangenome_class"] == c), key=lambda r: int(r["Profile_bin"]))
    assert len(d) == 140 and [int(r["Profile_bin"]) for r in d] == list(range(140))
    y = smooth([float(r["Mean_material_TE_percent"]) for r in d])
    lo = np.minimum(np.clip(smooth([float(r["CI95_low"]) for r in d]), 0, 100), y)
    hi = np.maximum(np.clip(smooth([float(r["CI95_high"]) for r in d]), 0, 100), y)
    curves[c] = (d, y, lo, hi)

X = AX_X0 + (np.arange(140) + 0.5) * AX_W / 140
fy = lambda v: AX_Y0 + v * K
ymax = max(c[3].max() for c in curves.values())
assert fy(ymax) < AX_Y0 + AX_H, f"data exceed y-axis: {ymax}"

doc = fitz.open(SRC)
s = get_stream(doc)
res = doc.xref_object(MAIN_FORM)
for c, (col, fm) in CLASSES.items():
    d, y, lo, hi = curves[c]
    # ribbon: widen the Illustrator clip to the axes box, rewrite the form XObject path
    m = re.search(r"q ([\d.]+ [\d.]+ [\d.]+ -?[\d.]+) re W n (q 0 g 0 G/GS7 gs 0 G /%s Do Q Q)" % fm, s)
    s = s[:m.start()] + f"q {num(AX_X0)} {num(AX_Y0, 4)} {num(AX_W, 5)} {num(AX_H, 5)} re W n {m.group(2)}" + s[m.end():]
    fx = int(re.search(r"/%s (\d+) 0 R" % fm, res).group(1))
    old = get_stream(doc, fx)
    assert old.startswith(f"q 0 623.789 519 -621 re W n /Perceptual ri {col} rg/GS0 gs")
    poly = [(X[i], fy(hi[i])) for i in range(140)] + [(X[i], fy(lo[i])) for i in range(139, -1, -1)]
    p = " ".join(f"{num(px)} {num(py)} {'m' if k == 0 else 'l'}" for k, (px, py) in enumerate(poly))
    put_stream(doc, f"q 0 623.789 519 -621 re W n /Perceptual ri {col} rg/GS0 gs {p} h f Q", fx)
    # mean line
    m = re.search(r"q 1 0 0 1 371\.742 ([\d.]+) cm %s RG \.95 w (0 0 m .*? l) S Q" % re.escape(col), s)
    y0 = fy(y[0])
    body = " ".join(f"{num(X[i] - X[0])} {num(fy(y[i]) - y0)} {'m' if i == 0 else 'l'}" for i in range(140))
    s = s[:m.start()] + f"q 1 0 0 1 {num(X[0], 4)} {num(y0, 4)} cm {col} RG .95 w {body} S Q" + s[m.end():]
# pre-existing glitch: the "+2 kb" tick was drawn 3 pt too low (280.6314 instead of 283.6333),
# so it showed as a stray mark after "kb"; put it on the axis like the other three x ticks
s = replace_once(s, "q 1 0 0 1 510.805 280.6314 cm 0 0 m 0 -2 l S Q", "q 1 0 0 1 510.805 283.6333 cm 0 0 m 0 -2 l S Q")
put_stream(doc, s)
doc.save(DST, garbage=0, deflate=True)

with open(FIX / "sd/Fig4i_curve.tsv", "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["Pangenome_class", "Profile_bin", "Region", "Region_position", "Material_count",
                "Mean_material_TE_percent", "SD_material_TE_percent", "SEM_material_TE_percent",
                "CI95_low", "CI95_high", "Display_mean_SG", "Display_CI_low_SG", "Display_CI_high_SG"])
    for c, (d, y, lo, hi) in curves.items():
        for i, r in enumerate(d):
            w.writerow([r[k] for k in ("Pangenome_class", "Profile_bin", "Region", "Region_position", "Material_count",
                                       "Mean_material_TE_percent", "SD_material_TE_percent", "SEM_material_TE_percent",
                                       "CI95_low", "CI95_high")] + [f"{y[i]:.4f}", f"{lo[i]:.4f}", f"{hi[i]:.4f}"])
for c, (d, y, lo, hi) in curves.items():
    v = np.array([float(r["Mean_material_TE_percent"]) for r in d])
    print(c, "upstream/body/downstream mean %:", v[:20].mean().round(2), v[20:120].mean().round(2), v[120:].mean().round(2))
print("max displayed CI high %:", round(float(ymax), 2), "axis top %:", round(AX_H / K, 2))
