#!/usr/bin/env python3
"""600 dpi pixel diff + word diff between two consecutive steps; prints panels that changed."""
import sys, collections
import numpy as np, fitz
PAN = {"a": (0, 0, 172, 112), "b": (172, 0, 343, 116), "c": (343, 0, 519, 116), "d": (0, 112, 172, 228),
       "e": (172, 116, 343, 228), "f": (343, 116, 519, 233), "g": (0, 228, 172, 348), "h": (172, 228, 343, 348),
       "i": (343, 233, 519, 348), "j": (0, 348, 226, 621), "k": (226, 348, 519, 621)}
a, b = fitz.open(sys.argv[1])[0], fitz.open(sys.argv[2])[0]
pa, pb = a.get_pixmap(dpi=600, alpha=False), b.get_pixmap(dpi=600, alpha=False)
A = np.frombuffer(pa.samples, np.uint8).reshape(pa.height, pa.width, 3).astype(int)
B = np.frombuffer(pb.samples, np.uint8).reshape(pb.height, pb.width, 3).astype(int)
D = np.abs(A - B).max(axis=2) > 8
sc = 600 / 72
ch = {p: int(D[int(y0*sc):int(y1*sc), int(x0*sc):int(x1*sc)].sum()) for p, (x0, y0, x1, y1) in PAN.items()}
wa = collections.Counter(w[4] for w in a.get_text("words")); wb = collections.Counter(w[4] for w in b.get_text("words"))
print(sys.argv[2].split("/")[-1], "changed px:", {k: v for k, v in ch.items() if v}, "| -", dict(wa - wb), "| +", dict(wb - wa))
