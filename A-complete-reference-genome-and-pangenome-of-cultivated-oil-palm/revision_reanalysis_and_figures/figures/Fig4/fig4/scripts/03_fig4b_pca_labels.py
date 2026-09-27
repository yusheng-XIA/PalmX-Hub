#!/usr/bin/env python3
"""Fig. 4b: PC axis labels -> share of total genetic variance (eigenvalue / trace of the GRM).
plink 1.9 --make-rel square on all.LDfilter.vcf (trace/Fig4/local/ldrel.rel, 308 x 308):
trace = 314.434173; eigenvalues 15.4830 / 10.5136 (identical to the published eigenval file)
-> PC1 4.9241 %, PC2 3.3436 %. Old labels (23.99 % / 16.29 %) were shares of the top-10 sum.
Same font resource (Arial TT0), size 6.5 pt, colour and baseline; label re-centred on the old centre.
in: steps/step2_4i.pdf  out: steps/step3_4b.pdf  + sd/Fig4b_PCA_variance.tsv
"""
import csv
import numpy as np, fitz
from common import *

SRC, DST = FIX / "steps/step2_4i.pdf", FIX / "steps/step3_4b.pdf"
import sys as _s
if len(_s.argv) > 2: SRC, DST = Path(_s.argv[1]), Path(_s.argv[2])
A = np.loadtxt(SCR / "trace/Fig4/local/ldrel.rel")
w = np.linalg.eigvalsh(A)[::-1]
tr = float(np.trace(A))
pct = w / tr * 100
assert abs(w[0] - 15.483) < 1e-3 and abs(w[1] - 10.5136) < 1e-3
p1, p2 = f"{pct[0]:.2f}", f"{pct[1]:.2f}"
assert (p1, p2) == ("4.92", "3.34")

doc = fitz.open(SRC)
wd = tt0_widths(doc)
s = get_stream(doc)
for old_s, new_s, horiz, (x, y) in (("PC1 (23.99%)", f"PC1 ({p1}%)", True, (249.7295, 508.9395)),
                                    ("PC2 (16.29%)", f"PC2 ({p2}%)", False, (182.0273, 544.1074))):
    shift = (wd(old_s) - wd(new_s)) * 6.5 / 2
    esc = lambda t: t.replace("(", r"\(").replace(")", r"\)")
    if horiz:
        old = f"6.5 0 0 6.5 {x} {y} Tm ({esc(old_s)})Tj"
        new = f"6.5 0 0 6.5 {num(x + shift, 4)} {y} Tm ({esc(new_s)})Tj"
    else:
        old = f"0 6.5 -6.5 0 {x} {y} Tm ({esc(old_s)})Tj"
        new = f"0 6.5 -6.5 0 {x} {num(y + shift, 4)} Tm ({esc(new_s)})Tj"
    s = replace_once(s, old, new)
put_stream(doc, s)
doc.save(DST, garbage=0, deflate=True)
with open(FIX / "sd/Fig4b_PCA_variance.tsv", "w", newline="") as fh:
    wr = csv.writer(fh, delimiter="\t")
    wr.writerow(["PC", "Eigenvalue", "Pct_of_total_variance(GRM_trace=%.6f)" % tr, "Pct_of_top10_sum(old_label)"])
    for i in range(10):
        wr.writerow([f"PC{i+1}", f"{w[i]:.4f}", f"{pct[i]:.4f}", f"{w[i] / w[:10].sum() * 100:.4f}"])
print("trace", tr, "PC1", p1, "PC2", p2)
