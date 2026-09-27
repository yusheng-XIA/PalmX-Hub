#!/usr/bin/env python3
"""Fig. 4c: the panel itself never says "weighted" -- its only F_ST wording is the colour-bar
title "F" + subscript "ST" (upright ArialMT, 6.5 / 4.55 pt). To carry the author's decision
("mean F_ST") into the panel, "Mean" is added as a first title line above F_ST
(Arial TT0 6.5 pt, black, left-aligned with F, baseline 7.4 pt above). F_ST glyphs, style
and all values are untouched.  in: steps/step3_4b.pdf  out: steps/step4_4c.pdf
"""
import fitz
from common import *

SRC, DST = FIX / "steps/step3_4b.pdf", FIX / "steps/step4_4c.pdf"
import sys as _s
if len(_s.argv) > 2: SRC, DST = Path(_s.argv[1]), Path(_s.argv[2])
doc = fitz.open(SRC)
s = get_stream(doc)
old = "6.5 0 0 6.5 490.9639 586.6758 Tm (F)Tj 4.55 0 0 4.55 494.9336 585.5645 Tm (ST)Tj"
new = "6.5 0 0 6.5 490.9639 594.0758 Tm (Mean)Tj " + old
s = replace_once(s, old, new)
put_stream(doc, s)
doc.save(DST, garbage=0, deflate=True)

# Source Data: same numbers, honest column name (values = N_VAR-weighted average of the
# per-window Weir & Cockerham MEAN_FST, vcftools 100-kb windows / 50-kb step; recompute_abc.log)
import csv, openpyxl
WEIGHTED = {("K4_Pop1", "K4_Pop2"): 0.153126, ("K4_Pop1", "K4_Pop3"): 0.076072, ("K4_Pop1", "K4_Pop4"): 0.049124,
            ("K4_Pop2", "K4_Pop3"): 0.111356, ("K4_Pop2", "K4_Pop4"): 0.158345, ("K4_Pop3", "K4_Pop4"): 0.088406}
wb = openpyxl.load_workbook(SD_XLSX, read_only=True)
rows = list(wb["Fig.4c_fst_long"].iter_rows(values_only=True))
assert rows[0] == ("pop1", "pop2", "FST_weighted", "n_used_weighted", "n_windows")
with open(FIX / "sd/Fig4c_fst.tsv", "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["pop1", "pop2", "FST_mean", "n_windows_used", "n_windows",
                "for_reference_WEIGHTED_FST_not_plotted"])
    for r in rows[1:]:
        w.writerow([r[0], r[1], r[2], r[3], r[4], WEIGHTED[(r[0], r[1])]])
print("mean F_ST range", min(r[2] for r in rows[1:]), max(r[2] for r in rows[1:]))
