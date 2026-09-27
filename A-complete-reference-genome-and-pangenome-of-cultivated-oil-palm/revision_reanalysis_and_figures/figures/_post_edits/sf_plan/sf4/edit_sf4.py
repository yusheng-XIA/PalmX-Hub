"""SF4 (sfplan S25): stage tick labels in d, '0d' -> '0 d' etc. (space between number and unit, as in all other figures)."""
import sys, collections, fitz
sys.path.insert(0, sys.path[0] + "/..")
import vecedit as V
src, dst = sys.argv[1], sys.argv[2]
d = fitz.open(src); p = d[0]; before = V.words(p)
labs = ["0d", "15d", "35d", "50d", "65d", "80d", "95d", "110d", "125d", "140d", "155d", "170d", "185d",
        "12h", "24h", "36h", "48h", "60h", "72h"]
n = 0
for t in labs:
    n += V.replace(p, t, t[:-1] + " " + t[-1], keep="end", expect=1)
assert n == 19
after = V.words(p)
exp = collections.Counter(before); exp.subtract(labs)
for t in labs: exp[t[:-1]] += 1; exp[t[-1]] += 1
got = collections.Counter(after); diff = {k: v for k, v in (got - +exp).items()} | {k: -v for k, v in (+exp - got).items()}
assert not diff, diff
d.save(dst, garbage=3, deflate=True); print("SF4 labels replaced:", n)

# S25 (second half): mark the Level 2 candidate annotations (Source Data SF4d_metabolites evidence_level) with a dagger
d = fitz.open(dst); p = d[0]
p.insert_font(fontname="AR", fontfile=V.FONTS["ArialMT"])
lvl2 = {"FA 18:1+3O", "Azelaic acid", "Hyperoside", "-Coumaric acid"}
k = 0
for s, dr in list(V.spans(p)):
    if s["text"] in lvl2 and 600 < s["bbox"][1] < 680 and s["bbox"][0] < 75:
        p.insert_text((s["bbox"][2] + 0.3, s["origin"][1] - 2.0), "†", fontname="AR", fontsize=4.5, color=(0.133,) * 3)
        k += 1
assert k == 4, k
d.save(dst + ".tmp.pdf", garbage=3, deflate=True)
import os; os.replace(dst + ".tmp.pdf", dst); print("daggers:", k)
