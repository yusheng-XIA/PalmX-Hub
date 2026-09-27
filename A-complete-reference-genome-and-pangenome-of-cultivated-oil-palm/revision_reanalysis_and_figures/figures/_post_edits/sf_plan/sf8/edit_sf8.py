"""SF8 (sfplan S20): the right k panel title was a clipped fragment ('of best-matching cluster'); replace it with a complete
short title ('Best-matching cluster'); the y label ('FL share of 185-d nuclei') is unchanged."""
import sys, fitz
sys.path.insert(0, sys.path[0] + "/..")
import vecedit as V
src, dst = sys.argv[1], sys.argv[2]
d = fitz.open(src); p = d[0]
n = V.replace(p, "of best-matching cluster", "Best-matching cluster", keep="start", expect=1)
# centre the new title over the right k axes (same x-span as the old fragment)
d.save(dst, garbage=3, deflate=True); print("replaced", n)
