"""Fig. 1b telomere column: update 7 counts to the final-assembly recount and re-bin the colour squares.

Input (read-only): deliver/Main_Figures_revised/Figure1.pdf
Output:            fix/fig1b_telo/Figure1_candidate.pdf

Vector edit in the panel's own Form XObject content stream (the one holding the telomere squares
'488.982 <y> 4.535 -4.535 re B'): for each target row, the fill operator in front of the square and the
first literal string of the label '(NN)Tj 0 Tc 0 Tw (/)Tj ... (3)Tj ... (2)Tj' are swapped in place.
Font (ArialMT subset TT1), size, colour, glyph positions and baseline are untouched (all digits 556 units wide).

Colour bins read from the original figure (fill rgb, PDF operands):
  32  -> .196 .549 .6      31 -> .42 .667 .694      30 -> .655 .804 .82      <=29 -> .831 .902 .91
  (checked: FL 32/32 darkest; 31/32 rows; 30/32 rows; NS-Hap1 29/32 and EO12 18/32 both use .831 .902 .91)
"""
import re, fitz
from pathlib import Path

HERE = Path(__file__).resolve().parent
S = HERE.parents[1]
SRC = S / "deliver/Main_Figures_revised/Figure1.pdf"
OUT = HERE / "Figure1_candidate.pdf"

BINS = {32: b".196 .549 .6", 31: b".42 .667 .694", 30: b".655 .804 .82"}
LOW = b".831 .902 .91"
def colour(n): return BINS.get(n, LOW if n <= 29 else None)

# stream-space y of each row's telomere square (top edge), old value, new value
EDITS = [  # row, square_y, old, new
    ("Nigerian-Hap1",    b"526.967", 30, 22),
    ("Nigerian-Hap2",    b"512.229", 31, 27),
    ("TK-Hap2",          b"482.748", 31, 29),
    ("NS-Hap1",          b"468.006", 29, 28),
    ("NS-Hap2",          b"453.268", 31, 30),
    ("TN-Hap2",          b"409.045", 30, 32),
    ("E. oleifera-Hap1", b"364.826", 30, 29),
]

doc = fitz.open(SRC)
page = doc[0]
# locate the Form XObject holding the telomere squares
cands = [x for (x, *_r) in page.get_xobjects()
         if b"4.5350039 -4.5350039 re B" in (doc.xref_stream(x) or b"")]
assert len(cands) == 1, cands
xref = cands[0]
s = doc.xref_stream(xref)

log = []
for row, y, old, new in EDITS:
    pat = re.compile(rb"(?P<rgb>[\d.]+ [\d.]+ [\d.]+) rg (?P<w>\.188 w )?488\.982 " + re.escape(y) +
                     rb" 4\.5350039 -4\.53\d+ re B (?P<mid>\.133 \.133 \.133 rg BT -\.055 Tc \.055 Tw 6 0 0 6 497\.4854 [\d.]+ Tm )"
                     rb"\((?P<num>\d\d)\)Tj 0 Tc 0 Tw \(/\)Tj 1\.334 0 Td \(3\)Tj \.501 0 Td \(2\)Tj")
    ms = list(pat.finditer(s))
    assert len(ms) == 1, (row, len(ms))
    m = ms[0]
    assert int(m["num"]) == old, (row, m["num"], old)
    assert m["rgb"] == colour(old), (row, m["rgb"], colour(old))   # original follows the binning
    rep = (colour(new) + b" rg " + (m["w"] or b"") + b"488.982 " + y +
           m.group(0)[m.group(0).index(b" 4.5350039"):m.start("num") - m.start()] +
           str(new).encode() + m.group(0)[m.end("num") - m.start():])
    s = s[:m.start()] + rep + s[m.end():]
    log.append((row, old, new, m["rgb"].decode(), colour(new).decode()))

doc.update_stream(xref, s)
doc.save(OUT, garbage=3, deflate=True)
for r in log:
    print("%-17s %d/32 -> %d/32   fill %s -> %s" % r)
print("xref", xref, "->", OUT)
