#!/usr/bin/env python3
"""Figure 5h/5i candidate: update favourable-locus capture to the current W = 0, P = 15 path.

Base (read-only): scratchpad/deliver/Main_Figures_revised/Figure5.pdf.
Vector-only stream edits in the two copies of the original page content (xref 63, clipped
panel-a copy; xref 8270, visible b-i copy):
  * triangle fills: '.851 .643 .255 scn' (#D9A441 captured) <-> '.722 .773 .8 scn' (#B8C5CC missed)
    for positions whose status changed (position_marks.tsv, one triangle per unique coordinate,
    as in 85_plot_african35_legacy_v5_gwas_vector.py: captured if any locus at that coordinate is);
  * 5h right-hand labels '* c/t': digit-for-digit substitution inside the existing Tj strings
    (Arial digits are equal-width, and every new value has the same number of digits), so font,
    size, kerning and positions are unchanged.
Triangles are matched to positions per row by rank order of x vs genomic coordinate; the script
asserts a linear x-Mb fit (max residual < 0.05 pt) and that every legacy colour equals the
legacy W = 4 mark before changing anything.
"""
import csv, re, sys
from collections import defaultdict
from pathlib import Path
import fitz
import numpy as np

HERE = Path(__file__).resolve().parent
BASE = HERE.parent.parent / "deliver/Main_Figures_revised/Figure5.pdf"
OUT = HERE / "Figure5_candidate.pdf"
GOLD, GREY = b".851 .643 .255", b".722 .773 .8"
ROW_Y = {247.28: "chr01B", 237.41: "chr02B", 227.53: "chr03B", 217.65: "chr04B", 207.78: "chr05B",
         197.9: "chr06B", 188.03: "chr07B", 178.15: "chr08B", 168.27: "chr09B", 158.4: "chr10B",
         148.52: "chr11B", 138.64: "chr12B", 128.77: "chr13B", 118.9: "chr14B", 109.02: "chr15B",
         99.14: "chr16B"}
ZOOM_Y = 253.96   # 5i marker row (chr01B)

marks = defaultdict(list)
for r in csv.DictReader(open(HERE / "position_marks.tsv"), delimiter="\t"):
    marks[r["Chrom"]].append((int(r["Pos"]), int(r["Mark_W0"]), int(r["Mark_Legacy_W4"])))
for c in marks: marks[c].sort()
per = {r["Chrom"]: r for r in csv.DictReader(open(HERE / "per_chrom_capture.tsv"), delimiter="\t")}
OLD_LABEL = {"chr01B": "121/209", "chr02B": "5/5", "chr03B": "3/4", "chr04B": "1/1", "chr05B": "25/27",
             "chr06B": "14/14", "chr07B": "2/2", "chr08B": "0/0", "chr09B": "2/3", "chr10B": "2/4",
             "chr11B": "1/3", "chr12B": "3/4", "chr13B": "2/3", "chr14B": "0/0", "chr15B": "3/4",
             "chr16B": "0/1"}
CHROMS = [f"chr{i:02d}B" for i in range(1, 17)]
NEW_LABEL = {c: f"{per[c]['Captured_W0']}/{per[c]['Fav_Total']}" for c in CHROMS}

TRI = re.compile(rb"q 1 0 0 1 ([\d.\-]+) ([\d.\-]+) cm (\.851 \.643 \.255|\.722 \.773 \.8) scn")


def recolor(s, log):
    rows = defaultdict(list)
    for m in TRI.finditer(s):
        x, y = float(m.group(1)), round(float(m.group(2)), 2)
        rows[y].append((x, m.start(3), m.group(3)))
    subs = []
    for y, chrom, panel in [(yy, c, "5h") for yy, c in ROW_Y.items()] + [(ZOOM_Y, "chr01B", "5i")]:
        tri = sorted(rows.get(y, []))
        pos = marks.get(chrom, [])
        assert len(tri) == len(pos), (panel, chrom, len(tri), len(pos))
        if not tri: continue
        if len(tri) >= 2:
            X = np.array([t[0] for t in tri]); P = np.array([p[0] for p in pos]) / 1e6
            a, b = np.polyfit(P, X, 1); res = np.abs(X - (a * P + b)).max()
            fit = (a, b, res)
        for (x, off, col), (p, new, old) in zip(tri, pos):
            assert col == (GOLD if old else GREY), (panel, chrom, p, col, old)   # legacy colour check
            want = GOLD if new else GREY
            if want != col:
                subs.append((off, len(col), want))
                log.append((panel, chrom, p, round(x, 3), "gold" if old else "grey", "gold" if new else "grey"))
    # global linear-fit check for the 5h rows (shared x axis)
    X, P = [], []
    for y, chrom in ROW_Y.items():
        for (x, *_), (p, *__) in zip(sorted(rows.get(y, [])), marks.get(chrom, [])):
            X.append(x); P.append(p / 1e6)
    a, b = np.polyfit(P, X, 1); res = np.abs(np.array(X) - (a * np.array(P) + b)).max()
    assert res < 0.05, ("5h x-fit residual", res)
    Xz = [t[0] for t in sorted(rows[ZOOM_Y])]; Pz = [p[0] / 1e6 for p in marks["chr01B"]]
    az, bz = np.polyfit(Pz, Xz, 1); rz = np.abs(np.array(Xz) - (az * np.array(Pz) + bz)).max()
    assert rz < 0.05, ("5i x-fit residual", rz)
    for off, n, want in sorted(subs, reverse=True):
        s = s[:off] + want + s[off + n:]
    return s, (a, b, res, az, bz, rz)


TOK = re.compile(rb"\(((?:\\.|[^\\)])*)\)Tj|([\-\d.]+) ([\-\d.]+) Td|ET|Tm")


def relabel(s, log):
    """Rewrite the 16 '* c/t' labels. The block starts at the Tm that draws '4 / 2091'."""
    start = s.find(b"(2091)Tj")
    assert start > 0 and s.count(b"(2091)Tj") == 1
    i = start; row = 0; out = bytearray(s[:start]); pos = start
    state = None; pieces = []
    for m in TOK.finditer(s, start):
        if m.group(0) == b"ET":
            end = m.start(); break
        if m.group(0) == b"Tm":
            if state == "star" and row < 16: _flush(row, pieces, log); row += 1
            state = None; continue
        if m.group(1) is not None:
            t = m.group(1)
            if t == b"*": state = "star"; pieces = []; continue
            if state == "star":
                if t.strip() == b"" and not pieces: continue      # the space after '*'
                pieces.append((m.start(1), m.end(1), t))
        else:
            if abs(float(m.group(3))) > 1:                          # line change
                if state == "star" and row < 16: _flush(row, pieces, log); row += 1
                state = None
    else:
        raise RuntimeError("no ET")
    assert row == 16, row
    s2 = bytearray(s)
    for (a, b, new) in sorted(REPL, reverse=True):
        s2[a:b] = new
    return bytes(s2)


REPL = []


def _flush(row, pieces, log):
    chrom = CHROMS[row]
    old = b"".join(p[2] for p in pieces).decode()
    assert old == OLD_LABEL[chrom], (chrom, old)
    new = NEW_LABEL[chrom]
    assert len(new) == len(old) and all((a == "/") == (b == "/") for a, b in zip(old, new)), (old, new)
    k = 0
    for a, b, t in pieces:
        seg = new[k:k + len(t)].encode(); k += len(t)
        if seg != t: REPL.append((a, b, seg))
    log.append(("5h label", chrom, old, new))


def main():
    doc = fitz.open(BASE)
    tri_log, lab_log = [], []
    for xref in (63, 8270):
        s = doc.xref_stream(xref)
        tl = []
        s, fit = recolor(s, tl)
        REPL.clear(); ll = []
        s = relabel(s, ll)
        doc.update_stream(xref, s)
        print(f"xref {xref}: triangles changed {len(tl)}; labels {sum(1 for l in ll if l[2]!=l[3])} changed; "
              f"5h fit slope {fit[0]:.4f} pt/Mb res {fit[2]:.4f}; 5i res {fit[5]:.4f}")
        if xref == 8270: tri_log, lab_log = tl, ll
    doc.save(OUT, garbage=3, deflate=True)
    with open(HERE / "figure_edit_log.tsv", "w") as fh:
        fh.write("panel\tchrom\tpos_or_old\tstream_x_or_new\told\tnew\n")
        for r in tri_log: fh.write("\t".join(map(str, r)) + "\n")
        for r in lab_log: fh.write("\t".join(map(str, r)) + "\t\t\n")
    print("saved", OUT)


if __name__ == "__main__":
    main()
