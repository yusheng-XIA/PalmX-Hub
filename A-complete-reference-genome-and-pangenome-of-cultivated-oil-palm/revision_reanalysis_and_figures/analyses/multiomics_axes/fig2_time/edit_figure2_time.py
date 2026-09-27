#!/usr/bin/env python3
"""Figure 2: put every divergence time on the final MCMCTree run (user, job 847506, 2026-06-17; PAML 4.10.9,
seven calibrations). Vector edits with PyMuPDF on the current Figure 2 (deliver/Main_Figures_revised/Figure2.pdf).

2a  node labels  ~48.43 / ~42.37 / ~33.85 Ma (old run, 2026-02-23)  ->  ~40.40 / ~33.69 / ~22.97 Ma (final run).
    ~4.92 and ~1.37 Ma are already final-run values and are left untouched.
    Same font (Arial Bold), size (5.93 pt), colour and baseline; digits have equal advance width in Arial, so the
    left edge is unchanged.
2b  (--2b-geometry) the time-scaled tree was drawn from the old run (node x = X0 - K * age, fitted below).
    Every internal node is moved to its final-run posterior-mean age: branch paths are redrawn with each point
    shifted by its parent node's shift (vertical stem and rounded corner) and the last point by the child node's
    shift (0 for tips); node circles move with their node; the +/- CAFE5 labels of each internal node move with it.
    Axis, highlight box, legend and tip labels are redrawn/kept unchanged. CAFE5 numbers are not changed here.

--node-circles TSV (with --2b-geometry): re-size and re-colour the 11 node circles by the rule of the authors' plotting
    script (phylo_divtime_cafe/plot_palm_v6.R): net = expanded - contracted; ggplot size = min(6, max(2.5,
    2.2 * log10(|net| + 1))); fill/border red if net > 0, blue if net < 0, grey if 0. Size -> radius in this PDF was
    calibrated on the circles whose sizes the current figure already has right (r = 0.2855 * size + 0.122 pt,
    fitted on <2>, <4>, <6>, <8>, <10>, <14>, <19>; residual <= 0.03 pt). TSV = Source Data Fig2b_gene_families layout.

usage: edit_figure2_time.py SRC.pdf OUT.pdf [--2b-geometry] [--node-circles TSV] [--log LOG.tsv]
"""
import sys
from pathlib import Path

import fitz

SRC, OUT = sys.argv[1], sys.argv[2]
GEOM = "--2b-geometry" in sys.argv
LOG = sys.argv[sys.argv.index("--log") + 1] if "--log" in sys.argv else None
CIRC = sys.argv[sys.argv.index("--node-circles") + 1] if "--node-circles" in sys.argv else None
TN2CAFE = {"t_n43": 2, "t_n42": 4, "t_n41": 6, "t_n40": 8, "t_n39": 10, "t_n38": 14, "t_n37": 15, "t_n36": 18,
           "t_n35": 19, "t_n34": 20, "t_n33": 22}
STYLE = {1: ((0.733, 0.302, 0.318), (0.988, 0.914, 0.918)),    # expansion (colours of the current figure)
         -1: ((0.278, 0.451, 0.627), (0.886, 0.929, 0.969)),   # contraction
         0: ((0.498, 0.498, 0.498), (0.878, 0.878, 0.878))}    # no net change: grey50 / #E0E0E0 as in plot_palm_v6.R
NET = {}
if CIRC:
    import csv as _csv, math as _m
    for r in _csv.DictReader(open(CIRC), delimiter="\t"):
        NET[int(r["node_id"])] = int(r["expanded_families"]) - int(r["contracted_families"])


def circle_items(cx, cy, r):
    k = 0.5523 * r
    P = fitz.Point
    return [("c", P(cx + r, cy), P(cx + r, cy + k), P(cx + k, cy + r), P(cx, cy + r)),
            ("c", P(cx, cy + r), P(cx - k, cy + r), P(cx - r, cy + k), P(cx - r, cy)),
            ("c", P(cx - r, cy), P(cx - r, cy - k), P(cx - k, cy - r), P(cx, cy - r)),
            ("c", P(cx, cy - r), P(cx + k, cy - r), P(cx + r, cy - k), P(cx + r, cy))]
AR = "/System/Library/Fonts/Supplemental/Arial.ttf"
ARB = "/System/Library/Fonts/Supplemental/Arial Bold.ttf"

# posterior means (Ma). OLD = FigTree.tre node heights of the 2026-02-23 run (what the current 2b geometry encodes);
# NEW = out.txt posterior means of the 2026-06-17 run (= node_ages_new_run.tsv), t_n labels of the new run.
NODES = [  # name, t_n, old age, new age
    ("Oryza vs (Musa + palms) [root of 2b]", "t_n33", 111.38, 111.01),
    ("Musa vs palms", "t_n34", 97.07, 96.98),
    ("Musa acuminata vs M. balbisiana", "t_n35", 19.54, 19.12),
    ("Palm crown", "t_n36", 72.92, 72.74),
    ("Calamus vs Daemonorops", "t_n37", 19.43, 18.57),
    ("Nypa vs other palms", "t_n38", 52.64, 44.78),
    ("Phoenix vs (Areca + Cocos + Elaeis)", "t_n39", 48.43, 40.40),
    ("Areca vs (Cocos + Elaeis)", "t_n40", 42.37, 33.69),
    ("Cocos vs Elaeis", "t_n41", 33.85, 22.97),
    ("E. oleifera vs E. guineensis", "t_n42", 5.14, 4.92),
    ("Dura vs pisifera", "t_n43", 1.40, 1.37),
]
# node positions in the current PDF (centre of node circle / start of branch paths)
NODE_XY = {"t_n33": (5.94, 330.58), "t_n34": (24.97, 310.27), "t_n35": (128.09, 331.37), "t_n36": (57.09, 289.18),
           "t_n37": (128.24, 305.34), "t_n38": (84.06, 273.01), "t_n39": (89.66, 260.20), "t_n40": (97.72, 247.60),
           "t_n41": (109.05, 235.40), "t_n42": (147.24, 224.01), "t_n43": (152.21, 214.25)}
# CAFE5 +/- labels belonging to each internal node: (text, x0, y0) of the spans in the current PDF
NODE_LABELS = {"t_n34": [("+13", 10.6, 302.4), ("−237", 10.6, 311.6)],
               "t_n35": [("+5058", 111.0, 323.4), ("−751", 111.0, 332.7)],
               "t_n36": [("+56", 38.6, 281.4), ("−2302", 38.6, 289.2)],
               "t_n37": [("+1412", 110.8, 298.3), ("−1680", 110.8, 306.4)],
               "t_n38": [("+13", 68.3, 265.1), ("−1049", 68.1, 273.5)],
               "t_n39": [("+35", 79.3, 251.7), ("−178", 91.1, 261.3)],
               "t_n40": [("+66", 86.9, 239.6), ("−561", 98.3, 248.9)],
               "t_n41": [("+192", 97.0, 228.3), ("−169", 98.2, 236.4)],
               "t_n42": [("+478", 131.6, 215.9), ("−1780", 131.6, 224.1)],
               "t_n43": [("+270", 138.9, 207.1), ("−349", 152.3, 207.8)]}
TIP_X = 154.10
PANEL_2B = fitz.Rect(0, 199, 262, 361)

doc = fitz.open(SRC)
pg = doc[0]
fR, fB = fitz.Font(fontfile=AR), fitz.Font(fontfile=ARB)
spans = [s for b in pg.get_text("dict")["blocks"] for l in b.get("lines", []) for s in l["spans"]]


def find(text, near, tol=1.5):
    out = [s for s in spans if s["text"] == text and abs(s["bbox"][0] - near[0]) < tol and abs(s["bbox"][1] - near[1]) < tol]
    if len(out) != 1:
        raise SystemExit(f"span {text!r} near {near}: {len(out)} hits")
    return out[0]


text_edits, log = [], []


def rewrite(s, new, dx=0.0, why=""):
    bold = "Bold" in s["font"]
    ox, oy = s["origin"]
    yc = oy - 0.33 * s["size"]
    pg.add_redact_annot(fitz.Rect(s["bbox"][0] + 0.3, yc - 0.12, s["bbox"][2] - 0.3, yc + 0.12), fill=False)
    c = s["color"]
    rgb = ((c >> 16) & 255) / 255, ((c >> 8) & 255) / 255, (c & 255) / 255
    text_edits.append((ox + dx, oy, new, "ArialBT" if bold else "ArialRT", s["size"], rgb))
    log.append(("text", s["text"], new, f"{ox:.2f},{oy:.2f}", f"{ox + dx:.2f},{oy:.2f}", why))


# ---------------- 2a node labels
for old, new, near, tn in [("48.43 Ma", "40.40 Ma", (220.9, 7.7), "t_n39"), ("42.37 Ma", "33.69 Ma", (286.9, 13.8), "t_n40"),
                           ("33.85 Ma", "22.97 Ma", (356.7, 30.5), "t_n41")]:
    rewrite(find(old, near), new, why=f"2a {tn}: old run -> final run posterior mean")
for keep, near in [("4.92 Ma", (405.9, 40.1)), ("1.37 Ma", (436.6, 58.5))]:
    find(keep, near)  # asserts presence; already final-run values

# ---------------- 2b geometry
if GEOM:
    # fit x = X0 - K * age on the current drawing (old ages)
    xs = [NODE_XY[t][0] for _, t, _, _ in NODES]
    ages = [a for _, _, a, _ in NODES]
    n = len(xs); ma = sum(ages) / n; mx = sum(xs) / n
    K = -sum((a - ma) * (x - mx) for a, x in zip(ages, xs)) / sum((a - ma) ** 2 for a in ages)
    X0 = mx + K * ma
    resid = max(abs(x - (X0 - K * a)) for x, a in zip(xs, ages))
    log.append(("fit", f"X0={X0:.3f}", f"K={K:.4f} pt/Ma", f"max residual {resid:.3f} pt", "", "old ages vs node x"))
    if resid > 0.3:
        raise SystemExit(f"2b axis fit residual {resid:.3f} pt: current drawing is not the old tree")
    SHIFT = {t: K * (old - new) for _, t, old, new in NODES}   # +x = younger
    for name, t, old, new in NODES:
        log.append(("node", t, name, f"{old:.2f} Ma @x={NODE_XY[t][0]:.2f}", f"{new:.2f} Ma @x={NODE_XY[t][0] + SHIFT[t]:.2f}",
                    f"dx={SHIFT[t]:+.2f} pt"))

    def node_at(x, y=None, tol=0.35):
        hits = [t for t, (nx, ny) in NODE_XY.items() if abs(nx - x) < tol and (y is None or abs(ny - y) < tol)]
        return hits[0] if len(hits) == 1 else None

    drawings = [g for g in pg.get_drawings() if PANEL_2B.contains(g["rect"])]
    if len(drawings) != 43:
        raise SystemExit(f"expected 43 vector objects in 2b, found {len(drawings)}")
    redrawn = []
    n_edge = n_circ = 0
    for g in drawings:
        pts_items = g["items"]
        dx_map = None
        if g["type"] == "s" and all(it[0] == "l" for it in pts_items) and g["lineCap"][0] == 1:
            # branch: first point = parent node; last point = child node or tip
            p0 = pts_items[0][1]; pl = pts_items[-1][2]
            par = node_at(p0.x, p0.y)
            chi = None if abs(pl.x - TIP_X) < 0.1 else node_at(pl.x, pl.y)
            if par is None or (chi is None and abs(pl.x - TIP_X) >= 0.1):
                raise SystemExit(f"unmatched branch {p0} -> {pl}")
            dp, dc = SHIFT[par], (SHIFT[chi] if chi else 0.0)
            new_items = []
            for k, it in enumerate(pts_items):
                a, b = it[1], it[2]
                last = k == len(pts_items) - 1
                new_items.append(("l", fitz.Point(a.x + dp, a.y), fitz.Point(b.x + (dc if last else dp), b.y)))
            # the stem + corner must stay left of the child end
            if max(it[2].x for it in new_items[:-1]) >= new_items[-1][2].x:
                raise SystemExit(f"branch {par}->{chi or 'tip'} too short for its corner after the shift")
            n_edge += 1
            log.append(("branch", par, chi or "tip", f"x {p0.x:.2f}->{pl.x:.2f}", f"x {p0.x + dp:.2f}->{pl.x + dc:.2f}", ""))
            g = dict(g, items=new_items)
        elif g["type"] == "fs" and all(it[0] == "c" for it in pts_items):
            r = g["rect"]; cx, cy = (r.x0 + r.x1) / 2, (r.y0 + r.y1) / 2
            t = node_at(cx, cy, tol=0.35)
            d = SHIFT[t] if t else 0.0
            new_items = [("c",) + tuple(fitz.Point(p.x + d, p.y) for p in it[1:]) for it in pts_items]
            if t:
                n_circ += 1
                note = ""
                if CIRC:
                    net = NET[TN2CAFE[t]]
                    size = min(6.0, max(2.5, 2.2 * _m.log10(abs(net) + 1)))
                    rad = 0.2855 * size + 0.122
                    sign = (net > 0) - (net < 0)
                    new_items = circle_items(cx + d, cy, rad)
                    g = dict(g, color=STYLE[sign][0], fill=STYLE[sign][1])
                    note = f"CAFE <{TN2CAFE[t]}> net {net:+d}: size {size:.2f}, r {(r.x1 - r.x0) / 2:.2f} -> {rad:.2f} pt, {['grey', 'red', 'blue'][sign]}"
                log.append(("circle", t, "", f"cx {cx:.2f}", f"cx {cx + d:.2f}", note))
            g = dict(g, items=new_items)
        redrawn.append(g)
    if n_edge != 22 or n_circ != 11:
        raise SystemExit(f"matched {n_edge} branches / {n_circ} node circles (expected 22 / 11)")
    # node labels
    for t, labs in NODE_LABELS.items():
        for text, x0, y0 in labs:
            rewrite(find(text, (x0, y0)), text, dx=SHIFT[t], why=f"2b label of {t} moved with its node")

# pass 1: text (thin redaction rectangles through glyph centres, text only)
pg.apply_redactions(images=fitz.PDF_REDACT_IMAGE_NONE, graphics=fitz.PDF_REDACT_LINE_ART_NONE,
                    text=fitz.PDF_REDACT_TEXT_REMOVE)
if GEOM:
    # pass 2: remove the 43 vector objects of 2b (all fully inside PANEL_2B; nothing from other panels is)
    pg.add_redact_annot(PANEL_2B, fill=False)
    pg.apply_redactions(images=fitz.PDF_REDACT_IMAGE_NONE, graphics=fitz.PDF_REDACT_LINE_ART_REMOVE_IF_COVERED,
                        text=fitz.PDF_REDACT_TEXT_NONE)
    left = [g for g in pg.get_drawings() if PANEL_2B.contains(g["rect"])]
    if left:
        raise SystemExit(f"{len(left)} 2b vector objects not removed")
    # redraw underneath all page content, in the original order (box first, then branches, circles, axis)
    sh = pg.new_shape()
    for g in redrawn:
        for it in g["items"]:
            if it[0] == "l":
                sh.draw_line(it[1], it[2])
            elif it[0] == "c":
                sh.draw_bezier(it[1], it[2], it[3], it[4])
            elif it[0] == "re":
                sh.draw_rect(it[1])
        typ = g["type"]
        sh.finish(color=g.get("color") if "s" in typ else None, fill=g.get("fill") if "f" in typ else None,
                  width=g.get("width") or 0.5, lineCap=(g.get("lineCap") or (0,))[0], lineJoin=int(g.get("lineJoin") or 0),
                  closePath=bool(g.get("closePath")), even_odd=bool(g.get("even_odd")),
                  fill_opacity=g.get("fill_opacity") or 1, stroke_opacity=g.get("stroke_opacity") or 1)
    sh.commit(overlay=False)

# new resource names: the page already carries subsetted "ArialR"/"ArialB" from fix_fig2.py (digits missing)
for x0, y, new, fn, sz, rgb in text_edits:
    pg.insert_text((x0, y), new, fontname=fn, fontfile=(ARB if fn == "ArialBT" else AR), fontsize=sz, color=rgb)
doc.subset_fonts()
# Arial maps U+002D and U+00AD to the same glyph (GID 16); keep the ToUnicode entry a normal hyphen (as fix_fig2.py)
for fx, *_rest in pg.get_fonts(full=True):
    base = _rest[2]
    if base.startswith("Arial") and "+" not in base:
        tu = doc.xref_get_key(fx, "ToUnicode")
        if tu[0] == "xref":
            tx = int(tu[1].split()[0])
            cmap = doc.xref_stream(tx).decode("latin1").replace("<0010> <00ad>", "<0010> <002d>")
            doc.update_stream(tx, cmap.encode("latin1"))
doc.save(OUT, garbage=4, deflate=True)
if LOG:
    with open(LOG, "w") as f:
        f.write("kind\told\tnew_or_item\tbefore\tafter\tnote\n")
        for r in log:
            f.write("\t".join(map(str, r)) + "\n")
for r in log:
    print(*r, sep=" | ")
