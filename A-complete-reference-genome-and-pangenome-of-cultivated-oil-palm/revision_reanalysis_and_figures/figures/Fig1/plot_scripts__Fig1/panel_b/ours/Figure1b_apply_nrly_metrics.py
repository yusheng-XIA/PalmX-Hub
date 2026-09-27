#!/usr/bin/env python3
"""Write new Nigerian-Hap1 / Nigerian-Hap2 metric values into Fig. 1b (vector edit, PyMuPDF).

Usage
  python3 apply_nrly_metrics.py VALUES.tsv [--pdf IN.pdf] [--out OUT.pdf] [--redraw-all] [--locate]

VALUES.tsv: tab-separated, header 'metric<TAB>Hap1<TAB>Hap2'. metric (case-insensitive) one of
  assembly size | scaffold N50 | contig N50 | LAI | QV | BUSCO | gaps | telomere
  Values are formatted with the column's own format (size 2 dp, N50 1 dp, LAI 1 dp, QV 2 dp, BUSCO 1 dp,
  gaps integer with thousands comma, telomere 'n/32'); 'NA' or empty = leave that cell unchanged.
  'contig N50' is accepted but Fig. 1b has no contig-N50 column -> reported and skipped.

Default --pdf is fix/fig1b_telo/Figure1_candidate.pdf (telomere-corrected), falling back to the
delivered Figure1.pdf. Default --out is fix/fig1b_telo/Figure1_candidate_nrly.pdf. The input is never modified.

How each Fig. 1b column is encoded (all inside one Form XObject of panel b, stream space, y up):
  assembly size  blue bar  '.306 .475 .655 rg X0 Y W -5.668 re f', X0=196.728, W = k*value (0-2 Gb axis);
                 label left-aligned at X0+W+1.843
  scaffold N50   orange bar '.776 .357 .157 rg', X0=243.216, W = k*value (0-150 Mb axis); label at bar end+1.843
  LAI            green bar  '.341 .616 .286 rg', X0=289.704, W = k*value; label at bar end+1.843
  QV             lollipop: stem (own Form XObject, '336.1929 y cm 0 0 m L 0 l S', L = k*QV) + dot
                 (d=3.684 circle, '.71 .604 .125 rg', centre x = 336.1929 + L); label at fixed x
  BUSCO          pie: pale disc (.945 .898 .937) + sector polygon (.635 .365 .624, 80 arc vertices,
                 r=3.684, from 12 o'clock counter-clockwise, sweep = BUSCO% x 360 deg); label at fixed x
  gaps           circle (.682 .549 .482), centre x fixed 444.759, diameter d with d^2 = a + b*log10(n+1)
                 (area ~ log gaps); label at fixed x
  telomere       square fill binned: 32 / 31 / 30 / <=29 (4 colours); label 'n/32' at fixed x
  Scales (k, a, b) are re-fitted from all 14 rows of the current figure against the values printed in it,
  and the fit residuals are reported.

Unchanged values (new text == printed text) are left byte-identical unless --redraw-all is given, in which
case every Nigerian element is regenerated from the fitted scale (used to validate the geometry model).
"""
import argparse, math, re, sys
from pathlib import Path
import fitz

HERE = Path(__file__).resolve().parent
S = HERE.parents[1]
NUM = rb"-?\d*\.?\d+"

# ---------------------------------------------------------------- column definitions
COLS = {
    "size":  dict(kind="bar", rgb=b".306 .475 .655", x0=196.728, fmt="{:.2f}", lab=(215, 240)),
    "n50":   dict(kind="bar", rgb=b".776 .357 .157", x0=243.216, fmt="{:.1f}", lab=(262, 285)),
    "lai":   dict(kind="bar", rgb=b".341 .616 .286", x0=289.704, fmt="{:.1f}", lab=(305, 335)),
    "qv":    dict(kind="qv", fmt="{:.2f}", lab=(367.0, 368.5)),
    "busco": dict(kind="pie", fmt="{:.1f}", lab=(404.0, 405.0)),
    "gaps":  dict(kind="gap", fmt="{:,.0f}", lab=(450.5, 451.5)),
    "telo":  dict(kind="telo", fmt="{:.0f}/32", lab=(497.0, 498.0)),
}
ALIASES = {"assembly size": "size", "size": "size", "scaffold n50": "n50", "n50": "n50",
           "contig n50": "contig", "lai": "lai", "qv": "qv", "busco": "busco", "gaps": "gaps",
           "telomere": "telo", "telomeres": "telo"}
TELO_BINS = {32: b".196 .549 .6", 31: b".42 .667 .694", 30: b".655 .804 .82"}
TELO_LOW = b".831 .902 .91"
WIDTH = {c: 556 for c in "0123456789"}; WIDTH.update({".": 278, ",": 278, "/": 278, "-": 333})


def f3(v):
    s = "%.3f" % v
    s = s.rstrip("0").rstrip(".") if "." in s else s
    if s.startswith("0."): s = s[1:]
    if s.startswith("-0."): s = "-" + s[2:]
    return s if s not in ("", "-") else "0"


# ---------------------------------------------------------------- text segments
TEXT_SEG = re.compile(rb"(?:(?P<tc>" + NUM + rb") Tc (?P<tw>" + NUM + rb") Tw )?(?P<tm>6 0 0 6) (?P<x>" + NUM + rb") (?P<y>" + NUM + rb") Tm (?P<body>(?:(?!Tm|ET).)*?)(?= ET|(?: " + NUM + rb"){6} Tm)", re.S)


def parse_segment(body):
    """glyph list [(char, x_offset_in_em)] from a Tj/Td/Tc run starting at a Tm (Tc/Tw state at start given)."""
    return re.findall(rb"\(([^)]*)\)Tj", body)


def glyph_text(body):
    return b"".join(parse_segment(body)).decode("latin1")


def substitute(body, new):
    """same-shape replacement: same length and non-digit characters in the same places."""
    old = glyph_text(body)
    if len(old) != len(new) or any((a.isdigit() != b.isdigit()) or (not a.isdigit() and a != b) for a, b in zip(old, new)):
        return None
    it = iter(new)
    return re.sub(rb"\(([^)]*)\)Tj", lambda m: b"(" + "".join(next(it) for _ in m.group(1).decode("latin1")).encode("latin1") + b")Tj", body)


def rebuild(new, adv):
    """new-shape text: one glyph per Tj, absolute offsets via Td, advance per char class taken from the
    column's existing labels (adv dict), falling back to font width - 0.011 em."""
    out, pos = [b"0 Tc 0 Tw"], 0.0
    for i, c in enumerate(new):
        if i:
            step = adv.get((new[i - 1], c), adv.get(new[i - 1], WIDTH[new[i - 1]] / 1000 - 0.011))
            out.append(("%s 0 Td" % f3(step)).encode()); pos += step
        out.append(b"(" + c.encode() + b")Tj")
    return b" ".join(out)


def learn_advances(cells_):
    """per-character advance (em) seen in existing labels: position of glyph i+1 minus glyph i.
    Tracks Tc for multi-glyph Tj runs and Td for explicit moves."""
    adv = {}
    for seg in cells_:
        toks = re.findall(rb"(" + NUM + rb") Tc|(" + NUM + rb") 0 Td|\(([^)]*)\)Tj", seg["body"])
        tc, line, cur, glyphs = seg["tc0"], 0.0, 0.0, []
        for tcv, td, tj in toks:
            if tcv: tc = float(tcv)
            elif td: line += float(td); cur = line
            else:
                for ch in tj.decode("latin1"):
                    glyphs.append((ch, cur)); cur += WIDTH[ch] / 1000 + tc
        for (a, xa), (b, xb) in zip(glyphs, glyphs[1:]):
            adv.setdefault(a, round(xb - xa, 3))
    return adv


# ---------------------------------------------------------------- locate everything
def locate(doc):
    page = doc[0]
    xs = [x for (x, *_r) in page.get_xobjects() if b"cm .635 .365 .624 rg" in (doc.xref_stream(x) or b"")]
    if len(xs) != 1: sys.exit("panel-b XObject not unique: %s" % xs)
    xref = xs[0]; s = doc.xref_stream(xref)

    # row centres from the BUSCO pies (one per row)
    pies = [(m.start(), m.end(), float(m["y"])) for m in re.finditer(
        rb"q 1 0 0 1 398\.2709 (?P<y>" + NUM + rb") cm \.635 \.365 \.624 rg .*? h f Q", s, re.S)]
    rows_y = sorted({round(p[2], 2) for p in pies}, reverse=True)
    assert len(rows_y) == 14, rows_y

    # map stream y -> page y with the BUSCO labels, then name rows with the haplotype labels on the page
    segs = [dict(x=float(m["x"]), y=float(m["y"]), xs=m["x"], ys=m["y"], a=m.start("tm"), be=m.end("body"),
                 tc0=float(m["tc"] or 0), tw0=float(m["tw"] or 0), body=m["body"]) for m in TEXT_SEG.finditer(s)]
    bl = sorted([g for g in segs if 367.0 < g["x"] < 368.5], key=lambda g: -g["y"])   # QV labels
    spans = [sp for b in page.get_text("dict")["blocks"] for l in b.get("lines", []) for sp in l["spans"]]
    pb = sorted([sp for sp in spans if abs(sp["origin"][0] - 341.27) < .1 and re.fullmatch(r"\d\d\.\d\d", sp["text"].strip())],
                key=lambda sp: sp["origin"][1])
    assert len(bl) == len(pb) == 14, (len(bl), len(pb))
    n = 14; sx = sum(g["y"] for g in bl) / n; sy = sum(sp["origin"][1] for sp in pb) / n
    b_ = sum((g["y"] - sx) * (sp["origin"][1] - sy) for g, sp in zip(bl, pb)) / sum((g["y"] - sx) ** 2 for g in bl)
    a_ = sy - b_ * sx
    resid = max(abs(a_ + b_ * g["y"] - sp["origin"][1]) for g, sp in zip(bl, pb))
    assert resid < .05, resid
    labels = [sp for sp in spans if abs(sp["origin"][0] - 60.78) < .1]
    names = {}
    for ry in rows_y:
        py = a_ + b_ * (ry - 1.98)   # row centre -> label baseline
        cand = [sp for sp in labels if abs(sp["origin"][1] - py) < 2.5]
        names[ry] = "".join(sp["text"] for sp in sorted(cand, key=lambda sp: sp["origin"][0])).replace("\xa0", " ").strip()

    def row_of(y):
        r = min(rows_y, key=lambda r: abs(r - y)); assert abs(r - y) < 3, y; return r

    cells = {}  # (row_y, col) -> dict(elements)
    def put(ry, col, key, val):
        cells.setdefault((ry, col), {})[key] = val

    for g in segs:
        if min(abs(r - g["y"] - 1.98) for r in rows_y) > 1.0: continue      # axis ticks, headers
        for col, c in COLS.items():
            if c["lab"][0] < g["x"] < c["lab"][1]:
                put(row_of(g["y"] + 1.98), col, "text", g)
    for col in ("size", "n50", "lai"):
        c = COLS[col]
        for m in re.finditer(re.escape(c["rgb"]) + rb" rg (?:" + NUM + rb" w )?(?P<x>" + NUM + rb") (?P<y>" + NUM +
                             rb") (?P<w>" + NUM + rb") (?P<h>-" + NUM + rb") re f", s):
            if abs(float(m["x"]) - c["x0"]) < .01:
                put(row_of(float(m["y"]) + float(m["h"]) / 2), col, "bar",
                    dict(a=m.start("w"), b=m.end("w"), w=float(m["w"])))
    for m in re.finditer(rb"q 1 0 0 1 (?P<x>" + NUM + rb") (?P<y>" + NUM + rb") cm \.71 \.604 \.125 rg 0 0 m 0 -(?P<k>" + NUM +
                         rb") -(?P<d>" + NUM + rb") ", s):
        put(row_of(float(m["y"])), "qv", "dot", dict(a=m.start("x"), b=m.end("x"), x=float(m["x"]), d=float(m["d"])))
    # QV stems: '/FmN Do' right before each dot; stem XObject holds '... cm 0 0 m L 0 l S'
    res = doc.xref_object(xref)
    for m in re.finditer(rb"/(?P<fm>Fm\d+) Do Q q 1 0 0 1 " + NUM + rb" (?P<y>" + NUM + rb") cm \.71 \.604 \.125 rg", s):
        fx = int(re.search(r"/%s (\d+) 0 R" % m["fm"].decode(), res).group(1))
        fs = doc.xref_stream(fx)
        mm = re.search(rb"q 1 0 0 1 (?P<x0>" + NUM + rb") (?P<y>" + NUM + rb") cm 0 0 m (?P<L>" + NUM + rb") 0 l S", fs)
        put(row_of(float(m["y"])), "qv", "stem", dict(xref=fx, a=mm.start("L"), b=mm.end("L"), L=float(mm["L"]), x0=float(mm["x0"]), bbox=doc.xref_get_key(fx, "BBox")))
    for m in re.finditer(rb"q 1 0 0 1 398\.2709 (?P<y>" + NUM + rb") cm \.635 \.365 \.624 rg (?P<path>.*? h) f Q", s, re.S):
        pts = [(float(a), float(b)) for a, b in re.findall(rb"(" + NUM + rb") (" + NUM + rb") [ml]", m["path"])]
        put(row_of(float(m["y"])), "busco", "pie", dict(a=m.start("path"), b=m.end("path"), pts=pts))
    for m in re.finditer(rb"q 1 0 0 1 (?P<x>" + NUM + rb") (?P<y>" + NUM + rb") cm \.682 \.549 \.482 rg (?P<path>0 0 m 0 -" + NUM +
                         rb" -(?P<d>" + NUM + rb") .*? c -" + NUM + rb" " + NUM + rb" 0 " + NUM + rb" 0 0 c)", s):
        put(row_of(float(m["y"])), "gaps", "circle", dict(a=m.start("x"), b=m.end("path"), x=float(m["x"]), y=m["y"], d=float(m["d"])))
    for m in re.finditer(rb"(?P<rgb>" + NUM + rb" " + NUM + rb" " + NUM + rb") rg (?:" + NUM + rb" w )?488\.982 (?P<y>" + NUM +
                         rb") 4\.5350039 (?P<h>-" + NUM + rb") re B", s):
        put(row_of(float(m["y"]) + float(m["h"]) / 2), "telo", "square", dict(a=m.start("rgb"), b=m.end("rgb"), rgb=m["rgb"]))
    # what the reader sees: page-level text per cell (some labels were re-set outside the panel XObject)
    PX = dict(size=(200, 216), n50=(245, 260), lai=(285, 305), qv=(340, 342), busco=(375, 376), gaps=(418, 419), telo=(461, 462))
    for ry in rows_y:
        py = a_ + b_ * (ry - 1.98)
        for col, (x0, x1) in PX.items():
            hit = [sp for sp in spans if x0 <= sp["origin"][0] <= x1 and abs(sp["origin"][1] - py) < 1.5]
            assert len(hit) == 1, (names[ry], col, [h["text"] for h in hit])
            put(ry, col, "page_text", hit[0]["text"].strip())
    return xref, s, rows_y, names, cells


# ---------------------------------------------------------------- scales
def printed_value(col, cell):
    t = cell["page_text"].replace(",", "")
    return float(t.split("/")[0]) if col == "telo" else float(t)


def fit_scales(rows_y, cells):
    sc, rep = {}, []
    for col in ("size", "n50", "lai"):
        pairs = [(printed_value(col, cells[(r, col)]), cells[(r, col)]["bar"]["w"]) for r in rows_y]
        k = sum(v * w for v, w in pairs) / sum(v * v for v, _ in pairs)
        off = [cells[(r, col)]["text"]["x"] - (COLS[col]["x0"] + cells[(r, col)]["bar"]["w"]) for r in rows_y]
        sc[col] = dict(k=k, lab_off=sum(off) / len(off))
        rep.append("%-5s bar width = %.5f x value; max |resid| %.3f pt; label offset %.3f (spread %.4f)" % (
            col, k, max(abs(w - k * v) for v, w in pairs), sc[col]["lab_off"], max(off) - min(off)))
    stems = [(printed_value("qv", cells[(r, "qv")]), cells[(r, "qv")]["stem"]) for r in rows_y if "stem" in cells[(r, "qv")]]
    k = sum(v * st["L"] for v, st in stems) / sum(v * v for v, _ in stems)
    sc["qv"] = dict(k=k, x0=stems[0][1]["x0"])
    rep.append("qv    stem length = %.5f x QV from x0=%.4f; max |resid| %.3f pt (%d stems; rows without stem: %d)" % (
        k, stems[0][1]["x0"], max(abs(st["L"] - k * v) for v, st in stems), len(stems), 14 - len(stems)))
    ds = [(math.log10(printed_value("gaps", cells[(r, "gaps")]) + 1), cells[(r, "gaps")]["circle"]["d"] ** 2) for r in rows_y]
    n = len(ds); mx = sum(x for x, _ in ds) / n; my = sum(y for _, y in ds) / n
    b = sum((x - mx) * (y - my) for x, y in ds) / sum((x - mx) ** 2 for x, _ in ds); a = my - b * mx
    cx = [cells[(r, "gaps")]["circle"]["x"] - cells[(r, "gaps")]["circle"]["d"] / 2 for r in rows_y]
    sc["gaps"] = dict(a=a, b=b, cx=sum(cx) / n)
    rep.append("gaps  diameter^2 = %.4f + %.4f x log10(n+1); max |d resid| %.3f pt; centre x %.3f (spread %.4f)" % (
        a, b, max(abs(math.sqrt(a + b * x) - math.sqrt(y)) for x, y in ds), sc["gaps"]["cx"], max(cx) - min(cx)))
    sw = []
    for r in rows_y:
        pts = cells[(r, "busco")]["pie"]["pts"]; ex, ey = pts[-1]
        ang = (2 * math.pi - math.atan2(ex, ey)) % (2 * math.pi)
        sw.append(abs(math.degrees(ang) / 360 * 100 - printed_value("busco", cells[(r, "busco")])))
    rep.append("busco sector sweep = BUSCO%% x 360 deg (CCW from 12 o'clock); max |%%-resid| %.3f (label rounding 0.05)" % max(sw))
    return sc, rep


def pie_path(pct, r=3.684, n=80):
    sweep = 2 * math.pi * pct / 100
    pts = [(0.0, 0.0)] + [(-r * math.sin(sweep * i / (n - 1)), r * math.cos(sweep * i / (n - 1))) for i in range(n)]
    return (" ".join("%s %s %s" % (f3(x), f3(y), "m" if i == 0 else "l") for i, (x, y) in enumerate(pts)) + " h").encode()


def circle_path(d):
    k = 2 * d / 3
    return ("0 0 m 0 -{k} -{d} -{k} -{d} 0 c -{d} {k} 0 {k} 0 0 c".format(k=f3(k), d=f3(d))).encode()


# ---------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("tsv")
    ap.add_argument("--pdf", default=None)
    ap.add_argument("--out", default=str(HERE / "Figure1_candidate_nrly.pdf"))
    ap.add_argument("--redraw-all", action="store_true")
    ap.add_argument("--locate", action="store_true", help="only print located cells and exit")
    A = ap.parse_args()
    src = Path(A.pdf) if A.pdf else (HERE / "Figure1_candidate.pdf" if (HERE / "Figure1_candidate.pdf").exists()
                                     else S / "deliver/Main_Figures_revised/Figure1.pdf")
    doc = fitz.open(src)
    xref, s, rows_y, names, cells = locate(doc)
    print("input:", src, "| panel-b XObject xref", xref)
    nig = {"Hap1": None, "Hap2": None}
    for r in rows_y:
        for h in nig:
            if names[r] == "Nigerian-" + h: nig[h] = r
    assert all(nig.values()), names
    for h, r in nig.items():
        got = {c: sorted(cells.get((r, c), {})) for c in COLS}
        print("Nigerian-%s row (stream y %.2f): " % (h, r) + "; ".join(
            "%s=%s[%s]" % (c, glyph_text(cells[(r, c)]["text"]["body"]), ",".join(k for k in got[c] if k != "text")) for c in COLS))
    sc, rep = fit_scales(rows_y, cells)
    for line in rep: print("  scale:", line)
    if A.locate: return

    vals = {}
    for i, line in enumerate(Path(A.tsv).read_text().splitlines()):
        if not line.strip() or i == 0 and line.lower().startswith("metric"): continue
        m, *v = [x.strip() for x in line.split("\t")]
        key = ALIASES.get(m.lower())
        if key is None: sys.exit("unknown metric: %r" % m)
        if key == "contig":
            print("  NOTE: contig N50 %s is not shown in Fig. 1b (no column) -> skipped" % v); continue
        vals[key] = dict(zip(("Hap1", "Hap2"), v + [""] * (2 - len(v))))

    edits = []                   # (start, end, bytes) in panel stream
    stem_edits = []              # (xref, start, end, bytes)
    changes = []
    for h, r in nig.items():
        for col, hv in vals.items():
            raw = hv[h]
            if raw in ("", "NA", "na", "-"): continue
            v = float(raw.replace(",", "").split("/")[0])
            new = COLS[col]["fmt"].format(v)
            cell = cells[(r, col)]; old = glyph_text(cell["text"]["body"])
            if new == old and not A.redraw_all:
                changes.append((h, col, old, new, "unchanged")); continue
            # --- graphics
            lab_x = None
            if COLS[col]["kind"] == "bar":
                if col == "size" and v > 2.0 or col == "n50" and v > 150 or col == "lai" and v > 30:
                    print("  WARNING %s %s=%s exceeds the drawn axis range" % (h, col, new))
                w = sc[col]["k"] * v
                edits.append((cell["bar"]["a"], cell["bar"]["b"], f3(w).encode()))
                lab_x = COLS[col]["x0"] + w + sc[col]["lab_off"]
            elif col == "qv":
                if v > 80: print("  WARNING %s QV %s exceeds the 0-80 axis" % (h, new))
                L = sc["qv"]["k"] * v; d = cell["dot"]["d"]
                edits.append((cell["dot"]["a"], cell["dot"]["b"], f3(sc["qv"]["x0"] + L + d / 2).encode()))
                st = cell["stem"]; stem_edits.append((st["xref"], st["a"], st["b"], f3(L).encode()))
            elif col == "busco":
                edits.append((cell["pie"]["a"], cell["pie"]["b"], pie_path(v)))
            elif col == "gaps":
                d = math.sqrt(sc["gaps"]["a"] + sc["gaps"]["b"] * math.log10(v + 1))
                c = cell["circle"]
                edits.append((c["a"], c["b"], (f3(sc["gaps"]["cx"] + d / 2) + " " + c["y"].decode() +
                                              " cm .682 .549 .482 rg ").encode() + circle_path(d)))
            elif col == "telo":
                n = int(v)
                edits.append((cell["square"]["a"], cell["square"]["b"], TELO_BINS.get(n, TELO_LOW if n <= 29 else b"0 0 0")))
            # --- label
            t = cell["text"]
            body = substitute(t["body"], new)
            if body is None:
                adv = learn_advances([cells[(rr, col)]["text"] for rr in rows_y])
                tc_end = re.findall(rb"(" + NUM + rb") Tc", t["body"]); tw_end = re.findall(rb"(" + NUM + rb") Tw", t["body"])
                body = rebuild(new, adv)
                # restore the text state the original segment left behind
                body += b" " + (tc_end[-1] if tc_end else f3(t["tc0"]).encode()) + b" Tc " + \
                    (tw_end[-1] if tw_end else f3(t["tw0"]).encode()) + b" Tw"
                print("  NOTE %s %s: '%s' -> '%s' changes shape; label rebuilt glyph-by-glyph" % (h, col, old, new))
            xs = f3(lab_x).encode() if lab_x is not None else t["xs"]
            edits.append((t["a"], t["be"], b"6 0 0 6 " + xs + b" " + t["ys"] + b" Tm " + body))
            changes.append((h, col, old, new, "edited"))

    edits.sort()
    for (a1, b1, _), (a2, b2, _) in zip(edits, edits[1:]):
        assert b1 <= a2, "overlapping edits"
    out = bytearray(s)
    for a, b, rep in reversed(edits): out[a:b] = rep
    doc.update_stream(xref, bytes(out))
    by = {}
    for fx, a, b, rep in stem_edits: by.setdefault(fx, []).append((a, b, rep))
    for fx, lst in by.items():
        fs = bytearray(doc.xref_stream(fx))
        for a, b, rep in sorted(lst, reverse=True): fs[a:b] = rep
        doc.update_stream(fx, bytes(fs))
        # widen the stem XObject BBox if the new stem is longer than the old one
        m = re.search(rb"q 1 0 0 1 (" + NUM + rb") " + NUM + rb" cm 0 0 m (" + NUM + rb") 0 l", bytes(fs))
        x1 = float(m.group(1)) + float(m.group(2)) + 0.5
        bb = [float(v) for v in doc.xref_get_key(fx, "BBox")[1].strip("[]").split()]
        if x1 > max(bb[0], bb[2]):
            if bb[2] >= bb[0]: bb[2] = x1
            else: bb[0] = x1
            doc.xref_set_key(fx, "BBox", "[%s]" % " ".join(f3(v) for v in bb))
    doc.save(A.out, garbage=3, deflate=True)
    for c in changes: print("  %-5s %-6s %-8s -> %-8s %s" % c)
    print("wrote", A.out)


if __name__ == "__main__":
    main()
