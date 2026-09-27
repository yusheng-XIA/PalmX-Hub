#!/usr/bin/env python3
"""Figure 5h/5i candidate: switch the donor path from W = 0 to W = 4 (P = 15).

Base (read-only): scratchpad/deliver/Main_Figures_revised/Figure5.pdf (current W = 0 version).
Vector-only stream edits, applied identically to the two copies of the page content
(xref 63 = clipped panel-a copy, xref 8270 = visible b-i copy):
  1. 5h donor-mosaic bars: the segment rectangles of the three chromosomes whose path changed
     (chr01B, chr02B, chr13B) are regenerated from segments_W4.tsv. Geometry comes from a linear
     fit (x = x0 + k * Mb; segment ends at window multiples, as in the plotting script) to the 362 existing W = 0 rectangles; every donor keeps its existing colour
     (donor -> colour map read from the W = 0 bars), so the legend is unchanged.
  2. 5h right-hand labels 'dSV / dSNP   * c/t': digit-for-digit substitution inside the existing Tj
     strings (every new value has the same number of digits; Arial digits are equal-width).
  3. 5h and 5i triangles: fill colour switched gold <-> grey where the W = 4 mark differs.
  4. 5i red boxes (selected chr01B segments): the box list is rebuilt from the W = 4 chr01B
     segments (47 boxes); rows from the 5i donor order, x from a fit to the existing boxes.
  5. 5i red dSV/dSNP labels: labels of unchanged segments keep their exact (hand-placed) positions;
     labels of removed segments are dropped; labels of new segments are placed at the median
     offset (right of the box, same baseline convention) measured on the existing labels.
Asserts check every W = 0 element against the W = 0 tables before anything is changed.
"""
import csv, re, statistics
from collections import defaultdict
from pathlib import Path
import fitz
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
FIX = HERE.parent
BASE = FIX.parent / "deliver/Main_Figures_revised/Figure5.pdf"
OUT = HERE / "Figure5_candidate_W4.pdf"
GOLD, GREY = b".851 .643 .255", b".722 .773 .8"
ROW_Y = {247.28: "chr01B", 237.41: "chr02B", 227.53: "chr03B", 217.65: "chr04B", 207.78: "chr05B",
         197.9: "chr06B", 188.03: "chr07B", 178.15: "chr08B", 168.27: "chr09B", 158.4: "chr10B",
         148.52: "chr11B", 138.64: "chr12B", 128.77: "chr13B", 118.9: "chr14B", 109.02: "chr15B",
         99.14: "chr16B"}
BAR_TOP = {247.406: "chr01B", 237.531: "chr02B", 227.652: "chr03B", 217.777: "chr04B", 207.902: "chr05B",
           198.023: "chr06B", 188.148: "chr07B", 178.273: "chr08B", 168.394: "chr09B", 158.519: "chr10B",
           148.644: "chr11B", 138.765: "chr12B", 128.89: "chr13B", 119.015: "chr14B", 109.136: "chr15B",
           99.261: "chr16B"}
ZOOM_Y = 253.96
CHROMS = [f"chr{i:02d}B" for i in range(1, 17)]
ZC = "chr01B"

# ---------------- data
seg0 = pd.read_csv(FIX / "fig5hi/src/ideal_loadonly_segments.tsv", sep="\t")
seg0 = seg0[seg0.Panel == "African35"]
seg4 = pd.read_csv(HERE / "segments_W4.tsv", sep="\t")
by0 = pd.read_csv(FIX / "enh_C/data/ideal_loadonly_by_chrom.tsv", sep="\t").set_index("Chrom")
by4 = pd.read_csv(HERE / "by_chrom_W4.tsv", sep="\t").set_index("Chrom")
per0 = {r["Chrom"]: r for r in csv.DictReader(open(FIX / "fig5hi/per_chrom_capture.tsv"), delimiter="\t")}
per4 = {r["Chrom"]: r for r in csv.DictReader(open(HERE / "per_chrom_capture_W4.tsv"), delimiter="\t")}
marks = defaultdict(list)
for r in csv.DictReader(open(HERE / "position_marks_W4.tsv"), delimiter="\t"):
    marks[r["Chrom"]].append((int(r["Pos"]), int(r["Mark_W4"]), int(r["Mark_W0"])))
for c in marks: marks[c].sort()
CHANGED = sorted({c for c in CHROMS if not seg0[seg0.Chrom == c][["Start_Window_Index", "End_Window_Index_Exclusive", "Donor_ID"]]
                  .reset_index(drop=True).equals(seg4[seg4.Chrom == c][["Start_Window_Index", "End_Window_Index_Exclusive", "Donor_ID"]]
                  .reset_index(drop=True))}, key=CHROMS.index)
M = pd.read_csv(FIX / "enh_C/data/load_matrix_500kb.tsv", sep="\t")
chromlen = M.groupby("Chrom").Window_End_0based.max().to_dict()


def segs(df, c):
    return [tuple(r) for r in df[df.Chrom == c][["Start_Window_Index", "End_Window_Index_Exclusive", "Donor_ID",
                                                  "DSV_Count", "DSNP_Count"]].itertuples(index=False)]


def fmt(v):
    s = f"{v:.3f}".rstrip("0").rstrip(".")
    return s[1:] if s.startswith("0.") else ("-" + s[2:] if s.startswith("-0.") else s)


# ---------------- 1. 5h bars
SEG = re.compile(rb"(?:(?P<col>[\d.]+ [\d.]+ [\d.]+) scn )?(?:\.075 w )?(?P<x>[\d.]+) (?P<y>[\d.]+) (?P<w>[\d.]+) "
                 rb"(?P<h>-?[\d.]+) re (?P<op>[BfS]) ?")


def bars(s, log):
    a = s.find(b".396 .51 .702 scn .075 w 30.633 247.406"); assert a > 0 and s.count(b".075 w 30.633 247.406") == 1
    e = s.find(b"q 1 0 0 1", a); blk = s[a:e]
    items, pos = [], 0
    for m in SEG.finditer(blk):
        assert m.start() == pos; pos = m.end(); items.append(m)
    assert pos == len(blk)
    drawn, col, i = [], None, 0
    while i < len(items):
        m = items[i]
        if m["col"]: col = m["col"]
        if m["op"] == b"B":
            drawn.append(dict(top=float(m["y"]), x=float(m["x"]), w=float(m["w"]), h=m["h"], col=col,
                              raw=blk[m.start():m.end()], has_col=bool(m["col"]))); i += 1
        else:
            n = items[i + 1]; assert m["op"] == b"f" and n["op"] == b"S" and n["col"] is None
            drawn.append(dict(top=float(n["y"]), x=float(n["x"]), w=float(n["w"]), h=n["h"], col=col,
                              raw=blk[m.start():n.end()], has_col=bool(m["col"]))); i += 2
    rows = defaultdict(list)
    for d in drawn: rows[BAR_TOP[round(d["top"], 3)]].append(d)
    order = [BAR_TOP[round(d["top"], 3)] for d in drawn]
    assert order == sorted(order, key=CHROMS.index)          # rows are contiguous, in chromosome order
    colour, X, P = {}, [], []
    for c in CHROMS:
        s0 = segs(seg0, c); dr = rows[c]; assert len(s0) == len(dr), (c, len(s0), len(dr))
        for (st, en, don, *_), d in zip(s0, dr):
            assert colour.setdefault(don, d["col"]) == d["col"], don
            X.append(d["x"]); P.append(st * 0.5)
            X.append(d["x"] + d["w"]); P.append(en * 0.5)
    k, x0 = np.polyfit(P, X, 1); res = np.abs(np.array(X) - (k * np.array(P) + x0)).max()
    assert res < 0.01, res
    assert len(set(colour.values())) == len(colour) == 35
    out = []
    for c in CHROMS:
        if c not in CHANGED:
            for d in rows[c]:
                out.append(d["raw"] if d["has_col"] else d["col"] + b" scn " + d["raw"])
            continue
        top = rows[c][0]["top"]; h = rows[c][0]["h"].decode()
        for st, en, don, *_ in segs(seg4, c):
            xa = x0 + k * st * 0.5; xb = x0 + k * en * 0.5
            out.append(colour[don] + f" scn {fmt(xa)} {fmt(top)} {fmt(xb - xa)} {h} re B ".encode())
        log.append(("5h bars", c, len(rows[c]), len(segs(seg4, c))))
    new = b".075 w " + b"".join(o if o.endswith(b" ") else o + b" " for o in out)
    return s[:a] + new + s[e:], (k, x0, res), colour


# ---------------- 2. 5h labels
def relabel(s, log):
    start = s.find(b"5.8988 0 0 6 183.6191 242.708 Tm"); assert start > 0
    end = s.find(b"ET", start)
    TOK = re.compile(rb"\(((?:\\.|[^\\)])*)\)Tj|([\-\d.]+) ([\-\d.]+) Td")
    rows, cur = [], []
    for m in TOK.finditer(s, start, end):
        if m.group(1) is not None: cur.append((m.start(1), m.end(1), m.group(1)))
        elif abs(float(m.group(3))) > 1: rows.append(cur); cur = []
    rows.append(cur); assert len(rows) == 16, len(rows)
    repl = []
    for c, pieces in zip(CHROMS, rows):
        # the last row continues into the header 'dSV/dSNP  * captured/total' in the same BT block:
        # keep the shortest run of pieces that forms a complete 'a / b   * c/t' label
        for n_ in range(1, len(pieces) + 1):
            j_ = b"".join(p[2] for p in pieces[:n_]).decode()
            nxt = pieces[n_][2].decode() if n_ < len(pieces) else ""
            if re.fullmatch(r"\d+ / \d+\s+\*\s+\d+/\d+", j_) and not nxt[:1].isdigit():
                pieces = pieces[:n_]; break
        old = b"".join(p[2] for p in pieces).decode()
        m = re.fullmatch(r"(\d+) / (\d+)(\s+)\*(\s+)(\d+)/(\d+)", old); assert m, old
        o = (int(by0.loc[c, "Residual_DSV"]), int(by0.loc[c, "Residual_DSNP"]), int(per0[c]["Captured_W0"]), int(per0[c]["Fav_Total"]))
        assert tuple(int(m.group(i)) for i in (1, 2, 5, 6)) == o, (c, old, o)
        n = (int(by4.loc[c, "Residual_DSV"]), int(by4.loc[c, "Residual_DSNP"]), int(per4[c]["Captured_W4"]), int(per4[c]["Fav_Total"]))
        new = f"{n[0]} / {n[1]}{m.group(3)}*{m.group(4)}{n[2]}/{n[3]}"
        assert len(new) == len(old) and [ch.isdigit() for ch in new] == [ch.isdigit() for ch in old], (old, new)
        k = 0
        for a, b, t in pieces:
            seg = new[k:k + len(t)].encode(); k += len(t)
            if seg != t: repl.append((a, b, seg))
        log.append(("5h label", c, old, new))
    s2 = bytearray(s)
    for a, b, n in sorted(repl, reverse=True): s2[a:b] = n
    return bytes(s2)


# ---------------- 3. triangles
TRI = re.compile(rb"q 1 0 0 1 ([\d.\-]+) ([\d.\-]+) cm (\.851 \.643 \.255|\.722 \.773 \.8) scn")


def recolor(s, log):
    rows = defaultdict(list)
    for m in TRI.finditer(s):
        rows[round(float(m.group(2)), 2)].append((float(m.group(1)), m.start(3), m.group(3)))
    subs = []
    for y, chrom, panel in [(yy, c, "5h") for yy, c in ROW_Y.items()] + [(ZOOM_Y, ZC, "5i")]:
        tri = sorted(rows.get(y, [])); pos = marks.get(chrom, [])
        assert len(tri) == len(pos), (panel, chrom, len(tri), len(pos))
        if len(tri) >= 2:
            X = np.array([t[0] for t in tri]); P = np.array([p[0] for p in pos]) / 1e6
            a, b = np.polyfit(P, X, 1); assert np.abs(X - (a * P + b)).max() < 0.05
        for (x, off, col), (p, new, old) in zip(tri, pos):
            assert col == (GOLD if old else GREY), (panel, chrom, p)
            want = GOLD if new else GREY
            if want != col:
                subs.append((off, len(col), want))
                log.append((panel, chrom, p, round(x, 3), "gold" if old else "grey", "gold" if new else "grey"))
    for off, n, want in sorted(subs, reverse=True):
        s = s[:off] + want + s[off + n:]
    return s


# ---------------- 4/5. 5i boxes and labels
def zoom_rows():
    z = M[M.Chrom == ZC]
    donors = sorted(z.Sample_ID.unique())
    load = z.pivot(index="Sample_ID", columns="Window_Index", values="Total_Load").reindex(index=donors).to_numpy()
    order = np.argsort(load.sum(1))
    return {donors[i]: r for r, i in enumerate(order)}


BOX = re.compile(rb"([\d.]+) ([\d.]+) ([\d.]+) (-[\d.]+) re S ")
LAB_TOK = re.compile(rb"\(((?:\\.|[^\\)])*)\)Tj|([\-\d.]+) ([\-\d.]+) Td|([\-\d.]+) Tc|([\-\d.]+) Tw")


# new 5i labels that would collide with an existing label when put right of their box are placed
# above the box, at the baseline offset of the label of the named reference segment (row above)
LABEL_LIKE = {(13, 15, "EG_025"): (8, 9, "EG_025")}


def zoom(s, log, strict=True):
    hdr = b".784 .282 .282 SCN 0 J 0 j .5 w "
    a = s.find(hdr); assert a > 0 and s.count(hdr) == 1
    b0 = a + len(hdr); lab = s.find(b".722 .247 .243 scn BT", b0)
    boxes = []; pos = b0
    for m in BOX.finditer(s, b0, lab):
        assert m.start() == pos; pos = m.end()
        boxes.append(tuple(float(m.group(i)) for i in (1, 2, 3)) + (m.group(4),))
    assert pos == lab
    s0 = segs(seg0, ZC); s4 = segs(seg4, ZC); assert len(boxes) == len(s0) == 45
    row = zoom_rows()
    X, P, Yr, Yv = [], [], [], []
    for (st, en, don, *_), (x, y, w, h) in zip(s0, boxes):
        X += [x, x + w]; P += [st, en]
        Yr.append(row[don]); Yv.append(y)
    # window edge -> x: existing box edges are reused exactly; other edges from a linear fit
    # (the right edge of the last window, 513.599, is the trimmed axis end and is kept as drawn)
    edge = defaultdict(list)
    for p_, x_ in zip(P, X): edge[p_].append(x_)
    fitP = [p_ for p_ in P if p_ != 354]; fitX = [x_ for p_, x_ in zip(P, X) if p_ != 354]
    kx, bx = np.polyfit(fitP, fitX, 1); rx = np.abs(np.array(fitX) - (kx * np.array(fitP) + bx)).max()
    ky, by = np.polyfit(Yr, Yv, 1); ry = np.abs(np.array(Yv) - (ky * np.array(Yr) + by)).max()
    assert rx < 0.02 and ry < 0.01, (rx, ry)
    hmap = {round(y, 3): h for x, y, w, h in boxes}
    yrow = {}
    for r_, y_ in zip(Yr, Yv): yrow.setdefault(r_, y_)

    def ex(p_):
        return float(np.mean(edge[p_])) if p_ in edge else bx + kx * p_

    def box_of(st, en, don):
        xa = ex(st); xb = ex(en)
        y = yrow.get(row[don], by + ky * row[don])
        return xa, y, xb - xa
    newbox = b""
    for st, en, don, *_ in s4:
        xa, y, w = box_of(st, en, don)
        h = hmap.get(round(y, 3), b"-6.227005")
        newbox += f"{fmt(xa)} {fmt(y)} {fmt(w)} ".encode() + h + b" re S "
    # ---- labels
    bt = lab + len(b".722 .247 .243 scn BT ")
    mtm = re.match(rb"([\d.]+) 0 0 ([\d.]+) ([\d.]+) ([\d.]+) Tm ", s[bt:]); assert mtm
    sx, sy, tx, ty = (float(mtm.group(i)) for i in (1, 2, 3, 4))
    et = s.find(b"ET", bt)
    units = []; ux = uy = 0.0; tc = tw = 0.0
    units.append(dict(u=(0.0, 0.0), tj=[]))
    for m in LAB_TOK.finditer(s, bt + mtm.end(), et):
        if m.group(1) is not None: units[-1]["tj"].append((tc, tw, m.group(1)))
        elif m.group(2) is not None:
            ux += float(m.group(2)); uy += float(m.group(3)); units.append(dict(u=(ux, uy), tj=[]))
        elif m.group(4) is not None: tc = float(m.group(4))
        else: tw = float(m.group(5))
    labels = []
    for u in units:
        prev = labels[-1][-1] if labels else None
        exp = sum(0.333 if ch == "/" else 0.501 for t in prev["tj"] for ch in t[2].decode()) if prev else 0
        if prev and abs(u["u"][1] - prev["u"][1]) < 1e-6 and abs(u["u"][0] - prev["u"][0] - exp) < 0.06:
            labels[-1].append(u)
        else:
            labels.append([u])
    def ltext(L): return b"".join(t[2] for u in L for t in u["tj"]).decode()
    def lxy(L): return tx + sx * L[0]["u"][0], ty + sy * L[0]["u"][1]
    # xref 8270 (visible) has 45 labels; the clipped copy xref 63 (not visible: its form BBox shows only
    # panel a) is an older state lacking '0/103' and '0/17' -- tolerated there, those two stay absent.
    assert len(labels) == 45 or (strict is False and len(labels) == 43), len(labels)
    # match labels to W0 segments (same text, nearest box)
    seginfo0 = [(sg, box_of(*sg[:3])) for sg in s0]
    free = set(range(len(seginfo0))); lab_of = {}
    cand = []
    for li, L in enumerate(labels):
        lx, ly = lxy(L)
        for si, (sg, (xa, y, w)) in enumerate(seginfo0):
            if f"{sg[3]}/{sg[4]}" == ltext(L):
                cand.append((abs(lx - (xa + w)) + abs(ly - (y - 6.227)), li, si))
    used_l = set()
    for d, li, si in sorted(cand):
        if li in used_l or si not in free: continue
        used_l.add(li); free.discard(si); lab_of[si] = li
    assert len(lab_of) == len(labels) and (not free or strict is False)
    offs = []
    for si, li in lab_of.items():
        sg, (xa, y, w) = seginfo0[si]; lx, ly = lxy(labels[li])
        offs.append((lx - (xa + w), ly - y))
    # the "right of box" convention: label starts after the box right edge, baseline within the row
    right = [o for o in offs if 0 <= o[0] < 3 and -6.3 < o[1] < 0]
    dx = statistics.median(o[0] for o in right); dy = statistics.median(o[1] for o in right)
    keep = {tuple(sg[:3]): li for si, li in lab_of.items() for sg in [seginfo0[si][0]]}
    s4keys = {tuple(sg[:3]) for sg in s4}
    out = []   # (abs_x, abs_y, [(tc,tw,str)...]) per unit
    for sg in s4:
        k3 = tuple(sg[:3])
        unchanged = any(tuple(x) == tuple(sg) for x in s0)
        if unchanged:
            for u in (labels[keep[k3]] if k3 in keep else []):
                out.append((tx + sx * u["u"][0], ty + sy * u["u"][1], u["tj"]))
        else:
            xa, y, w = box_of(*sg[:3]); lx, ly = xa + w + dx, y + dy
            ref = LABEL_LIKE.get(k3)
            if ref is not None:          # place like the label of a neighbouring segment (above the box)
                rx_, ry_, _ = box_of(*ref); rl = lxy(labels[keep[ref]])
                rtxt = ltext(labels[keep[ref]]); rend = rl[0] + sx * (sum(0.333 if ch == "/" else 0.556 for ch in rtxt))
                lx, ly = max(xa, rend + 1.2), y + (rl[1] - ry_)   # after the reference label, same baseline offset
            txt = f"{sg[3]}/{sg[4]}"
            adv = {"/": 0.333}; cx = 0.0
            for ch in txt:
                out.append((lx + sx * cx, ly, [(0.0, 0.0, ch.encode())])); cx += adv.get(ch, 0.501)
            log.append(("5i label added", ZC, f"{sg[0]}-{sg[1]} {sg[2]}", txt, round(lx, 2), round(ly, 2)))
    for sg in s0:
        if tuple(sg[:3]) not in s4keys or next(x for x in s4 if tuple(x[:3]) == tuple(sg[:3]))[3:] != sg[3:]:
            log.append(("5i label removed", ZC, f"{sg[0]}-{sg[1]} {sg[2]}", f"{sg[3]}/{sg[4]}", "", ""))
    # re-emit the BT block: first unit via Tm, then relative Td in text units (rounded, drift-free)
    body = b""; cu = (0.0, 0.0); ctc = ctw = 0.0; first = True
    for axx, ayy, tj in out:
        ux_, uy_ = (axx - tx) / sx, (ayy - ty) / sy
        if first:
            nt = f"{fmt(sx)} 0 0 {fmt(sy)} {fmt(axx)} {fmt(ayy)} Tm ".encode(); cu0 = (ux_, uy_)
            body += nt; cu = (round(ux_, 3), round(uy_, 3)); base = (axx, ayy); first = False
            ux_r, uy_r = 0.0, 0.0; cur = (0.0, 0.0)
        else:
            rx_ = (axx - base[0]) / sx; ry_ = (ayy - base[1]) / sy
            ddx = round(rx_ - cur[0], 3); ddy = round(ry_ - cur[1], 3)
            body += f"{fmt(ddx)} {fmt(ddy)} Td ".encode(); cur = (cur[0] + ddx, cur[1] + ddy)
        for tc_, tw_, t in tj:
            if tc_ != ctc: body += f"{fmt(tc_)} Tc ".encode(); ctc = tc_
            if tw_ != ctw: body += f"{fmt(tw_)} Tw ".encode(); ctw = tw_
            body += b"(" + t + b")Tj "
    if ctc: body += b"0 Tc "
    if ctw: body += b"0 Tw "
    s = s[:b0] + newbox + b".722 .247 .243 scn BT " + body + s[et:]
    log.append(("5i boxes", ZC, len(boxes), len(s4), f"dx={dx:.3f} dy={dy:.3f}", ""))
    return s


def main():
    doc = fitz.open(BASE)
    logs = {}
    for xref in (63, 8270):
        s = doc.xref_stream(xref); lg = []
        s, fit, colour = bars(s, lg)
        s = relabel(s, lg)
        s = recolor(s, lg)
        s = zoom(s, lg, strict=(xref == 8270))
        doc.update_stream(xref, s); logs[xref] = lg
        print(f"xref {xref}: bar fit k={fit[0]:.5f} pt/Mb x0={fit[1]:.3f} res={fit[2]:.4f}; "
              f"triangles {sum(1 for l in lg if l[0] in ('5h', '5i'))}; labels "
              f"{sum(1 for l in lg if l[0] == '5h label' and l[2] != l[3])}")
    assert [l for l in logs[63]] == [l for l in logs[8270]]
    doc.save(OUT, garbage=3, deflate=True)
    with open(HERE / "figure_edit_log_W4.tsv", "w") as fh:
        fh.write("item\tchrom\tpos_or_old\tx_or_new\told\tnew\n")
        for r in logs[8270]: fh.write("\t".join(map(str, r)) + "\n")
    print("changed chromosomes:", CHANGED, "saved", OUT)


if __name__ == "__main__":
    main()
