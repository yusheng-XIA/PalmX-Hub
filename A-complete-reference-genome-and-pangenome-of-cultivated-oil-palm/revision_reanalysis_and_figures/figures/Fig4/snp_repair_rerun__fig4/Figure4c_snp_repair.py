#!/usr/bin/env python3
"""Build Figure4_candidate.pdf from deliver/Main_Figures_revised/Figure4.pdf (read-only) changing ONLY panel 4c:
 - pi labels (3 significant digits x10^-3)
 - node circle radius: r^2 = r2min + (r2max-r2min)*norm(pi)   (decoded from the original: P2 10.398, P4 14.312; fits P1/P3 to 1e-3)
 - edge width: w = 0.85 + 2.056*norm(FST); edge colour = colour-bar strip at position 0.16+0.84*norm(FST)
   (decoded from the original: widths 0.85..2.906 linear in FST; colours match bar strips at 0.161..0.998)
 - colour-bar range = [FSTmin, FSTmax] (original: bar 54.707 pt spans 0.04365-0.12933); tick marks/labels moved
usage: f1_figure4c.py values.json out.pdf"""
import fitz, re, json, sys
from pathlib import Path
SCR = Path('${WORK_DIR}')
SRC = SCR / 'deliver/Main_Figures_revised/Figure4.pdf'
V = json.load(open(sys.argv[1])); OUT = sys.argv[2]
doc = fitz.open(SRC); X = 78
s = doc.xref_stream(X).decode('latin1')
def once(old, new):
    global s
    n = s.count(old)
    if n != 1: raise SystemExit('count %d for %r' % (n, old[:90]))
    s = s.replace(old, new)
pi = V['pi']; F = V['fst']            # pi in units of 1e-3 (floats); F keys 'Pop1_Pop2' ...
# ---- pi labels: original digit ops "(5)Tj .553 0 Td (.)Tj (6)Tj ... <37d7>Tj" ----
old_lab = {'Pop1': ('356.9502 605.5889', '5', '6', '37'), 'Pop4': ('453.0449 605.5889', '5', '9', '37'),
           'Pop3': ('453.0449 508.0762', '5', '7', '37'), 'Pop2': ('356.9502 508.0762', '4', '7', '30')}
for p, (pos, a, b, c) in old_lab.items():
    t = '%.2f' % pi[p]; assert len(t) == 4
    old = '%s Tm <028c>Tj/TT0 1 Tf .833 0 Td (=)Tj .747 0 Td (%s)Tj .553 0 Td (.)Tj (%s)Tj .163 Tc -.163 Tw .831 0 Td <%sd7>Tj' % (pos, a, b, c)
    new = '%s Tm <028c>Tj/TT0 1 Tf .833 0 Td (=)Tj .747 0 Td (%s)Tj .553 0 Td (.)Tj (%s)Tj .163 Tc -.163 Tw .831 0 Td <%02xd7>Tj' % (pos, t[0], t[2], ord(t[3]))
    once(old, new)
# ---- node circles ----
node = {'Pop1': (376.629, 586.7674, 13.469), 'Pop4': (472.723, 586.7724, 14.312), 'Pop3': (472.723, 527.2444, 13.773), 'Pop2': (376.629, 527.2444, 10.398)}
R2MIN, R2MAX = 10.398 ** 2, 14.312 ** 2
pmin, pmax = min(pi.values()), max(pi.values())
newr = {p: (R2MIN + (R2MAX - R2MIN) * (pi[p] - pmin) / (pmax - pmin)) ** 0.5 for p in pi}
def fmt(v):
    t = ('%.3f' % v).rstrip('0').rstrip('.')
    return t.replace('0.', '.', 1) if t.startswith('0.') else (t.replace('-0.', '-.', 1) if t.startswith('-0.') else t)
for p, (cx, cy, r) in node.items():
    top = '%s %s cm' % (('%.3f' % cx).rstrip('0'), ('%.4f' % (cy + r)).rstrip('0'))
    pat = re.compile(r'q 1 0 0 1 ' + re.escape(top) + r'( [^Q]*?)0 0 m ([-\d. c]+?) h? ?(S|f) Q')
    ms = list(pat.finditer(s)); assert len(ms) == 2, (p, len(ms), top)
    k = newr[p] / r
    for m in reversed(ms):
        path = re.sub(r'-?\d*\.?\d+', lambda q: fmt(float(q.group()) * k), m.group(2))
        rep = 'q 1 0 0 1 %s %s cm%s0 0 m %s %s%s Q' % (fmt(cx), fmt(cy + newr[p]), m.group(1), path, 'h ' if ' h ' in m.group(0) else '', m.group(3))
        s = s[:m.start()] + rep + s[m.end():]
# ---- colour-bar strips (read colours, keep bar untouched) ----
page = doc[0]
strips = sorted([((d['rect'].y0 + d['rect'].y1) / 2, d['fill']) for d in page.get_drawings()
                 if d['type'] == 'f' and abs(d['rect'].x0 - 494.0) < .2 and d['rect'].width < 5.5 and d['rect'].y1 < 100])
Y0, Y1 = 39.97, 94.68      # page coords top/bottom of bar
def colour_at(frac):       # frac 0 = bottom, 1 = top
    y = Y1 - frac * (Y1 - Y0)
    return min(strips, key=lambda t: abs(t[0] - y))[1]
fmin, fmax = min(F.values()), max(F.values())
edges = {'Pop1_Pop4': ('376.629 586.7714', '.957 .694 .635', '.85', '96.094 0'),
         'Pop1_Pop3': ('376.629 586.7714', '.929 .51 .455', '1.38', '96.094 -59.527'),
         'Pop1_Pop2': ('376.629 586.7714', '.525 .208 .286', '2.834', '0 -59.527'),
         'Pop3_Pop4': ('472.723 586.7714', '.878 .388 .392', '1.77', '0 -59.527'),
         'Pop2_Pop4': ('472.723 586.7714', '.498 .2 .282', '2.906', '-96.094 -59.527'),
         'Pop2_Pop3': ('376.629 527.2444', '.847 .31 .353', '2.032', '96.094 0')}
for k, (pos, col, w, to) in edges.items():
    nrm = (F[k] - fmin) / (fmax - fmin)
    c = colour_at(0.16 + 0.84 * nrm)
    nc = ' '.join(fmt(round(v, 3)) for v in c); nw = fmt(0.85 + 2.056 * nrm)
    once('q 1 0 0 1 %s cm %s RG 1 j %s w 0 0 m %s l S Q' % (pos, col, w, to), 'q 1 0 0 1 %s cm %s RG 1 j %s w 0 0 m %s l S Q' % (pos, nc, nw, to))
# ---- colour-bar ticks: form y of bar bottom/top ----
H = 620.787; yb, yt = H - Y1, H - Y0     # 526.107 .. 580.817  (orig bar re: 493.984 580.818 h -54.707)
yb, yt = 580.818 - 54.70703, 580.818
ticks = V['ticks']
old_ticks = ''.join('q 1 0 0 1 498.945 %s cm 2 J 1 j 0 0 m 1.7 0 l S Q ' % y for y in ('530.1624', '542.9354', '555.7054', '568.4784'))
new_ticks = ''.join('q 1 0 0 1 498.945 %s cm 2 J 1 j 0 0 m 1.7 0 l S Q ' % fmt(yb + (t - fmin) / (fmax - fmin) * (yt - yb)) for t in ticks)
once(old_ticks, new_ticks)
# labels: original baseline = tick y - 1.4983 (530.1624 -> 528.6641); step Td in text space (5.8 pt font)
lab_old = 'BT 5.8 0 0 5.8 501.6377 528.6641 Tm (0.05)Tj 0 2.202 Td (0.07)Tj 0 2.202 Td (0.09)Tj 0 2.202 Td (0.)Tj (1)Tj 1.317 0 Td (1)Tj '
ty = [yb + (t - fmin) / (fmax - fmin) * (yt - yb) - 1.4983 for t in ticks]
lab_new = 'BT 5.8 0 0 5.8 501.6377 %s Tm (%s)Tj ' % (fmt(ty[0]), '%.2f' % ticks[0])
for a, b, t in zip(ty, ty[1:], ticks[1:]):
    lab_new += '0 %s Td (%s)Tj ' % (fmt((b - a) / 5.8), '%.2f' % t)
once(lab_old, lab_new)
doc.update_stream(X, s.encode('latin1'))
doc.save(OUT, garbage=0, deflate=True)
print('radii', {p: round(v, 3) for p, v in newr.items()})
