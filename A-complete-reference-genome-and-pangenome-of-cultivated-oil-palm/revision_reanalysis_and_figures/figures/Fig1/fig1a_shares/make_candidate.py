#!/usr/bin/env python3
"""Figure 1a candidate: update the vegetable-oil bar (39.5% -> 40.1%) under scheme A.
Edits the vector content stream of form XObject 142 in place (original embedded
Arial-BoldMT subset, same size/baseline; digits share one advance width, so the
right edge of the label is unchanged). Bar left edge fixed; width = share x track width."""
import fitz, pandas as pd
SRC = '../../deliver/Main_Figures_revised/Figure1.pdf'
t = pd.read_csv('Fig1a_shares_final.tsv', sep='\t')
r = t[(t.scheme == 'A_9crops_recommended') & (t.crop == 'Oil palm')].iloc[0]
area, oil = r.share_of_area_pct, r.share_of_oil_pct
TRACK = 142.56
doc = fitz.open(SRC); s = doc.xref_stream(142)
old_bar_a = b'364.757 712.962 12.320007 7.5599977 re f'
old_bar_o = b'364.757 689.395 56.329988 7.5599977 re f'
assert s.count(old_bar_a) == 1 and s.count(old_bar_o) == 1
assert s.count(b'(8.6%)Tj') == 1 and s.count(b'(39.5%)Tj') == 1
lab_a, lab_o = '%.1f%%' % area, '%.1f%%' % oil
print('area', area, lab_a, 'bar', TRACK*area/100, '(old 12.320007)')
print('oil ', oil, lab_o, 'bar', TRACK*oil/100, '(old 56.329988)')
if lab_a != '8.6%':
    s = s.replace(old_bar_a, b'364.757 712.962 %.6f 7.5599977 re f' % (TRACK*area/100))
    s = s.replace(b'(8.6%)Tj', b'(%s)Tj' % lab_a.encode())
if lab_o != '39.5%':
    s = s.replace(old_bar_o, b'364.757 689.395 %.6f 7.5599977 re f' % (TRACK*oil/100))
    s = s.replace(b'(39.5%)Tj', b'(%s)Tj' % lab_o.encode())
doc.update_stream(142, s)
doc.save('Figure1_candidate.pdf', garbage=0, deflate=True)
# verification + 600 dpi renders
for f in [SRC, 'Figure1_candidate.pdf']:
    d = fitz.open(f); p = d[0]
    for w in ['8.6%', '39.5%', '40.1%']:
        for b in p.search_for(w): print(f.split('/')[-1], w, [round(v, 3) for v in b])
    for dr in p.get_drawings():
        rr = dr['rect']
        if dr.get('fill') and abs(dr['fill'][0]-0.776) < 1e-3 and rr.y0 < 70 and rr.x0 > 330:
            print('  bar', [round(v, 3) for v in rr], 'width %.3f = %.3f%% of track 132.310' % (rr.width, 100*rr.width/132.31))
clip = fitz.Rect(330, 0, 481.44, 75)
for f, tag in [(SRC, 'original'), ('Figure1_candidate.pdf', 'candidate')]:
    fitz.open(f)[0].get_pixmap(dpi=600, clip=clip).save('fig1a_bars_%s_600dpi.png' % tag)
    fitz.open(f)[0].get_pixmap(dpi=600).save('Figure1_%s_600dpi.png' % tag)
