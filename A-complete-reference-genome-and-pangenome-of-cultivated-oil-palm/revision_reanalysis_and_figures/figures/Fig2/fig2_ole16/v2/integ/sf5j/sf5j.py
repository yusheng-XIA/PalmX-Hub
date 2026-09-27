#!/usr/bin/env python3
"""Supplementary Fig. 5j: OLE16a upstream intergenic region across 39 haplotypes (five structural classes),
drawn from fix/ole16_cis/prom (summary.tsv, insertions_vs_FLHap2.tsv); also writes the Source Data table."""
from pathlib import Path
import csv
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Polygon
from matplotlib.lines import Line2D
HERE = Path(__file__).resolve().parent
S = HERE.parents[4]
P = S / 'fix/ole16_cis/prom'
plt.rcParams.update({'font.family': 'Arial', 'font.size': 6, 'pdf.fonttype': 42,
                     'mathtext.fontset': 'custom', 'mathtext.rm': 'Arial', 'mathtext.it': 'Arial:italic',
                     'mathtext.bf': 'Arial:bold', 'mathtext.sf': 'Arial:italic', 'axes.linewidth': 0.6, 'xtick.major.width': 0.6,
                     'xtick.major.size': 2.2, 'text.color': '#222222', 'axes.edgecolor': '#222222'})
NAME = {'African_hap2': 'FL-Hap2', 'American_hap1': 'FL-Hap1', 'BK_hap1': 'TN-Hap1', 'BK_hap2': 'TN-Hap2',
        'dura_hap1': 'TK-Hap1', 'dura_hap2': 'TK-Hap2', 'pisifera_hap1': 'NS-Hap1', 'pisifera_hap2': 'NS-Hap2',
        'nrly_hap1': 'Nigerian-Hap1', 'nrly_hap2': 'Nigerian-Hap2', 'MZ4_hap1': 'E. oleifera-Hap1',
        'MZ4_hap2': 'E. oleifera-Hap2'}
summ = {r['genome']: r for r in csv.DictReader(open(P / 'summary.tsv'), delimiter='\t')}
ins = {}
for r in csv.DictReader(open(P / 'insertions_vs_FLHap2.tsv'), delimiter='\t'):
    ins.setdefault(r['haplotype'], []).append(r)


def klass(g):
    v = ins.get(g, [])
    if not v:
        return 'EG-A'
    if len(v) == 2:
        return 'EG-C'
    s = int(v[0]['ins_start_rel_ATG'])
    if s < -18000:
        return 'EO'
    if -11400 < s < -11300:
        return 'EG-B'
    return 'EG-D'


rows = []
for g in summ:
    k = klass(g)
    iv = ins.get(g, [])
    rows.append([NAME.get(g, g), 'E. oleifera' if g.startswith('MZ4') or g == 'American_hap1' and False else
                 ('E. oleifera' if g.startswith('MZ4') else ('E. oleifera-derived (FL)' if g == 'American_hap1' else 'E. guineensis')),
                 k, summ[g]['upstream_intergenic_bp'],
                 ';'.join(f"{r['ins_start_rel_ATG']}..{r['ins_end_rel_ATG']}" for r in iv),
                 ';'.join(r['ins_len'] for r in iv)])
cnt = {}
for r in rows:
    cnt[r[2]] = cnt.get(r[2], 0) + 1
assert cnt == {'EO': 3, 'EG-A': 20, 'EG-B': 12, 'EG-C': 3, 'EG-D': 1}, cnt
order = {'EO': 0, 'EG-A': 1, 'EG-B': 2, 'EG-C': 3, 'EG-D': 4}
rows.sort(key=lambda r: (order[r[2]], r[0]))
with open(HERE / 'SF5j_OLE16a_upstream.tsv', 'w', newline='') as fh:
    w = csv.writer(fh, delimiter='\t')
    w.writerow(['Haplotype', 'Species origin', 'Upstream class', 'Upstream intergenic length (bp; start codon to the upstream neighbour gene)',
                'Insertion(s), start..end relative to the start codon (bp)', 'Insertion length(s) (bp)'])
    w.writerows(rows)

EOI = r'$\it{E.\ oleifera}$'; EGI = r'$\it{E.\ guineensis}$'
classes = [
    (EOI + '-type (n = 3)', r'both ' + EOI + ' haplotypes, FL-Hap1', -23.03, [(-18.45, -7.40, 'eo')]),
    (EGI + ' A (n = 20)', 'incl. FL-Hap2, TK-Hap1/2, TN-Hap1, Nigerian-Hap1', -11.95, []),
    (EGI + ' B (n = 12)', 'incl. TN-Hap2, NS-Hap1, Nigerian-Hap2', -21.38, [(-11.35, -1.91, 'b')]),
    (EGI + ' C (n = 3)', 'NS-Hap2, EG_071, EG_072', -34.03, [(-26.68, -15.60, 'rel'), (-12.91, -1.91, 'rel')]),
    (EGI + ' D (n = 1)', 'EG_025', -14.42, [(-5.82, -3.29, 'd')]),
]
col = {'eo': '#1b7f5f', 'b': '#9aa3ad', 'rel': '#6fb59b', 'd': '#c7ccd1'}
W, H = 510.0, 100.0
fig = plt.figure(figsize=(W / 72, H / 72))
ax = fig.add_axes([150 / W, 22 / H, 250 / W, 74 / H])
ax.set_xlim(-36, 1.4); ax.set_ylim(-len(classes) + 0.5, 0.55)
y = 0
for name, members, L, iv in classes:
    ax.plot([L, 0], [y, y], color='#444444', lw=0.7, zorder=1)
    ax.add_patch(Rectangle((L - 1.6, y - 0.16), 1.6, 0.32, fc='#d9d9d9', ec='#666666', lw=0.4))
    for a, b, c in iv:
        ax.add_patch(Rectangle((a, y - 0.2), b - a, 0.4, fc=col[c], ec='none', zorder=2))
        if c in ('eo', 'rel'):
            for x0 in (a, b - 1.9):
                ax.add_patch(Rectangle((x0, y - 0.2), 1.9, 0.4, fc='none', ec='#222222', lw=0.3, hatch='////', zorder=3))
    ax.add_patch(Polygon([[0, y - 0.22], [0, y + 0.22], [1.0, y]], fc='#c0392b', ec='none', zorder=3))
    for p in (-0.97, -0.16):
        ax.plot([p, p], [y + 0.2, y + 0.34], color='#b8860b', lw=0.6)
    fig.text(10 / W, ax.transData.transform((0, y))[1] / fig.dpi / H * 72 + 1.2 / H, name, ha='left', va='bottom',
             fontsize=5.5)
    fig.text(10 / W, ax.transData.transform((0, y))[1] / fig.dpi / H * 72 - 0.6 / H, members, ha='left', va='top',
             fontsize=5, color='#555555')
    y -= 1
ax.set_yticks([])
ax.set_xticks([-30, -20, -10, 0], ['−30', '−20', '−10', 'ATG'])
ax.spines['bottom'].set_bounds(-36, 0)
ax.set_xlabel('Distance upstream of the ' + r'$\mathsf{OLE16a}$' + ' start codon (kb)', labelpad=1.5)
for s in ('left', 'right', 'top'):
    ax.spines[s].set_visible(False)
h = [Rectangle((0, 0), 1, 1, fc=col['eo']), Rectangle((0, 0), 1, 1, fc=col['rel']), Rectangle((0, 0), 1, 1, fc=col['b']),
     Rectangle((0, 0), 1, 1, fc=col['d']), Rectangle((0, 0), 1, 1, fc='white', ec='#222222', lw=0.3, hatch='////'),
     Line2D([0], [0], color='#b8860b', lw=0.8), Rectangle((0, 0), 1, 1, fc='#d9d9d9', ec='#666666', lw=0.4),
     Line2D([0], [0], marker='>', color='none', mfc='#c0392b', mec='none', ms=4)]
lab = [EOI + '-type LTR-RT (11.1 kb)', 'Related LTR-RTs (91% identity)', 'Other insertion (9.4 kb)',
       'Other insertion (2.5 kb)', 'LTR (~2 kb)', 'RY elements', 'Upstream neighbour gene', r'$\mathsf{OLE16a}$']
fig.legend(h, lab, loc='upper left', bbox_to_anchor=(408 / W, 1.0), frameon=False, fontsize=5, handlelength=1.3,
           handleheight=0.8, labelspacing=0.35, borderaxespad=0.2)
fig.text(0, 1, 'j', fontsize=8, fontweight='bold', va='top', ha='left')
fig.savefig(HERE / 'SF5j.pdf'); fig.savefig(HERE / 'SF5j.png', dpi=300)
print(cnt)
