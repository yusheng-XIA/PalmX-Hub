#!/usr/bin/env python3
"""Candidate ED/Supplementary panel: dSV robustness within E. guineensis (African35)."""
import sys
import numpy as np, pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import stats

D = sys.argv[1]          # local dir with pulled results
TAG = sys.argv[2]        # callability source used for the oleifera state (mm2 or syri)
MAIN = sys.argv[3] if len(sys.argv) > 3 else 'EG35_c_Phoenix_and_Oleifera_k1'
plt.rcParams.update({'font.family': 'Arial', 'font.size': 6, 'axes.linewidth': 0.5, 'xtick.major.width': 0.5,
                     'ytick.major.width': 0.5, 'xtick.major.size': 2, 'ytick.major.size': 2, 'pdf.fonttype': 42,
                     'axes.spines.top': False, 'axes.spines.right': False})
S = pd.read_csv(f'{D}/coloc_{TAG}/coloc_summary.tsv', sep='\t').set_index('label')
V = pd.read_csv(f'{D}/coloc_syri/coloc_summary.tsv', sep='\t').set_index('label')
fig = plt.figure(figsize=(180 / 25.4, 105 / 25.4))
gs = fig.add_gridspec(2, 6, height_ratios=[1, 1.05], hspace=0.62, wspace=1.3, left=0.07, right=0.98, top=0.93, bottom=0.1)
ax = [fig.add_subplot(gs[0, 0:2]), fig.add_subplot(gs[0, 2:4]), fig.add_subplot(gs[0, 4:6]), fig.add_subplot(gs[1, 1:4])]
ax_d2 = fig.add_subplot(gs[1, 4:6], sharey=ax[3])
sch = [('a_Phoenix', 'Phoenix', '#4C72B0'), ('b_Oleifera', 'E. oleifera', '#DD8452'),
       ('c_Phoenix_and_Oleifera', 'Phoenix + E. oleifera', '#55A868')]

# a: numbers of dSVs
a = ax[0]; x = np.arange(3); wd = 0.26
for i, (s, lab, col) in enumerate(sch):
    y = [S.loc[f'EG35_{s}_k{k}', 'n_dSV'] for k in (1, 2, 3)]
    b = a.bar(x + (i - 1) * wd, y, wd, color=col, label=lab, lw=0)
    for xx, yy in zip(x + (i - 1) * wd, y):
        a.text(xx, yy + 25, f'{yy:,}', ha='center', va='bottom', fontsize=4.5, rotation=90)
a.axhline(1480, color='0.4', ls='--', lw=0.5, label='All38 published (1,480)')
a.set_xticks(x); a.set_xticklabels(['1/35', '≤2/35', '≤3/35']); a.set_xlabel('Maximum carriers (E. guineensis)')
a.set_ylabel('Candidate dSVs'); a.set_ylim(0, 2700); a.set_xlim(-0.5, 2.5)
leg = a.legend(frameon=False, fontsize=5, loc='upper left', handlelength=0.8, title='Derived-state outgroup', title_fontsize=5, ncol=2, bbox_to_anchor=(0, 1.04), columnspacing=0.8)
leg._legend_box.align = 'left'

# b: 5e scatter for the main set
b = ax[1]
cnt = pd.read_csv(f'{D}/coloc_{TAG}/{MAIN}.fig_e_counts.tsv', sep='\t')
b.scatter(cnt.dSV_Count, cnt.dSNP_Count, s=3, color='#55A868', alpha=0.6, lw=0)
sl = stats.linregress(cnt.dSV_Count, cnt.dSNP_Count); xx = np.linspace(0, cnt.dSV_Count.max(), 10)
b.plot(xx, sl.intercept + sl.slope * xx, color='k', lw=0.6)
r = S.loc[MAIN]
b.text(0.03, 0.97, f"r = {r.e_r:.2f}\nP = {r.e_P:.1e}\nn = {int(r.e_N)}", transform=b.transAxes, va='top', fontsize=5)
b.set_xlabel('dSV count'); b.set_ylabel('dSNP count')
b.set_title('E. guineensis haplotype × chromosome', fontsize=6)

# c: 5f profile
c = ax[2]
pr = pd.read_csv(f'{D}/coloc_{TAG}/{MAIN}.fig_f_profile.tsv', sep='\t')
c.plot(pr.bin_mid_kb, pr.Observed, color='#C44E52', lw=0.6, label='dSV carrier')
c.plot(pr.bin_mid_kb, pr.Background, color='0.5', lw=0.6, label='Matched background')
c.set_yscale('symlog', linthresh=1); c.set_xlabel('Distance to focal dSV (kb)'); c.set_ylabel('dSNPs per 10-kb bin')
c.legend(frameon=False, fontsize=5, loc='upper right', handlelength=1, bbox_to_anchor=(1.02, 0.84))
c.set_ylim(1, 3000); c.text(0.03, 0.97, f"n = {int(r.f_pairs)} dSV–carrier pairs\nPaired Wilcoxon P = {r.f_wilcoxon_P:.1e}", transform=c.transAxes, fontsize=5, va='top')

# d: summary across definitions
d = ax[3]
rows = [('All38 (published)', V.loc['V0_All38_formal_1480']), ('African35, Phoenix (924)', V.loc['V1_African35_924_dSNP70269'])]
for k in (1, 2, 3):
    for s, lab, col in sch:
        rows.append((f'{lab}, ≤{k}/35' if k > 1 else f'{lab}, 1/35', S.loc[f'EG35_{s}_k{k}']))
ys = np.arange(len(rows))[::-1]
rr = [q.e_r for _, q in rows]
ratio = [q.f_obs_mean / q.f_bg_mean for _, q in rows]
d.scatter(rr, ys, s=8, color='k', zorder=3)
d.set_yticks(ys); d.set_yticklabels([l for l, _ in rows], fontsize=5)
d.set_xlim(0.6, 1.0); d.set_xlabel('Pearson r, dSV vs dSNP (Fig. 5e)')
d.axvline(0.9009, color='0.6', lw=0.5, ls='--')
ax_d2.scatter(ratio, ys, s=8, color='#C44E52', zorder=3); ax_d2.set_xlim(0, 9.5)
ax_d2.axvline(1, color='0.6', lw=0.5, ls='--'); ax_d2.tick_params(labelleft=False)
ax_d2.set_xlabel('Observed / background dSNPs (±1 Mb; Fig. 5f)')
for y, q in zip(ys, rows):
    ax_d2.text(7.0, y, f"P={q[1].f_wilcoxon_P:.0e}", fontsize=4.5, va='center')
fig.text(0.015, 0.44, 'd', fontsize=8, fontweight='bold', va='top')
for i, lab in enumerate('abc'):
    ax[i].text(-0.22, 1.1, lab, transform=ax[i].transAxes, fontsize=8, fontweight='bold', va='top')
fig.savefig(f'{D}/enhB_candidate_panel_{TAG}.pdf'); fig.savefig(f'{D}/enhB_candidate_panel_{TAG}.png', dpi=450)
print('saved')
