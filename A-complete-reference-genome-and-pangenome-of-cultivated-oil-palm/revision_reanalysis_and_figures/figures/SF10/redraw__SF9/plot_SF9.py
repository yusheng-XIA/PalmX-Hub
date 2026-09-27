#!/usr/bin/env python3
"""Supplementary Fig. 9 redraw (SV-GWAS for shell weight and flesh thickness).

Data: a/d SF9a/SF9d_manhattan.tsv and b/e SF9b/SF9e_*_QQ.tsv, built from the full EMMAX
results by build_SF9_full_tables.py; c/f Source Data Supp_SVGWAS_c/f (REVIEW workbook).
"""
import sys
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from f8_style import BLUE, DARK, GREY, MM, ORANGE, RED, TEAL, clean, letter, sci, sheet  # noqa: E402

HERE = Path(__file__).resolve().parent
N_TESTS = 370706
BONF = 0.05 / N_TESTS
SV_COL = {"INS": BLUE, "DEL": ORANGE, "MNV/COMPLEX": "#8E8E8E"}
SV_LAB = {"INS": "Insertion", "DEL": "Deletion", "MNV/COMPLEX": "MNV/complex"}

fig = plt.figure(figsize=(180 * MM, 200 * MM))
rows = [(0.965, 0.785), (0.705, 0.535), (0.455, 0.275), (0.195, 0.045)]


def manhattan(ax, name, title):
    m = pd.read_csv(HERE / name, sep="\t", dtype={"p": str})
    m["pos"] = m.pos.astype(float)
    m["nl"] = m.neglog10_p.astype(float)
    chroms = sorted(m.chrom.unique(), key=lambda c: int(c[3:5]))
    offs, ticks, off = {}, [], 0.0
    for c in chroms:
        offs[c] = off
        L = m.loc[m.chrom == c, "pos"].max()
        ticks.append(off + L / 2)
        off += L + 8e6
    m["x"] = m.pos + m.chrom.map(offs)
    for k, c in enumerate(chroms):
        s = m[m.chrom == c]
        if k % 2:
            ax.axvspan(s.x.min() - 4e6, s.x.max() + 4e6, color="#F3F4F6", lw=0, zorder=0)
    for t in ("MNV/COMPLEX", "DEL", "INS"):
        s = m[m.svtype == t]
        ax.scatter(s.x / 1e6, s.nl, s=2.2, color=SV_COL[t], lw=0, rasterized=True, label=SV_LAB[t], zorder=2)
    ax.axhline(-np.log10(BONF), color=RED, ls="--", lw=0.7, zorder=3)
    ax.text(off / 1e6, -np.log10(BONF) + 0.15, rf"Bonferroni $P = {sci(BONF, 4)}$", ha="right", va="bottom",
            fontsize=5.5, color=RED)
    ax.set_xticks(np.array(ticks) / 1e6, [c.replace("chr", "") for c in chroms], fontsize=5.5)
    ax.set_xlim(-5, off / 1e6)
    ax.set_ylim(0, m.nl.max() * 1.08)
    ax.set_xlabel("Chromosome")
    ax.set_ylabel(r"$-\log_{10}P$")
    ax.set_title(title, loc="left", pad=3)
    clean(ax)
    return m


def qq(ax, name):
    q = pd.read_csv(HERE / name, sep="\t")
    e, o = q.expected_neglog10P, q.observed_neglog10P
    lam = float(q.lambda_GC_full.iloc[0])
    ax.scatter(e, o, s=3, color=DARK, lw=0, rasterized=True)
    lim = max(e.max(), o.max()) * 1.05
    ax.plot([0, lim], [0, lim], color=RED, lw=0.7, ls="--")
    ax.set_xlim(0, e.max() * 1.1)
    ax.set_ylim(0, o.max() * 1.08)
    ax.set_xlabel(r"Expected $-\log_{10}P$")
    ax.set_ylabel(r"Observed $-\log_{10}P$")
    ax.set_title(rf"{int(q.full_test_count.iloc[0]):,} SV tests; $\lambda_{{\mathrm{{GC}}}}$ = {lam:.3f}",
                 loc="left", pad=3)
    clean(ax)
    return lam


def box(ax, name, groups, labels, colors, ylab, test):
    g = sheet(name)
    g["Phenotype"] = g.Phenotype.astype(float)
    g["d"] = g.Minor_allele_dosage.astype(int)
    data = [g.loc[g.d.isin(gr), "Phenotype"].values for gr in groups]
    bp = ax.boxplot(data, widths=0.5, showfliers=False, patch_artist=True, whis=1.5,
                    medianprops=dict(color=DARK, lw=1.0), whiskerprops=dict(lw=0.6), capprops=dict(lw=0.6),
                    boxprops=dict(lw=0.6))
    rng = np.random.default_rng(7)
    for i, (v, col) in enumerate(zip(data, colors), start=1):
        bp["boxes"][i - 1].set_facecolor(mpl.colors.to_rgba(col, 0.35))
        ax.scatter(i + rng.uniform(-0.15, 0.15, len(v)), v, s=5, color=col, alpha=0.8, lw=0, zorder=3)
    ax.set_xticks(range(1, len(data) + 1), [f"{l}\n(n = {len(v)})" for l, v in zip(labels, data)])
    ax.set_ylabel(ylab)
    if test == "mw":
        p = stats.mannwhitneyu(data[0], data[1], alternative="two-sided").pvalue
        t = rf"Two-sided Mann–Whitney $P = {sci(p)}$"
    else:
        p = stats.kruskal(*data).pvalue
        t = rf"Kruskal–Wallis $P = {sci(p)}$"
    ax.set_title(f"{g.Variant_ID.iloc[0]}\n" + t, loc="left", pad=3, fontsize=6.5)
    clean(ax)
    return [len(v) for v in data], p, g


# a-c shell weight
ax_a = fig.add_axes([0.08, rows[0][1], 0.905, rows[0][0] - rows[0][1]])
ma = manhattan(ax_a, "SF9a_manhattan.tsv", "Shell weight (n = 88)")
ax_a.legend(frameon=False, loc="upper center", markerscale=3, ncol=3, bbox_to_anchor=(0.5, 1.06),
            handletextpad=0.1, columnspacing=1.0)
ax_b = fig.add_axes([0.08, rows[1][1], 0.36, rows[1][0] - rows[1][1]])
lam_b = qq(ax_b, "SF9b_shell_weight_QQ.tsv")
ax_c = fig.add_axes([0.58, rows[1][1], 0.405, rows[1][0] - rows[1][1]])
n_c, p_c, gc = box(ax_c, "Supp_SVGWAS_c_genotypes", [[0], [1, 2]],
                   ["Major-allele homozygote", "Minor-allele carrier"], [BLUE, ORANGE], "Shell weight (g)", "mw")

# d-f flesh thickness
ax_d = fig.add_axes([0.08, rows[2][1], 0.905, rows[2][0] - rows[2][1]])
md = manhattan(ax_d, "SF9d_manhattan.tsv", "Flesh thickness (n = 140)")
ax_e = fig.add_axes([0.08, rows[3][1], 0.36, rows[3][0] - rows[3][1]])
lam_e = qq(ax_e, "SF9e_flesh_thickness_QQ.tsv")
ax_f = fig.add_axes([0.58, rows[3][1], 0.405, rows[3][0] - rows[3][1]])
n_f, p_f, gf = box(ax_f, "Supp_SVGWAS_f_genotypes", [[0], [1], [2]], ["Dosage 0", "Dosage 1", "Dosage 2"],
                   [BLUE, ORANGE, RED], "Flesh thickness (mm)", "kw")

for ax, s in [(ax_a, "a"), (ax_b, "b"), (ax_c, "c"), (ax_d, "d"), (ax_e, "e"), (ax_f, "f")]:
    letter(fig, ax, s, dx=-0.07, dy=0.015)
for ext in ("pdf", "png"):
    fig.savefig(HERE / f"Supplementary_Fig_09.{ext}", dpi=600, facecolor="white")

sig_a = (ma.nl > -np.log10(BONF)).sum()
sig_d = (md.nl > -np.log10(BONF)).sum()
print(f"Bonferroni {BONF:.4e}; lambda b {lam_b:.3f} e {lam_e:.3f}")
print(f"c n={n_c} MW P={p_c:.3e}; f n={n_f} KW P={p_f:.3e}")
print(f"significant markers shown: shell {sig_a}, flesh {sig_d}; max -log10P {ma.nl.max():.2f} / {md.nl.max():.2f}")
