"""Fig. 4c candidates: pi / mean FST for the K = 3, 4 and 8 ADMIXTURE partitions (argmax ancestry; 308-based SNP set).
Variant A: three networks side by side.  Variant B: K = 3 and K = 4 networks + K = 8 FST heatmap with pi strip.
Shared FST colour scale across all partitions; node colours as the K bars of Fig. 4a. Panel box 175 x 100.5 pt."""
import sys, itertools, numpy as np, pandas as pd, matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path
H = Path(__file__).resolve().parent; D = H / "data"; OUT = H / "out"
mpl.rcParams.update({"font.family": "Arial", "pdf.fonttype": 42, "font.size": 5.5, "axes.linewidth": 0.4})
PT = 1 / 72; WC = (175.0, 100.5)
sys.path.insert(0, str(H)); from fig4_palette import KC
def load(K):
    if K == 4:
        pf = pd.read_csv(D / "pi_fst_genomewide.tsv", sep="\t"); pf = pf[pf.snp_set == "present"]
        pi = {int(g[-1]): 1000 * v for g, v in pf[(pf.metric == "pi") & pf.estimator.str.startswith("missing-aware")][["group_or_pair", "value"]].values}
        F = {(int(k[6]), int(k[-1])): v for k, v in pf[pf.metric == "fst"][["group_or_pair", "value"]].values}
    else:
        t = pd.read_csv(D / f"pi_fst_K{K}.tsv", sep="\t")
        pi = {int(g.split("Pop")[1]): 1000 * v for g, v in t[t.metric == "pi"][["group_or_pair", "value"]].values}
        F = {tuple(int(x) for x in re_split(k)): v for k, v in t[t.metric == "fst"][["group_or_pair", "value"]].values}
    n = pd.read_csv(D / f"groups_K{K}_assignments.tsv", sep="\t")["pop"].value_counts().to_dict()
    return pi, F, n
def re_split(k): return [p.split("_")[0] for p in k.split("Pop")[1:]]
DAT = {K: load(K) for K in (3, 4, 8)}
allF = [f for K in DAT for f in DAT[K][1].values()]; allpi = [p for K in DAT for p in DAT[K][0].values()]
lo, hi = np.floor(min(allF) * 100) / 100, np.ceil(max(allF) * 100) / 100
cmap = mpl.colors.LinearSegmentedColormap.from_list("f", ["#F8D4C8", "#F08A78", "#D64A58", "#7F3348"]); norm = mpl.colors.Normalize(lo, hi)
PLO, PHI = min(allpi), max(allpi)
def rad(p, rmin, rmax): return np.sqrt(rmin ** 2 + (rmax ** 2 - rmin ** 2) * (p - PLO) / (PHI - PLO))   # node area linear in pi
def network(ax, K, rmin, rmax, lwmax, label_pi):
    pi, F, n = DAT[K]
    if K == 4: pos = {1: (0.2, 0.8), 4: (0.8, 0.8), 2: (0.2, 0.2), 3: (0.8, 0.2)}
    elif K == 3: pos = {1: (0.2, 0.75), 2: (0.8, 0.75), 3: (0.5, 0.2)}
    else: pos = {i: (0.5 + 0.38 * np.sin(2 * np.pi * (i - 1) / K), 0.5 + 0.38 * np.cos(2 * np.pi * (i - 1) / K)) for i in range(1, K + 1)}
    for (a, b), f in sorted(F.items(), key=lambda x: x[1]):
        ax.plot(*zip(pos[a], pos[b]), color=cmap(norm(f)), lw=0.4 + lwmax * (f - lo) / (hi - lo), zorder=1, solid_capstyle="butt")
    for p, (x, y) in pos.items():
        ax.add_patch(mpl.patches.Circle((x, y), rad(pi[p], rmin, rmax), color=KC[K][p - 1], zorder=2, lw=0))
        ax.text(x, y, str(p), ha="center", va="center", fontsize=5.5, zorder=3)
        if label_pi:
            ax.text(x, y + (1 if y > 0.5 else -1) * (rmax + 0.07), f"{pi[p]:.2f}", ha="center", va="center", fontsize=5)
    ax.set_xlim(-0.08, 1.08); ax.set_ylim(-0.12, 1.12); ax.set_aspect("equal"); ax.axis("off")
    ax.set_title(f"K = {K}", fontsize=5.5, pad=1)
def cbar(fig, rect):
    cax = fig.add_axes(rect); cb = mpl.colorbar.ColorbarBase(cax, cmap=cmap, norm=norm, ticks=[t for t in (0.04, 0.08, 0.12, 0.16) if lo <= t <= hi])
    cb.outline.set_linewidth(0.4); cax.tick_params(labelsize=5, width=0.4, length=2); cax.set_title("Mean\n$F_{ST}$", fontsize=5, pad=2)
# ---- variant A
fig = plt.figure(figsize=(WC[0] * PT, WC[1] * PT))
for i, (K, r0, r1, lw) in enumerate([(3, 0.09, 0.15, 2.4), (4, 0.08, 0.14, 2.4), (8, 0.055, 0.09, 1.6)]):
    network(fig.add_axes([0.005 + i * 0.285, 0.10, 0.28, 0.80]), K, r0, r1, lw, K != 8)
cbar(fig, [0.885, 0.25, 0.03, 0.45])
fig.text(0.44, 0.03, "node size, π (× 10$^{-3}$; values shown for K = 3, 4); numbers, group", ha="center", fontsize=5)
fig.savefig(OUT / "Fig4c_multi_A_networks.pdf"); fig.savefig(OUT / "Fig4c_multi_A_networks.png", dpi=600); plt.close(fig)
# ---- variant B (chosen 27 Sep; also written as out/Fig4c_snp_repair.pdf)
fig = plt.figure(figsize=(WC[0] * PT, WC[1] * PT))
network(fig.add_axes([0.0, 0.12, 0.29, 0.78]), 3, 0.09, 0.15, 2.4, True)
network(fig.add_axes([0.28, 0.12, 0.29, 0.78]), 4, 0.08, 0.14, 2.4, True)
pi8, F8, n8 = DAT[8]
M = np.full((8, 8), np.nan)
for (a, b), f in F8.items(): M[a - 1, b - 1] = M[b - 1, a - 1] = f
ax = fig.add_axes([0.62, 0.20, 0.23, 0.23 * WC[0] / WC[1]])
ax.imshow(M, cmap=cmap, norm=norm); ax.set_xticks(range(8)); ax.set_yticks(range(8))
ax.set_xticklabels(range(1, 9), fontsize=5); ax.set_yticklabels(range(1, 9), fontsize=5); ax.tick_params(length=0, pad=1)
for i in range(8): ax.add_patch(mpl.patches.Rectangle((i - 0.5, i - 0.5), 1, 1, color=KC[8][i], lw=0))
ax.set_title("K = 8", fontsize=5.5, pad=1)
axp = fig.add_axes([0.62, 0.06, 0.23, 0.08]); axp.bar(range(8), [pi8[i] for i in range(1, 9)], color=KC[8], width=0.8)
axp.set_xlim(-0.5, 7.5); axp.set_ylim(0, max(pi8.values()) * 1.1); axp.set_xticks([]); axp.tick_params(labelsize=5, length=1.5, width=0.4)
axp.set_ylabel("π", fontsize=5, labelpad=1); [axp.spines[s].set_visible(False) for s in ("top", "right")]
cbar(fig, [0.875, 0.25, 0.03, 0.45])
fig.savefig(OUT / "Fig4c_multi_B_heatmap.pdf"); fig.savefig(OUT / "Fig4c_snp_repair.pdf"); fig.savefig(OUT / "Fig4c_multi_B_heatmap.png", dpi=600); plt.close(fig)
print("FST range", lo, hi, "pi range", PLO, PHI)
