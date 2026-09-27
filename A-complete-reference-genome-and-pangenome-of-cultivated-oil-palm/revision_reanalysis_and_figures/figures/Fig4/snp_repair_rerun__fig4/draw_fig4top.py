"""Fig. 4a-c from the repaired SNP call set, drawn at the size of the original panel boxes (Arial, 5-6 pt).
a: ADMIXTURE K = 3, 4, 8 (new Q, columns matched to the published colours) ordered by the published dendrogram;
b: PCA (PLINK, new LD-pruned set; signs oriented to the published PCs; % of GRM trace); c: pi / mean FST network."""
import sys, numpy as np, pandas as pd, matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse, Rectangle
from Bio import Phylo
from pathlib import Path
H = Path(__file__).resolve().parent; D = H / "data"; OUT = H / "out"; OUT.mkdir(exist_ok=True)
sys.path.insert(0, str(H.parents[1] / "beautify/common")); import palA
sys.path.insert(0, str(H))
mpl.rcParams.update({"font.family": "Arial", "pdf.fonttype": 42, "font.size": 5.5, "axes.linewidth": 0.5,
                     "xtick.major.width": 0.4, "ytick.major.width": 0.4, "xtick.major.size": 2, "ytick.major.size": 2,
                     "xtick.labelsize": 5.5, "ytick.labelsize": 5.5})
PT = 1 / 72
K4C = {f"K4_Pop{i}": c for i, c in enumerate(["#70B4E7", "#F3766A", "#82C7B8", "#9D7BCB"], 1)}
K3C = ["#70B4E7", "#F2B76A", "#82C7B8"]
from fig4_palette import K8C
ids = [l.split()[1] for l in open(D / "all.fam")]
Q = {k: pd.read_csv(D / f"all.{k}.Q", sep=r"\s+", header=None, index_col=None).set_index(pd.Index(ids)) for k in (3, 4, 8)}
A = pd.read_csv(D / "K4_dominant_assignments.tsv", sep="\t").set_index("sample")
order = [s for s in (l.strip() for l in open(D / "integrated_rooted_tree_structure_sample_order.txt")) if s in ids]
order += [s for s in ids if s not in set(order)]
grp = A.dominant_group.to_dict(); cnt = A.dominant_group.value_counts().to_dict()
W = {"a": (170.0, 100.5), "b": (170.5, 100.5), "c": (175.0, 100.5)}

# ---------------------------------------------------------------- a
fig = plt.figure(figsize=(W["a"][0] * PT, W["a"][1] * PT))
L, R = 0.155, 0.995
ax_t = fig.add_axes([L, 0.745, R - L, 0.235]); axs = [fig.add_axes([L, 0.60 - i * 0.135, R - L, 0.13]) for i in range(3)]
ax_s = fig.add_axes([L, 0.305, R - L, 0.028])
x = {s: i for i, s in enumerate(order)}; n = len(order)
tree = Phylo.read(D / "oilpalm_newick.txt", "newick")
for tip in list(tree.get_terminals()):
    if tip.name not in x: tree.prune(tip)
nx, ny = {}, {}
def rec(c):
    if c.is_terminal(): nx[c] = x[c.name]; ny[c] = 0.0; return
    for ch in c.clades: rec(ch)
    nx[c] = float(np.mean([nx[ch] for ch in c.clades])); ny[c] = max(ny[ch] for ch in c.clades) + 1
rec(tree.root)
for p in tree.find_clades(order="preorder"):
    for ch in p.clades:
        ax_t.plot([nx[ch], nx[ch], nx[p]], [ny[ch], ny[p], ny[p]], color="#6E7479", lw=0.25, solid_capstyle="butt")
# background: local majority of K4 groups along the tree order (as published)
labs = [grp.get(s) for s in order]
maj = []
for i in range(n):
    w = [l for l in labs[max(0, i - 20): i + 21] if l]; maj.append(max(set(w), key=w.count))
s0 = 0
for i in range(1, n + 1):
    if i == n or maj[i] != maj[s0]:
        ax_t.axvspan(s0 - 0.5, i - 0.5, color=K4C[maj[s0]], alpha=0.18, lw=0, zorder=-5); s0 = i
ax_t.set_xlim(-0.5, n - 0.5); ax_t.set_ylim(0, max(ny.values()) * 1.02); ax_t.axis("off")
for ax, k, cols in zip(axs, (3, 4, 8), (K3C, None, K8C)):
    q = Q[k].loc[order].to_numpy(); bottom = np.zeros(n)
    for j in range(k):
        c = K4C[f"K4_Pop{j + 1}"] if k == 4 else cols[j]
        ax.bar(np.arange(n), q[:, j], bottom=bottom, width=1.0, color=c, lw=0); bottom += q[:, j]
    ax.set_xlim(-0.5, n - 0.5); ax.set_ylim(0, 1); ax.set_xticks([]); ax.set_yticks([])
    for sp in ax.spines.values(): sp.set_visible(False)
    ax.text(-0.012, 0.5, f"K={k}", transform=ax.transAxes, ha="right", va="center", fontsize=5.5)
ax_s.bar(np.arange(n), 1, width=1.0, color=[K4C[g] for g in labs], lw=0); ax_s.set_xlim(-0.5, n - 0.5); ax_s.axis("off")
hs = [Rectangle((0, 0), 1, 1, color=K4C[f"K4_Pop{i}"]) for i in (1, 2, 3, 4)]
fig.legend([hs[0], hs[2], hs[1], hs[3]], [f"Pop{i} (n = {cnt[f'K4_Pop{i}']})" for i in (1, 3, 2, 4)], ncol=2, frameon=False,
           loc="lower left", bbox_to_anchor=(L - 0.01, 0.0), handlelength=0.9, handleheight=0.8, columnspacing=1.2, fontsize=5.5,
           handletextpad=0.4, labelspacing=0.2)
fig.savefig(OUT / "Fig4a_snp_repair.pdf"); plt.close(fig)

# ---------------------------------------------------------------- b
pc = pd.read_csv(D / "PCA_10_oriented.eigenvec", sep=r"\s+", header=None).set_index(1)
ev = [float(x) for x in open(D / "PCA_10.eigenval")][:2]; tr = float(open(D / "grm_trace.txt").read())   # % of GRM trace
fig = plt.figure(figsize=(W["b"][0] * PT, W["b"][1] * PT)); ax = fig.add_axes([0.2, 0.2, 0.77, 0.77])
X = pc.loc[A.index, 2].to_numpy(); Y = pc.loc[A.index, 3].to_numpy(); mq = A.max_Q.to_numpy(); g = A.dominant_group.to_numpy()
ax.axhline(0, color="#DDDDDD", lw=0.4, zorder=0); ax.axvline(0, color="#DDDDDD", lw=0.4, zorder=0)
lo = mq < 0.7
for gg in ["K4_Pop4", "K4_Pop3", "K4_Pop1", "K4_Pop2"]:
    m = (g == gg)
    ax.scatter(X[m & lo], Y[m & lo], s=4, color=K4C[gg], alpha=0.35, lw=0, zorder=2)
    ax.scatter(X[m & ~lo], Y[m & ~lo], s=5, color=K4C[gg], edgecolor="white", lw=0.2, zorder=3)
# Pop4 ellipse (2 s.d., core accessions)
m = (g == "K4_Pop4") & ~lo; P = np.c_[X[m], Y[m]]; mu0 = np.median(P, 0); dd = np.hypot(*(P - mu0).T); P = P[dd <= np.quantile(dd, 0.86)]
mu = P.mean(0); C = np.cov(P.T); w, v = np.linalg.eigh(C)
ang = np.degrees(np.arctan2(v[1, 1], v[0, 1])); ax.add_patch(Ellipse(mu, 2 * 2.2 * np.sqrt(w[1]), 2 * 2.2 * np.sqrt(w[0]), angle=ang,
                                                                  facecolor=K4C["K4_Pop4"], alpha=0.12, edgecolor=K4C["K4_Pop4"], lw=0.6, zorder=1))
ax.set_xlabel(f"PC1 ({100 * ev[0] / tr:.2f}%)", labelpad=1.5); ax.set_ylabel(f"PC2 ({100 * ev[1] / tr:.2f}%)", labelpad=1.5)
for sp in ("top", "right"): ax.spines[sp].set_visible(False)
hs = [plt.Line2D([], [], marker="o", ls="", ms=2.6, mfc=K4C[f"K4_Pop{i}"], mec="none") for i in (1, 2, 3, 4)]
hs.append(plt.Line2D([], [], marker="o", ls="", ms=2.6, mfc="#CCCCCC", mec="none"))
ax.legend(hs, [f"Pop{i} (n = {cnt[f'K4_Pop{i}']})" for i in (1, 2, 3, 4)] + ["max Q < 0.7"], frameon=False, loc="upper right",
          fontsize=5.5, handletextpad=0.1, labelspacing=0.15, borderaxespad=0.1)
fig.savefig(OUT / "Fig4b_snp_repair.pdf"); plt.close(fig)
print("grey n", int(lo.sum()))

# ---------------------------------------------------------------- c
pf = pd.read_csv(D / "pi_fst_genomewide.tsv", sep="\t"); pf = pf[pf.snp_set == "present"]
pi = {g[3:]: 1000 * v for g, v in pf[(pf.metric == "pi") & pf.estimator.str.startswith("missing-aware")][["group_or_pair", "value"]].values}
F = {(k[3:7], k[11:]): v for k, v in pf[pf.metric == "fst"][["group_or_pair", "value"]].values}
print("pi", pi, "F", F)
pos = {"Pop1": (0.25, 0.64), "Pop4": (0.95, 0.64), "Pop2": (0.25, 0.10), "Pop3": (0.95, 0.10)}
cmap = mpl.colors.LinearSegmentedColormap.from_list("f", ["#F8D4C8", "#F08A78", "#D64A58", "#7F3348"]); norm = mpl.colors.Normalize(0.04, 0.13)
fig = plt.figure(figsize=(W["c"][0] * PT, W["c"][1] * PT)); ax = fig.add_axes([0.0, 0.0, 0.82, 1.0]); ax.axis("off")
for (a, b), f in F.items():
    ax.plot(*zip(pos[a], pos[b]), color=cmap(norm(f)), lw=0.8 + 3.2 * (f - 0.04) / 0.09, zorder=1, solid_capstyle="butt")
for p_, (xx, yy) in pos.items():
    r = 0.08 + 0.06 * (pi[p_] - 4.5) / 1.5
    ax.add_patch(mpl.patches.Circle((xx, yy), r, color=K4C["K4_" + p_], zorder=2, transform=ax.transData))
    ax.text(xx, yy, p_, ha="center", va="center", fontsize=6, zorder=3)
    ax.text(xx, yy + (0.19 if yy > 0.4 else -0.19), f"π = {pi[p_]:.2f} × 10$^{{-3}}$", ha="center", va="center", fontsize=5.5)
ax.set_aspect("equal", adjustable="datalim"); ax.set_xlim(0.0, 1.2); ax.set_ylim(-0.12, 0.86)
cax = fig.add_axes([0.86, 0.2, 0.035, 0.5]); cb = mpl.colorbar.ColorbarBase(cax, cmap=cmap, norm=norm, ticks=[0.05, 0.07, 0.09, 0.11, 0.13])
cb.outline.set_linewidth(0.4); cax.tick_params(labelsize=5.5, width=0.4, length=2); cax.set_title("Mean\n$F_{ST}$", fontsize=5.5, pad=3)
fig.savefig(OUT / "Fig4c_snp_repair.pdf"); plt.close(fig)
