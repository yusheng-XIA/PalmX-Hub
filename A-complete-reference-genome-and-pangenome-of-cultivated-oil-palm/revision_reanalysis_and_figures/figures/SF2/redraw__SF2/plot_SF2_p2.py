#!/usr/bin/env python3
"""Supplementary Fig. 2, page 2 (panels c-f): maximum-likelihood lipid-enzyme gene trees.

Input: IQ-TREE ML trees (Fig2e_<family>_phylogeny_ml.treefile; node labels are
SH-aLRT/UFBoot, 1,000 ultrafast-bootstrap replicates) and tip metadata from
03_V3/02_figure/05_Figure2_evolution_multiomics_panels_flat_20260727/.
IQ-TREE trees are unrooted; they are drawn here midpoint-rooted (the previous
figure drew them from the first taxon, an arbitrary position). Supports are
attached to bipartitions, so rerooting does not change them.
"""
import re
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D

HERE = Path(__file__).resolve().parent
SRC = HERE / "src"
OUT = HERE / "out"
OUT.mkdir(exist_ok=True)
MM = 1 / 25.4
DARK = "#242A30"

mpl.rcParams.update({
    "font.family": "Arial", "font.size": 6, "axes.labelsize": 6, "xtick.labelsize": 5.5,
    "axes.linewidth": 0.5, "xtick.major.width": 0.5, "xtick.major.size": 2,
    "axes.unicode_minus": True, "mathtext.fontset": "custom", "mathtext.rm": "Arial",
    "mathtext.it": "Arial:italic", "pdf.fonttype": 42, "svg.fonttype": "none",
})

GENOME_COL = {
    "American_hap1": "#D55E00", "Dura": "#009E73", "Pisifera": "#009E73",
    "Cocos_nucifera": "#3C78A8", "Areca_catechu": "#3C78A8", "Phoenix_dactylifera": "#3C78A8",
    "Nypa_fruticans": "#3C78A8", "Calamus": "#3C78A8", "Daemonorops": "#3C78A8",
    "Arabidopsis_thaliana": "#7F7F7F",
}
TREES = [("c", "FAD", "FAD desaturases (FAD2/6 and FAD3/7/8)"),
         ("d", "FAT", "FATA/FATB acyl-ACP thioesterases"),
         ("e", "DGAT1", "DGAT1"),
         ("f", "DGAT2", "DGAT2")]
TIP_FS = 4.8


# ---------------------------------------------------------------- Newick handling
def parse_newick(s):
    """Return (adjacency {node: {nbr: length}}, tip names {node: name}, support {frozenset(edge): (sh, uf)})."""
    s = s.strip().rstrip(";")
    adj, names, sup = {}, {}, {}
    counter = [0]

    def new():
        counter[0] += 1
        adj[counter[0]] = {}
        return counter[0]

    pos = 0

    def parse_label():
        nonlocal pos
        m = re.match(r"[^,():;]*", s[pos:])
        pos += m.end()
        return m.group(0)

    def parse_len():
        nonlocal pos
        if pos < len(s) and s[pos] == ":":
            pos += 1
            m = re.match(r"[-\d.eE]+", s[pos:])
            pos += m.end()
            return float(m.group(0))
        return 0.0

    def parse_sub():
        nonlocal pos
        node = new()
        if s[pos] == "(":
            pos += 1
            while True:
                child, lab, ln = parse_sub()
                adj[node][child] = ln
                adj[child][node] = ln
                if lab:
                    sh, uf = (lab.split("/") + [None])[:2]
                    sup[frozenset((node, child))] = (float(sh), float(uf) if uf else None)
                if s[pos] == ",":
                    pos += 1
                    continue
                assert s[pos] == ")"
                pos += 1
                break
            lab = parse_label()
            ln = parse_len()
            return node, lab, ln
        lab = parse_label()
        names[node] = lab
        ln = parse_len()
        return node, None, ln

    root, _, _ = parse_sub()
    return adj, names, sup, root


def farthest(adj, start):
    dist = {start: 0.0}
    stack = [start]
    prev = {start: None}
    while stack:
        u = stack.pop()
        for v, w in adj[u].items():
            if v not in dist:
                dist[v] = dist[u] + w
                prev[v] = u
                stack.append(v)
    far = max(dist, key=dist.get)
    return far, dist, prev


def midpoint_root(adj, names, sup):
    tip0 = next(iter(names))
    a, _, _ = farthest(adj, tip0)
    b, dist, prev = farthest(adj, a)
    half = dist[b] / 2
    # walk from b towards a until the cumulative distance passes the midpoint
    path = [b]
    while prev[path[-1]] is not None:
        path.append(prev[path[-1]])
    acc = 0.0
    for u, v in zip(path, path[1:]):
        w = adj[u][v]
        if acc + w >= half:
            x = half - acc          # distance from u along the edge u-v
            r = max(adj) + 1
            adj[r] = {u: x, v: w - x}
            del adj[u][v], adj[v][u]
            adj[u][r] = x
            adj[v][r] = w - x
            key = frozenset((u, v))
            if key in sup:          # the bipartition is unchanged; keep support on both halves
                sup[frozenset((r, u))] = sup[key]
                sup[frozenset((r, v))] = sup[key]
                del sup[key]
            return r, dist[b]
        acc += w
    raise RuntimeError("midpoint not found")


def to_rooted(adj, root):
    children = {}
    blen = {root: 0.0}
    order = [root]
    seen = {root}
    i = 0
    while i < len(order):
        u = order[i]
        children[u] = [v for v in adj[u] if v not in seen]
        for v in children[u]:
            seen.add(v)
            blen[v] = adj[u][v]
            order.append(v)
        i += 1
    return children, blen


def ntips(children, u, memo={}):
    if not children[u]:
        return 1
    return sum(ntips(children, c) for c in children[u])


# ---------------------------------------------------------------- drawing
def draw_tree(ax, family, title, letter):
    raw = (SRC / f"Fig2e_{family}_phylogeny_ml.treefile").read_text()
    meta = pd.read_csv(SRC / f"Fig2e_{family}_phylogeny_tip_metadata.tsv", sep="\t").set_index("Tip")
    adj, names, sup, _ = parse_newick(raw)
    assert set(names.values()) == set(meta.index), family
    root, diam = midpoint_root(adj, names, sup)
    children, blen = to_rooted(adj, root)
    # ladderize (small clades first)
    size = {}

    def count(u):
        size[u] = 1 if not children[u] else sum(count(c) for c in children[u])
        return size[u]
    count(root)
    for u in children:
        children[u].sort(key=lambda c: size[c])
    x = {root: 0.0}
    y = {}
    tips = []

    def place(u):
        for c in children[u]:
            x[c] = x[u] + blen[c]
            place(c)
        if not children[u]:
            y[u] = len(tips)
            tips.append(u)
        else:
            y[u] = (y[children[u][0]] + y[children[u][-1]]) / 2
    place(root)
    n = len(tips)
    xmax = max(x.values())
    lw = 0.55
    for u, cs in children.items():
        if cs:
            ax.plot([x[u], x[u]], [y[cs[0]], y[cs[-1]]], color=DARK, lw=lw, solid_capstyle="butt")
            for c in cs:
                ax.plot([x[u], x[c]], [y[c], y[c]], color=DARK, lw=lw, solid_capstyle="butt")
    ax.plot([-xmax * 0.015, 0], [y[root], y[root]], color=DARK, lw=lw)
    # supports (UFBoot) on internal nodes
    n_strong = n_mid = 0
    for u, cs in children.items():
        if not cs or u == root:
            continue
        parent = next(p for p, cc in children.items() if u in cc)
        sh, uf = sup.get(frozenset((parent, u)), (None, None))
        if uf is None:
            continue
        if uf >= 95:
            ax.scatter(x[u], y[u], s=5, color=DARK, zorder=4, lw=0)
            n_strong += 1
        elif uf >= 70:
            ax.scatter(x[u], y[u], s=5, facecolor="white", edgecolor=DARK, lw=0.45, zorder=4)
            n_mid += 1
    for t in tips:
        name = names[t]
        g = meta.loc[name, "Genome"]
        bold = bool(meta.loc[name, "Oil_palm"])
        ax.text(x[t] + xmax * 0.012, y[t], name, va="center", ha="left", fontsize=TIP_FS,
                color=GENOME_COL.get(g, DARK), fontweight="bold" if bold else "normal")
    ax.set_ylim(n - 0.3, -0.7)
    ax.set_xlim(-xmax * 0.02, xmax * 1.28)
    for sp in ("top", "right", "left"):
        ax.spines[sp].set_visible(False)
    ax.set_yticks([])
    ax.set_xlabel("Substitutions per site", labelpad=1.5)
    ax.tick_params(axis="x", pad=1)
    ax.set_title(title, loc="left", fontsize=6.5, pad=2)
    # Newick with supports for Source Data
    def nwk(u):
        if not children[u]:
            return f"{names[u]}:{blen[u]:.6g}"
        inner = ",".join(nwk(c) for c in children[u])
        lab = ""
        if u != root:
            parent = next(p for p, cc in children.items() if u in cc)
            s2 = sup.get(frozenset((parent, u)))
            if s2:
                lab = f"{s2[0]:g}/{s2[1]:g}" if s2[1] is not None else f"{s2[0]:g}"
        return f"({inner}){lab}:{blen[u]:.6g}" if u != root else f"({inner});"
    (HERE / f"SF2{letter}_{family}.nwk").write_text(nwk(root) + "\n")
    return n, n_strong, n_mid, xmax, diam


fig = plt.figure(figsize=(180 * MM, 226 * MM))
# heights proportional to tip counts
counts = {f: len(pd.read_csv(SRC / f"Fig2e_{f}_phylogeny_tip_metadata.tsv", sep="\t")) for _, f, _ in TREES}
top_n = max(counts["FAD"], counts["FAT"])
bot_n = max(counts["DGAT1"], counts["DGAT2"])
usable = 0.80
h_top = usable * top_n / (top_n + bot_n)
h_bot = usable * bot_n / (top_n + bot_n)
y_top0 = 0.965 - h_top
y_bot0 = y_top0 - 0.062 - h_bot
boxes = {"FAD": [0.06, y_top0, 0.40, h_top], "FAT": [0.555, y_top0, 0.40, h_top],
         "DGAT1": [0.06, y_bot0, 0.40, h_bot], "DGAT2": [0.555, y_bot0, 0.40, h_bot]}
summary = []
for letter, fam, title in TREES:
    ax = fig.add_axes(boxes[fam])
    n, s95, s70, xmax, diam = draw_tree(ax, fam, title, letter)
    fig.text(boxes[fam][0] - 0.045, boxes[fam][1] + boxes[fam][3] + 0.012, letter, fontsize=9,
             fontweight="bold", va="bottom")
    summary.append((letter, fam, n, s95, s70, round(xmax, 4)))

handles = [
    Line2D([], [], color="#D55E00", lw=0, marker="s", ms=4, label="FL-Hap1 (Ame)"),
    Line2D([], [], color="#009E73", lw=0, marker="s", ms=4, label="TK, dura (Dur); NS, pisifera (Pis)"),
    Line2D([], [], color="#3C78A8", lw=0, marker="s", ms=4,
           label="Other palms (Coc, Are, Pho, Nyp, Cal, Dae)"),
    Line2D([], [], color="#7F7F7F", lw=0, marker="s", ms=4, label=r"$\mathit{Arabidopsis\ thaliana}$ (Ath)"),
    Line2D([], [], color=DARK, lw=0, marker="o", ms=2.8, label="UFBoot ≥ 95"),
    Line2D([], [], markerfacecolor="white", markeredgecolor=DARK, color=DARK, lw=0, marker="o",
           ms=2.8, markeredgewidth=0.5, label="UFBoot 70–94"),
]
fig.legend(handles=handles, loc="lower center", ncol=3, frameon=False, fontsize=5.8,
           bbox_to_anchor=(0.5, 0.005), handletextpad=0.3, columnspacing=1.2)
for ext in ("pdf", "png"):
    fig.savefig(OUT / f"Supplementary_Fig_02_p2.{ext}", dpi=600, facecolor="white")
for r in summary:
    print(r)
