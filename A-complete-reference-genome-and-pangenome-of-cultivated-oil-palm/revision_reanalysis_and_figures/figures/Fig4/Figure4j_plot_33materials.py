#!/usr/bin/env python3
"""Fig. 4j redrawn with the 33-material RGA data (author-corrected source, 2026-09-24), in the style and slot of
the current Figure 4 panel j.  Data semantics follow build_rga_tree_dotmatrix_redesign.py / 06_plot_rga_33.py:
dot size = raw RGA count, dot colour = column z-score (clipped to ±3), red dots = private RGA copies
(labelled when <= 12 or >= 100), cladogram with aligned tips (branch lengths not used).

Output: panel_j_33.pdf, exactly 224.5 x 273.8 pt, to be placed at (0, 347) on the Figure 4 page."""
import csv, re
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.patches import Rectangle

H = Path(__file__).resolve().parent
plt.rcParams.update({"font.family": "Arial", "pdf.fonttype": 42, "axes.unicode_minus": True})
W_PT, H_PT, Y0 = 224.5, 273.8, 347.0            # slot on the Figure 4 page (x 0-224.5, y 347-620.8)

# ---------------------------------------------------------------- data
tips = list(csv.DictReader(open(H / "RGA_tree_dotmatrix_tip_order.tsv"), delimiter="\t"))
tips.sort(key=lambda r: int(r["Y_order_top_to_bottom"]))
DISPLAY = {"dura": "TK", "pisifera": "NS", "BK": "TN", "FL": "FL", "meizhou4": "E. oleifera",
           "niriliya": "Nigerian", "houke": "EG_houke"}
GROUP = {"FL": "FL", "BK": "TN", "dura": "par", "pisifera": "par", "niriliya": "par", "houke": "par",
         "meizhou4": "ole"}
GCOL = {"FL": "#62B6B7", "TN": "#9C8AC7", "par": "#D89C35", "EG": "#86B68A", "ole": "#B8B8B8"}
disp = lambda g: DISPLAY.get(g, f"EG_{int(g):03d}" if g.isdigit() else g)

z = {}
for r in csv.DictReader(open(H / "RGA_major_group_counts_and_zscores.tsv"), delimiter="\t"):
    z[(r["Genome"], r["RGA_type"])] = (float(r["Count"]), float(r["Zscore_clipped"]))
priv = {r["Raw_genome"]: int(r["Private_RGA_copies"])
        for r in csv.DictReader(open(H / "RGA_global_private_synteny_block_audit.tsv"), delimiter="\t")}
COLS = ["NBS-class", "RLK", "RLP", "TM-CC", "Other"]
LAB = ["NBS", "RLK", "RLP", "TM-CC", "Other"]


# ---------------------------------------------------------------- tree (newick, ladderized like Bio.Phylo reverse=True)
def parse(s):
    s = re.sub(r"\)[0-9.]+:", "):", s.strip().rstrip(";"))
    pos = 0

    def node():
        nonlocal pos
        if s[pos] == "(":
            pos += 1
            kids = [node()]
            while s[pos] == ",":
                pos += 1
                kids.append(node())
            pos += 1                                  # ')'
            m = re.match(r"[^,():;]*(:[0-9.eE-]+)?", s[pos:]); pos += m.end()
            return {"kids": kids}
        m = re.match(r"([^,():;]+)(:[0-9.eE-]+)?", s[pos:]); pos += m.end()
        return {"name": m.group(1), "kids": []}
    return node()


def nt(n):
    return 1 if not n["kids"] else sum(nt(k) for k in n["kids"])


def ladder(n):
    for k in n["kids"]:
        ladder(k)
    n["kids"].sort(key=nt, reverse=True)


tree = parse(open(H / "SpeciesTree_33varieties_display.nwk").read())
ladder(tree)
term = []
def collect(n):
    if not n["kids"]:
        term.append(n)
    for k in n["kids"]:
        collect(k)
collect(tree)
order = [t["name"] for t in term]
assert order == [r["Tree_tip"] for r in tips], (order, [r["Tree_tip"] for r in tips])
N = len(term)

# ---------------------------------------------------------------- page geometry (pt, page coordinates)
TOP, BOT = 381.0, 559.0                          # first / last row centres (as in the current panel)
ROW = (BOT - TOP) / (N - 1)
ry = {t["name"]: TOP + i * ROW for i, t in enumerate(term)}
TX0, TX1 = 8.0, 61.5                              # tree
STRIP = (63.8, 67.0)
LABEL_X = 112.8
CX = dict(zip(COLS, [118.7, 137.2, 155.4, 174.0, 192.4]))
PX, PNUM_X = 211.3, 215.2

fig = plt.figure(figsize=(W_PT / 72, H_PT / 72))
ax = fig.add_axes([0, 0, 1, 1])
ax.set_xlim(0, W_PT); ax.set_ylim(Y0 + H_PT, Y0)   # page coordinates, y downwards
ax.axis("off")

# tree
depth = {}
def setx(n, d=0):
    depth[id(n)] = d
    for k in n["kids"]:
        setx(k, d + 1)
setx(tree)
maxd = max(depth[id(t)] for t in term)
xof = lambda n: TX1 if not n["kids"] else TX0 + (TX1 - TX0) * depth[id(n)] / maxd
def yof(n):
    return ry[n["name"]] if not n["kids"] else sum(yof(k) for k in n["kids"]) / len(n["kids"])
def draw(n):
    if not n["kids"]:
        return
    x, ys = xof(n), [yof(k) for k in n["kids"]]
    ax.plot([x, x], [min(ys), max(ys)], color="#3F3F3F", lw=0.5, solid_capstyle="butt")
    for k in n["kids"]:
        ax.plot([x, xof(k)], [yof(k), yof(k)], color="#3F3F3F", lw=0.5, solid_capstyle="butt")
        draw(k)
draw(tree)

cmap = LinearSegmentedColormap.from_list("z", ["#4B5AA7", "#A7AED3", "#F7F7F7", "#E8A39E", "#C84B4B"])
norm = Normalize(-3, 3)
maxc = max(v[0] for v in z.values())
dsz = lambda c: 0 if c <= 0 else (10.0 + (c / maxc) ** 0.72 * 250.0) * 0.19
maxp = max(priv.values())
psz = lambda c: 0 if c <= 0 else (9.0 + (c / maxp) ** 0.68 * 150.0) * 0.22

for r in tips:
    g, y = r["Genome"], ry[r["Tree_tip"]]
    ax.add_patch(Rectangle((STRIP[0], y - ROW / 2), STRIP[1] - STRIP[0], ROW, color=GCOL[GROUP.get(g, "EG")], lw=0))
    name = disp(g)
    ax.text(LABEL_X, y, name, ha="right", va="center", fontsize=6, color="#222222",
            style="italic" if name == "E. oleifera" else "normal")
    for c in COLS:
        cnt, zz = z[(g, c)]
        if cnt > 0:
            ax.scatter(CX[c], y, s=dsz(cnt), color=cmap(norm(zz)), edgecolors="white", linewidths=0.25, zorder=3)
    pc = priv.get(g, 0)
    if pc > 0:
        ax.scatter(PX, y, s=psz(pc), color="#D86C59", edgecolors="white", linewidths=0.3, alpha=0.9, zorder=3)
        if pc <= 12 or pc >= 100:
            ax.text(PNUM_X, y, str(pc), ha="left", va="center", fontsize=5.5, color="#9D3E34")

for c, lab in zip(COLS, LAB):
    ax.text(CX[c], 374.2, lab, ha="center", va="baseline", fontsize=6, color="#222222")
ax.text(PX, 371.3, "Private", ha="center", va="baseline", fontsize=6, color="#222222")
ax.text(PX, 377.1, "RGA", ha="center", va="baseline", fontsize=6, color="#222222")
ax.text(7.9, 568.0, "Cladogram; branch lengths not used.", fontsize=5.5, color="#555555", va="baseline")

# legends
ax.text(7.9, 577.6, "Genome group", fontsize=6, fontweight="bold", va="baseline")
entries = [("FL", "FL"), ("TN", "TN"), ("par", "TK, NS, Nigerian and EG_houke"), ("EG", "EG assemblies"),
           ("ole", "E. oleifera")]
for i, (k, lab) in enumerate(entries):
    y = 584.2 + i * 7.0
    ax.add_patch(Rectangle((7.9, y - 2.3), 4.6, 4.6, color=GCOL[k], lw=0))
    ax.text(15.0, y, lab, fontsize=5.5, va="center", style="italic" if k == "ole" else "normal")

ax.text(107.4, 577.6, "RGA count", fontsize=6, fontweight="bold", va="baseline")
for v, x in zip([50, 150, 300, 600], [112.7, 139.0, 164.0, 191.0]):
    ax.scatter(x, 582.0, s=dsz(v), facecolor="white", edgecolor="#333333", linewidths=0.4)
    ax.text(x, 589.2, str(v), ha="center", va="baseline", fontsize=5.5)
ax.text(106.9, 598.2, "Private RGA copies", fontsize=6, fontweight="bold", va="baseline")
for v, x in zip([10, 50, 100], [166.0, 183.5, 201.5]):
    ax.scatter(x, 596.2, s=psz(v), color="#D86C59", edgecolors="white", linewidths=0.3)
    ax.text(x + 3.0, 596.2, str(v), fontsize=5.5, va="center")
ax.text(106.8, 606.3, "Normalized abundance (column z-score)", fontsize=6, fontweight="bold", va="baseline")
cb = fig.add_axes([106.8 / W_PT, 1 - (612.0 - Y0) / H_PT, (209.5 - 106.8) / W_PT, 3.0 / H_PT])
import numpy as np
cb.imshow(np.linspace(-3, 3, 256)[None, :], cmap=cmap, norm=norm, aspect="auto", extent=(-3, 3, 0, 1))
cb.set_yticks([]); cb.set_xticks([])
for s in cb.spines.values():
    s.set_linewidth(0.3)
for v in (-3, -1.5, 0, 1.5, 3):
    x = 106.8 + (v + 3) / 6 * (209.5 - 106.8)
    ax.text(x, 618.6, ("−" if v < 0 else "") + (f"{abs(v):g}"), ha="center", va="baseline", fontsize=5.5)

fig.savefig(H / "panel_j_33.pdf", transparent=True)
# plotting table for Source Data
with open(H / "Fig4j_33_materials.tsv", "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["Row_top_to_bottom", "Material", "Display_name", "Genome_group"] +
               [f"{c}_count" for c in LAB] + [f"{c}_zscore_clipped" for c in LAB] + ["Private_RGA_copies"])
    for i, r in enumerate(tips, 1):
        g = r["Genome"]
        grp = {"FL": "FL", "TN": "TN", "par": "TK, NS, Nigerian and EG_houke", "EG": "EG assemblies",
               "ole": "E. oleifera"}[GROUP.get(g, "EG")]
        w.writerow([i, g, disp(g), grp] + [int(z[(g, c)][0]) for c in COLS] +
                   [round(z[(g, c)][1], 4) for c in COLS] + [priv.get(g, 0)])
print("rows", N)
