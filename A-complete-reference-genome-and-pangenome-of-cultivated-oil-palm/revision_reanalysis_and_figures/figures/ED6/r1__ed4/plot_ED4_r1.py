#!/usr/bin/env python3
"""[polish copy: tip labels EG-33/EG-8/Houke -> EG_033/EG_008/EG_houke; no other change]
Extended Data Fig. 4 redraw (population diversity, K selection and phylogeny).

Inputs (all tabular):
  sd.pkl                      Source Data sheets ED5a_LD_curve, ED5a_LD50, ED5b_pi, ED5b_fst, ED5c_tree_K4
  src/SpeciesTree_rooted.txt  OrthoFinder rooted species tree (RUN-PAN39-OF315-20260811-001), unchanged rooting
  src/PAN39_tree_tip_K4_mapping.tsv  tip -> display name (identical to ED5c_tree_K4)
  Supplementary Table 16 rows 312-321 (ADMIXTURE CV, K = 2-8, one run per K) -> data/ED4e_ADMIXTURE_CV.tsv
  data/ED4b_pi.tsv, data/ED4b_fst.tsv  panel b inputs (replaceable; columns group/pi/n_windows and
                                       pop1/pop2/FST_mean/n_windows). Initial values = Source Data ED4b_pi/ED4b_fst.

2026-09-24 (ed_fix): panel b reads the two TSVs above; F_ST panel title "Pairwise mean F_ST".
"""
import pickle
import re
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import openpyxl

# polish copy (2026-09-24): inputs read from the beautify work dir; outputs written next to this script
HERE = Path("${WORK_DIR}/fix/beautify/work/ED4")
OUT = Path(__file__).resolve().parent / "out"  # revision r1 copy
DATA = HERE / "data"
DATA_OUT = OUT
OUT.mkdir(exist_ok=True)
MM = 1 / 25.4
DARK = "#242A30"
import sys; sys.path.insert(0, "${WORK_DIR}/fix/beautify/common"); import palA
POP_COL = {f"K4_Pop{i + 1}": c for i, c in enumerate(palA.POP)}  # beautify: scheme A (same HEX as Fig. 4)
POPS = list(POP_COL)

mpl.rcParams.update({
    "font.family": "Arial", "font.size": 7, "axes.labelsize": 7, "axes.titlesize": 7,
    "xtick.labelsize": 6, "ytick.labelsize": 6, "legend.fontsize": 6,
    "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6,
    "xtick.major.size": 2.5, "ytick.major.size": 2.5, "axes.unicode_minus": True,
    "mathtext.fontset": "custom", "mathtext.rm": "Arial", "mathtext.it": "Arial:italic",
    "mathtext.bf": "Arial:bold", "pdf.fonttype": 42, "svg.fonttype": "none",
})


def clean(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def letter(fig, ax, s, dx=-0.05, dy=0.01):
    bb = ax.get_position()
    fig.text(bb.x0 + dx, bb.y1 + dy, s, fontsize=palA.LETTER_PT, fontweight="bold", va="bottom", ha="left")


# ------------------------------------------------------------------ data
d = pickle.load(open(HERE / "sd.pkl", "rb"))
ld, ld50 = d["ED5a_LD_curve"], d["ED5a_LD50"]
pi = pd.read_csv(DATA / "ED4b_pi.tsv", sep="\t").rename(columns={"pi_mean": "pi"})
fst = pd.read_csv(DATA / "ED4b_fst.tsv", sep="\t")
assert set(pi.group) == set(POP_COL) and len(fst) == 6
q = d["ED5c_tree_K4"].set_index("tree_tip")
# [revision r1, 2026-09-26] assembly-to-accession mapping corrected (TK = Y151, NS = Y271; no resequenced accession
# corresponds to the sequenced TN individual): K = 4 proportions from Supplementary Data 16
q = q.copy()
for t in ("dura_hap1", "dura_hap2"):
    q.loc[t, ["Q1", "Q2", "Q3", "Q4"]] = [0.99997, 1e-05, 1e-05, 1e-05]; q.loc[t, "dominant_group"] = "K4_Pop1"
for t in ("pisifera_hap1", "pisifera_hap2"):
    q.loc[t, ["Q1", "Q2", "Q3", "Q4"]] = [1e-05, 1e-05, 0.99997, 1e-05]; q.loc[t, "dominant_group"] = "K4_Pop3"
for t in ("bk_hap1", "bk_hap2"):
    q.loc[t, ["Q1", "Q2", "Q3", "Q4"]] = [np.nan] * 4; q.loc[t, "dominant_group"] = "no_corresponding_accession"
    q.loc[t, "k4_status"] = "no_corresponding_accession"

SCR = Path("${WORK_DIR}")
wb = openpyxl.load_workbook(SCR / "Supplementary_Tables_submission_114RNAseq.xlsx", read_only=True)
rows = list(wb["Supplementary Table 16"].iter_rows(values_only=True))
hdr_i = next(i for i, r in enumerate(rows) if r[0] == "K" and r[1] == "10-fold CV error")
cv = pd.DataFrame([r[:6] for r in rows[hdr_i + 1:hdr_i + 8]], columns=rows[hdr_i][:6])
cv["K"] = cv["K"].astype(int)
cv["10-fold CV error"] = cv["10-fold CV error"].astype(float)
cv.to_csv(DATA_OUT / "ED4e_ADMIXTURE_CV.tsv", sep="\t", index=False)
assert list(cv.K) == [2, 3, 4, 5, 6, 7, 8] and set(cv["Run record"]) == {"1 run"}


# newick parser (no external dependency)
def parse_newick(s):
    s = s.strip().rstrip(";")
    pos = 0

    def node():
        nonlocal pos
        n = {"children": [], "name": "", "len": 0.0, "support": None}
        if s[pos] == "(":
            pos += 1
            while True:
                n["children"].append(node())
                if s[pos] == ",":
                    pos += 1
                    continue
                assert s[pos] == ")"
                pos += 1
                break
        m = re.match(r"[^:,();]*", s[pos:])
        label = m.group(0)
        pos += len(label)
        if n["children"]:
            n["support"] = float(label) if label else None
        else:
            n["name"] = label
        if pos < len(s) and s[pos] == ":":
            pos += 1
            m = re.match(r"[0-9.eE+-]+", s[pos:])
            n["len"] = float(m.group(0))
            pos += len(m.group(0))
        return n

    return node()


tree = parse_newick((HERE / "src/SpeciesTree_rooted.txt").read_text())
names = pd.read_csv(HERE / "src/PAN39_tree_tip_K4_mapping.tsv", sep="\t").set_index("tree_tip")
tips = []


def layout(n, depth=0.0):
    n["x"] = depth + n["len"]
    if not n["children"]:
        n["y"] = len(tips)
        tips.append(n["name"])
    else:
        for c in n["children"]:
            layout(c, n["x"])
        n["y"] = np.mean([c["y"] for c in n["children"]])


tree["len"] = 0.0
layout(tree)
assert len(tips) == 39
ntip = len(tips)

# ------------------------------------------------------------------ figure
fig = plt.figure(figsize=(180 * MM, 168 * MM))
# top row: a LD decay | b pi + FST | e CV
ax_a = fig.add_axes([0.07, 0.715, 0.30, 0.245])
ax_b1 = fig.add_axes([0.45, 0.715, 0.10, 0.245])
ax_b2 = fig.add_axes([0.60, 0.735, 0.14, 0.205])
ax_e = fig.add_axes([0.83, 0.715, 0.155, 0.245])
# bottom: c tree | d ADMIXTURE
ax_c = fig.add_axes([0.03, 0.07, 0.52, 0.56])
ax_d = fig.add_axes([0.70, 0.07, 0.285, 0.56])

# a ---------------------------------------------------------------
for g in POPS:
    s = ld[ld.group_id == g]
    r = ld50.set_index("group_id").loc[g]
    lab = f"Pop{g[-1]} (n = {int(r.sample_n)}; LD50 {'> 500' if str(r.LD50_kb_smooth).startswith('>') else r.LD50_kb_smooth} kb)"
    ax_a.plot(s.dist_kb, s.r2_smooth, color=POP_COL[g], lw=1.1, label=lab)
ax_a.set_xlim(0, 500)
ax_a.set_ylim(0, 0.45)
ax_a.set_xlabel("Distance (kb)")
ax_a.set_ylabel(r"Mean $r^2$")
ax_a.legend(frameon=False, loc="upper right", handlelength=1.4, borderaxespad=0.1)
clean(ax_a)

# b ---------------------------------------------------------------
x = np.arange(4)
ax_b1.bar(x, pi.set_index("group").loc[POPS].pi * 1e3, color=[POP_COL[g] for g in POPS], width=0.7)
ax_b1.set_xticks(x, [f"Pop{i}" for i in range(1, 5)], rotation=45, ha="right")
ax_b1.set_ylabel(r"Nucleotide diversity $\pi$ ($\times10^{-3}$)")
ax_b1.set_ylim(0, 7)
clean(ax_b1)
M = np.full((4, 4), np.nan)
for _, r in fst.iterrows():
    i, j = POPS.index(r.pop1), POPS.index(r.pop2)
    M[i, j] = M[j, i] = r.FST_mean
im = ax_b2.imshow(M, cmap="Blues", vmin=0, vmax=0.14)
for i in range(4):
    for j in range(4):
        if i == j:
            ax_b2.add_patch(plt.Rectangle((j - .5, i - .5), 1, 1, color="#F2F2F2"))
            ax_b2.text(j, i, "–", ha="center", va="center", fontsize=6)
        else:
            ax_b2.text(j, i, f"{M[i, j]:.3f}", ha="center", va="center", fontsize=5.5,
                       color="white" if M[i, j] > 0.09 else DARK)
ax_b2.set_xticks(x, [f"Pop{i}" for i in range(1, 5)], rotation=45, ha="right")
ax_b2.set_yticks(x, [f"Pop{i}" for i in range(1, 5)])
ax_b2.set_title(r"Pairwise mean $F_{\mathrm{ST}}$", pad=3)
for s in ax_b2.spines.values():
    s.set_linewidth(0.4)

# e ---------------------------------------------------------------
ax_e.plot(cv.K, cv["10-fold CV error"], color=DARK, lw=0.9, marker="o", ms=3, zorder=2)
k4 = cv.set_index("K").loc[4, "10-fold CV error"]
kmin = cv.loc[cv["10-fold CV error"].idxmin()]
ax_e.scatter([4], [k4], s=30, facecolor="none", edgecolor=POP_COL["K4_Pop2"], lw=1.0, zorder=3)
ax_e.annotate("K = 4\n(working partition)", xy=(4, k4), xytext=(2.0, 0.2828), fontsize=5.5,
              arrowprops=dict(arrowstyle="-", lw=0.5, color="#555555"))
ax_e.annotate(f"minimum\n(K = {int(kmin.K)})", xy=(kmin.K, kmin["10-fold CV error"]),
              xytext=(7.3, 0.2835), fontsize=5.5, ha="center",
              arrowprops=dict(arrowstyle="-", lw=0.5, color="#555555"))
ax_e.set_xticks(cv.K)
ax_e.set_xlabel("K")
ax_e.set_ylabel("ADMIXTURE CV error")
ax_e.set_ylim(0.282, 0.299)
clean(ax_e)

# c tree ----------------------------------------------------------
xs = []


def draw(n):
    if n["children"]:
        ys = [c["y"] for c in n["children"]]
        ax_c.plot([n["x"], n["x"]], [min(ys), max(ys)], color="#4A4A4A", lw=0.55)
        for c in n["children"]:
            ax_c.plot([n["x"], c["x"]], [c["y"], c["y"]], color="#4A4A4A", lw=0.55)
            draw(c)
        if n["support"] is not None and n["support"] < 0.95 and n is not tree:
            ax_c.text(n["x"] - 0.00005, n["y"] - 0.35, f"{n['support']:.2f}", fontsize=5.0, color="#6E6E6E",
                      ha="right", va="bottom", zorder=5,
                      bbox=dict(fc="white", ec="none", pad=0.3, alpha=0.9))  # beautify: keep label legible over branches
    else:
        xs.append(n["x"])


draw(tree)
tipx = max(xs)
label_x = tipx + 0.0006
def relabel(lab):
    """polish: accession labels as in the main text and Figs 4-5 (EG_033, EG_008, EG_houke)"""
    m = re.fullmatch(r"EG-(\d+)", lab)
    if m:
        return f"EG_{int(m.group(1)):03d}"
    return "EG_houke" if lab == "Houke" else lab


for n_name, yv in zip(tips, range(ntip)):
    lab = relabel(names.loc[n_name, "display_name"])
    if lab.startswith("E. oleifera"):
        ax_c.text(label_x, yv, r"$\it{E.\ oleifera}$" + lab.replace("E. oleifera", ""), fontsize=5.5, va="center")
    else:
        ax_c.text(label_x, yv, lab, fontsize=5.5, va="center")
# dotted leaders from each tip to the label column
def leaders(n):
    if n["children"]:
        for c in n["children"]:
            leaders(c)
    else:
        ax_c.plot([n["x"], tipx + 0.0004], [n["y"], n["y"]], color="#D0D0D0", lw=0.35, ls=(0, (1, 1.5)))


leaders(tree)
ax_c.plot([0.0005, 0.0055], [ntip + 0.2, ntip + 0.2], color="black", lw=0.8)
ax_c.text(0.003, ntip + 0.6, "0.005 substitutions per site", ha="center", va="top", fontsize=5.5)
ax_c.set_ylim(ntip + 1.3, -1)
ax_c.set_xlim(-0.0003, label_x + 0.0075)
ax_c.axis("off")

# d ADMIXTURE bars aligned to tips ---------------------------------
for yv, t in enumerate(tips):
    r = q.loc[t]
    if r.k4_status != "reviewed_material_mapping":
        ax_d.add_patch(plt.Rectangle((0, yv - 0.38), 1, 0.76, facecolor="#F2F2F2", edgecolor="#9A9A9A", lw=0.4))
        ax_d.text(0.5, yv, "not in the 308-accession panel", ha="center", va="center", fontsize=5.0, color="#555555")
        continue
    left = 0
    for i, g in enumerate(POPS, start=1):
        w = r[f"Q{i}"]
        ax_d.barh(yv, w, left=left, height=0.76, color=POP_COL[g], lw=0)
        left += w
ax_d.set_ylim(ntip + 1.3, -1)
ax_d.set_xlim(0, 1)
ax_d.set_yticks([])
ax_d.set_xticks([0, 0.5, 1], ["0", "0.5", "1"])
ax_d.set_xlabel("Ancestry proportion (K = 4)")
for s in ("top", "right", "left"):
    ax_d.spines[s].set_visible(False)
handles = [plt.Rectangle((0, 0), 1, 1, color=POP_COL[g]) for g in POPS]
ax_d.legend(handles, [f"Pop{i}" for i in range(1, 5)], ncol=4, frameon=False, loc="lower center",
            bbox_to_anchor=(0.5, 1.0), handlelength=1.0, columnspacing=0.8, handletextpad=0.3)

for ax, s, dx in [(ax_a, "a", -0.055), (ax_b1, "b", -0.055), (ax_e, "e", -0.06), (ax_c, "c", -0.0),
                  (ax_d, "d", -0.02)]:
    letter(fig, ax, s, dx=dx, dy=0.012 if s in "abe" else 0.02)

for ext in ("pdf", "png"):
    fig.savefig(OUT / f"Extended_Data_Fig_04_candidate.{ext}", dpi=600, facecolor="white")

# source tables for redrawn/new panels
cv.to_csv(DATA_OUT / "ED4e_ADMIXTURE_CV.tsv", sep="\t", index=False)
pd.DataFrame({"tip_order_top_to_bottom": range(1, ntip + 1), "tree_tip": tips,
              "display_name": [relabel(names.loc[t, "display_name"]) for t in tips]}).merge(
    q.reset_index()[["tree_tip", "k4_status", "Q1", "Q2", "Q3", "Q4", "dominant_group"]], on="tree_tip").to_csv(
    DATA_OUT / "ED4cd_tree_tip_order_K4.tsv", sep="\t", index=False)
print("pi x1e3", (pi.pi * 1e3).round(2).tolist(), "FST range", fst.FST_mean.min(), fst.FST_mean.max())
print("CV", cv[["K", "10-fold CV error"]].values.tolist())
print("LD50", ld50[["label", "sample_n", "LD50_kb_smooth"]].values.tolist())
