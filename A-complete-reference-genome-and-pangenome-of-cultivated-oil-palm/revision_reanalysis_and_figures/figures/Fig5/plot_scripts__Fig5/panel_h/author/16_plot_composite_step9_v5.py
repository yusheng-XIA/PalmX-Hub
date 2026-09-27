#!/usr/bin/env python3
"""
组合 Fig5a (GWAS 增强): 左=IPHs 马赛克(SV3-blue)+优先 GWAS 位点标记(暖金=抓到有利等位/灰=未)+图例;
右=chr01B 放大(负荷热图+红框选择路径)。
用法: python3 16_plot_composite_step9_v5.py <dSNP_ALT_derived.tsv>
"""
import os, sys
from collections import defaultdict
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.patches import FancyBboxPatch, Patch, Rectangle
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.lines import Line2D
import matplotlib.font_manager as fm
import seaborn as sns
from PIL import Image

D    = "${ANALYSIS_DIR}/21_MS/06_result/dSVs"
DSV  = f"{D}/results/01_core_dsv/dsv_v5_candidates.tsv"
CHRFAI = "${ANALYSIS_DIR}/08_hifi_chromosome/05_fianl_all_genome_chrom_only/Africa_hap2.fasta.fai"
OUT  = "${ANALYSIS_DIR}/21_MS/06_result/dSVs/results/09_ideal_parent_haplotypes"
WINDOW = 500_000
EXCLUDE = {"American_hap1"}; NAME_MAP = {"dura": "EG_dura", "pisifera": "EG_pisifera"}
def norm(s): return NAME_MAP.get(s.strip(), s.strip())
DSNP_PATH = sys.argv[1]; ZC = "chr01B"
OUT_BASENAME = "Fig5a_composite_STEP9_v5"

# Palette harmonized with Fig5_SV3.jpg, v5:
# keep the successful pale/medium blues, remove green as a dominant hue,
# and use warm gold only for favourable GWAS hits.
MOSAIC_COLORS = [
    "#E9F3FA", "#D8E8F8", "#C8D8E8", "#AFCFE3",
    "#8FBAD6", "#68A8D8", "#5A94BF", "#7898C8",
    "#627FB0", "#4D6F9A", "#3D638E", "#2D536F", "#213F58"
]
HEATMAP_COLORS = [
    "#FBFCFE", "#EFF6FB", "#E2EFF7", "#CDE2F1",
    "#AFCFE6", "#88BBDC", "#68A8D8", "#3F83B0", "#245E7A"
]
CHROM_BG = "#F3F7FA"
CAPTURED_COLOR = "#D88838"
MISSED_COLOR = "#B8C5CC"
SELECT_RED = "#C84848"
SELECT_TEXT_RED = "#B83F3E"

for cand in ["Arial", "Helvetica", "DejaVu Sans"]:
    if any(cand == f.name for f in fm.fontManager.ttflist):
        plt.rcParams["font.family"] = cand; break
plt.rcParams.update({"axes.edgecolor": "#9aa3ab", "axes.linewidth": 0.6, "xtick.color": "#4d555c",
    "ytick.color": "#4d555c", "text.color": "#2b3136", "axes.labelcolor": "#2b3136",
    "xtick.major.width": 0.6, "mathtext.fontset": "dejavusans", "mathtext.default": "regular"})

chromlen = {ln.split("\t")[0]: int(ln.split("\t")[1]) for ln in open(CHRFAI)}
chroms = [f"chr{n:02d}B" for n in range(1, 17)]
path = {c: {} for c in chroms}
for ln in open(f"{OUT}/ideal_step9_ALL_path.tsv"):
    if ln.startswith("Chrom"): continue
    ch, wi, ws, don = ln.rstrip("\n").split("\t"); path[ch][int(wi)] = norm(don)
counts = {}
for ln in open(f"{OUT}/ideal_step9_ALL_by_chrom.tsv"):
    if ln.startswith(("Chrom", "TOTAL")): continue
    c = ln.rstrip("\n").split("\t"); counts[c[0]] = (c[1], c[2], c[4], c[5])  # dSV,dSNP,favcap,favtot
# GWAS 位点(去重 chrom:pos, 绿=任一性状抓到)
loci_mark = {}
for ln in open(f"{OUT}/ideal_step9_ALL_loci.tsv"):
    if ln.startswith("SV\t"): continue
    sv, ch, pos, tgt, sel, cap = ln.rstrip("\n").split("\t")
    k = (ch, int(pos)); loci_mark[k] = max(loci_mark.get(k, 0), int(cap))

contrib = {}
for ch in chroms:
    for d in path[ch].values(): contrib[d] = contrib.get(d, 0) + 1
used = sorted(contrib, key=lambda d: -contrib[d])
donor_cmap = LinearSegmentedColormap.from_list("sv3_mosaic_v5_blue", MOSAIC_COLORS)
lo, span = (0.06, 0.88)
shade = {d: donor_cmap(lo + span*(i/max(1, len(used)-1))) for i, d in enumerate(used)}

# chr01B 放大负荷
n_win = chromlen[ZC]//WINDOW + 1
dsv_w = defaultdict(lambda: np.zeros(n_win, int)); dsnp_w = defaultdict(lambda: np.zeros(n_win, int))
with open(DSNP_PATH) as f:
    h = f.readline().rstrip("\n").split("\t"); ci, pi, sm = h.index("Chrom"), h.index("Pos"), h.index("Samples")
    pol = h.index("ALT_Polarity_Phoenix") if "ALT_Polarity_Phoenix" in h else None
    for ln in f:
        c = ln.rstrip("\n").split("\t")
        if c[ci] != ZC: continue
        if pol is not None and c[pol] != "ALT_Derived": continue
        w = min(int(c[pi])//WINDOW, n_win-1)
        for s in c[sm].replace(";", ",").split(","):
            s = norm(s)
            if s and s not in EXCLUDE: dsnp_w[s][w] += 1
with open(DSV) as f:
    h = f.readline().rstrip("\n").split("\t"); ci, si, ei, sm = h.index("Chrom"), h.index("Start"), h.index("End"), h.index("Samples")
    for ln in f:
        c = ln.rstrip("\n").split("\t")
        if c[ci] != ZC: continue
        s = norm(c[sm])
        if s in EXCLUDE: continue
        dsv_w[s][min(((int(c[si])+int(c[ei]))//2)//WINDOW, n_win-1)] += 1
donors = sorted(set(dsnp_w) | set(dsv_w))
load = np.array([dsnp_w[d] + dsv_w[d] for d in donors]); order = np.argsort(load.sum(1))
donors = [donors[i] for i in order]; load = load[order]; row_of = {d: i for i, d in enumerate(donors)}

# ===== 画 =====
fig = plt.figure(figsize=(20, 8.8)); fig.patch.set_facecolor("white")
gs = GridSpec(2, 2, height_ratios=[8, 1.2], width_ratios=[1.02, 1.0], hspace=0.06, wspace=0.14,
              left=0.05, right=0.985, top=0.9, bottom=0.08)
axm = fig.add_subplot(gs[0, 0]); axl = fig.add_subplot(gs[1, 0]); axz = fig.add_subplot(gs[0, 1])
ymax = len(chroms); BH = 0.6; maxMb = max(chromlen.values())/1e6

for yi, ch in enumerate(chroms):
    y = ymax - yi; nw = chromlen[ch]//WINDOW + 1
    axm.add_patch(FancyBboxPatch((0, y-BH/2), chromlen[ch]/1e6, BH,
                  boxstyle="round,pad=0,rounding_size=0.2", linewidth=0, facecolor=CHROM_BG, zorder=1))
    segs = []; cur = None; s0 = 0
    for w in range(nw):
        d = path[ch].get(w, used[0])
        if d != cur:
            if cur is not None: segs.append((s0, w, cur))
            cur = d; s0 = w
    segs.append((s0, nw, cur))
    for a, b, d in segs:
        axm.barh(y, (b-a)*WINDOW/1e6, left=a*WINDOW/1e6, height=BH, color=shade[d],
                 edgecolor="white", linewidth=0.35, zorder=2)
    # GWAS 位点标记(在染色体条上方)
    for (mch, mpos), cap in loci_mark.items():
        if mch != ch: continue
        axm.plot(mpos/1e6, y+BH/2+0.16, marker="v", markersize=4.3,
                 color=(CAPTURED_COLOR if cap else MISSED_COLOR),
                 markeredgecolor="white", markeredgewidth=0.3, zorder=4)
    dsv, dsnp, fc, ft = counts.get(ch, ("0","0","0","0"))
    axm.text(chromlen[ch]/1e6 + maxMb*0.012, y, f"{dsv} / {dsnp}   $\\star$ {fc}/{ft}",
             va="center", ha="left", fontsize=7.8, color="#2b3136")
    axm.text(-maxMb*0.012, y, ch, va="center", ha="right", fontsize=8.5, color="#57616a")
axm.set_xlim(-maxMb*0.05, maxMb*1.17); axm.set_ylim(0.3, ymax+1.1); axm.set_yticks([])
for s in ["top","right","left"]: axm.spines[s].set_visible(False)
axm.spines["bottom"].set_bounds(0, maxMb); axm.set_xticks(np.arange(0, maxMb+1, 25))
axm.set_xlabel("Chromosomal position (Mb)", fontsize=9)
axm.text(maxMb + maxMb*0.012, ymax+0.62, r"dSVs / dSNPs   $\star$ fav captured/total", fontsize=7.8, color="#57616a", ha="left", style="italic")
axm.set_title("a   Ideal parental haplotype (IPHs): low deleterious load + favourable GWAS alleles",
              fontsize=12.5, loc="left", pad=12, color="#1c2126", fontweight="bold")

axl.axis("off")
handles = [Patch(facecolor=shade[d], edgecolor="white", linewidth=0.3, label=d) for d in used]
mk = [Line2D([0],[0], marker="v", color="w", markerfacecolor=CAPTURED_COLOR, markersize=6, label="fav. GWAS allele captured"),
      Line2D([0],[0], marker="v", color="w", markerfacecolor=MISSED_COLOR, markersize=6, label="GWAS locus missed")]
leg1 = axl.legend(handles=handles, loc="center left", bbox_to_anchor=(0.0, 0.5), ncol=7, fontsize=6.4,
                  frameon=False, handlelength=1.1, handleheight=1.1, labelspacing=0.3, columnspacing=0.9,
                  title=f"Donor haplotype (n={len(used)})", title_fontsize=7.2)
leg1.get_title().set_color("#57616a"); axl.add_artist(leg1)
axl.legend(handles=mk, loc="center right", bbox_to_anchor=(1.0, 0.5), fontsize=6.8, frameon=False)

teal = LinearSegmentedColormap.from_list("teal", HEATMAP_COLORS)
im = axz.imshow(np.log1p(load), aspect="auto", cmap=teal, interpolation="nearest",
                extent=[0, n_win*WINDOW/1e6, len(donors)-0.5, -0.5])
axz.set_yticks(range(len(donors))); axz.set_yticklabels(donors, fontsize=6.2, color="#57616a")
axz.set_xlabel(f"{ZC} position (Mb)", fontsize=9)
axz.set_title(f"     {ZC}: candidate donor haplotypes (load per {WINDOW//1000}kb) — red = selected into IPHs; label = dSV/dSNP",
              fontsize=10, loc="left", pad=12, color="#1c2126")
segs = []; cur = None; s0 = 0
for w in range(n_win):
    d = path[ZC].get(w)
    if d != cur:
        if cur is not None: segs.append((s0, w, cur))
        cur = d; s0 = w
segs.append((s0, n_win, cur))
for a, b, d in segs:
    r = row_of.get(d)
    if r is None: continue
    x0 = a*WINDOW/1e6; wdt = (b-a)*WINDOW/1e6
    axz.add_patch(Rectangle((x0, r-0.5), wdt, 1, fill=False, edgecolor=SELECT_RED, linewidth=1.3))
    sv = int(dsv_w[d][a:b].sum()); sn = int(dsnp_w[d][a:b].sum())
    axz.text(x0+wdt/2, r-0.62, f"{sv}/{sn}", color=SELECT_TEXT_RED, fontsize=4.8, ha="center", va="bottom")
# 在放大图顶部也标 chr01B 的 GWAS 位点
for (mch, mpos), cap in loci_mark.items():
    if mch != ZC: continue
    axz.plot(mpos/1e6, -0.9, marker="v", markersize=4.5, clip_on=False,
             color=(CAPTURED_COLOR if cap else MISSED_COLOR), markeredgecolor="white", markeredgewidth=0.3)
for s in ["top","right"]: axz.spines[s].set_visible(False)
cb = fig.colorbar(im, ax=axz, fraction=0.018, pad=0.01); cb.set_label("log(1+load / window)", fontsize=7.5)
cb.ax.tick_params(labelsize=6.5)

png_path = f"{OUT}/figures/{OUT_BASENAME}.png"
pdf_path = f"{OUT}/figures/{OUT_BASENAME}.pdf"
plt.savefig(png_path, dpi=270, facecolor="white", bbox_inches="tight")
Image.open(png_path).convert("RGB").save(pdf_path, "PDF", resolution=300.0)
print(f"wrote figures/{OUT_BASENAME}.png/.pdf", file=sys.stderr)
