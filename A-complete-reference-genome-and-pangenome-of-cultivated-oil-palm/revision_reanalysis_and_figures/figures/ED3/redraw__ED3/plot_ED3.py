#!/usr/bin/env python3
"""Extended Data Fig. 3 redraw (allele-specific expression) with robustness panels.

Panel sources (all tabular; nothing is pasted from rasters):
  a  src/source2_B_upset_intersections.tsv, src/source2_B_upset_set_sizes.tsv
     (identical to Source Data ED4a_upset_intersect / ED4a_upset_sizes)
  b,c src/source_03A_FL_DBA_DEBA_chromosomes.tsv, src/source_03B_TN_DBA_DEBA_chromosomes.tsv
     (DBA = eligible_genes and DEBA = ASE_genes of Source Data ED4bc_ASE_windows)
  d  src/source_02D_trait_gene_atlas_all.tsv (gene-level recurrent-ASE classes)
  e,f Supplementary Table 14, diagnostic-marker concordance and robustness blocks (rows 51-245)
Outputs: out/Extended_Data_Fig_03.{pdf,png} and data/ED3*_*.tsv source tables.
"""
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import openpyxl
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Patch, Rectangle

HERE = Path(__file__).resolve().parent
SRC = HERE / "src"
OUT = HERE / "out"
DATA = HERE / "data"
OUT.mkdir(exist_ok=True)
DATA.mkdir(exist_ok=True)
ST = HERE.parent.parent / "Supplementary_Tables_submission_114RNAseq.xlsx"

MM = 1 / 25.4
DARK, GREY = "#242A30", "#B9C0C7"
FL_COL, TN_COL = "#2A9D8F", "#D95F5F"
CLASS_COLOR = {"HapDom": "#F2D65C", "Sub": "#8DD3C7", "NoDiff": "#F8766D", "NoASE": "#4EA5DF"}
ATLAS_COLOR = {"A recurrent": "#E76F61", "B recurrent": "#4C9ED1",
               "Switching recurrent": "#F3C85B", "No/stage-limited ASE": "#D2D5D7"}

mpl.rcParams.update({
    "font.family": "Arial", "font.size": 6.5, "axes.labelsize": 6.5, "axes.titlesize": 6.5,
    "xtick.labelsize": 5.8, "ytick.labelsize": 5.8, "legend.fontsize": 5.8,
    "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6,
    "xtick.major.size": 2.2, "ytick.major.size": 2.2, "axes.unicode_minus": True,
    "mathtext.fontset": "custom", "mathtext.rm": "Arial", "mathtext.it": "Arial:italic",
    "mathtext.bf": "Arial:bold", "pdf.fonttype": 42, "svg.fonttype": "none",
})


def clean(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def letter(fig, x, y, s):
    fig.text(x, y, s, fontsize=9, fontweight="bold", va="top", ha="left")


# ------------------------------------------------------------------ data
ups = pd.read_csv(SRC / "source2_B_upset_intersections.tsv", sep="\t", dtype={"pattern": str})
sizes = pd.read_csv(SRC / "source2_B_upset_set_sizes.tsv", sep="\t")
chrom = {m: pd.read_csv(SRC / f"source_03{k}_{m}_DBA_DEBA_chromosomes.tsv", sep="\t")
         for k, m in (("A", "FL"), ("B", "TN"))}
atlas = pd.read_csv(SRC / "source_02D_trait_gene_atlas_all.tsv", sep="\t")
# A/B follow the source definition (as Fig. 3c): FL A = FL-Hap2 (Africa hap2), B = FL-Hap1 (American hap1);
# TN A = TK (dura)-like, B = NS (pisifera)-like. ASE_status is the source class, unchanged.
atlas["A_B_definition"] = np.where(atlas.analysis == "FL", "A=FL-Hap2; B=FL-Hap1",
                                   "A=TK (dura)-like; B=NS (pisifera)-like")

wb = openpyxl.load_workbook(ST, read_only=True)
rows = [list(r) for r in wb["Supplementary Table 14"].iter_rows(values_only=True)]


def block(start_label, stop_blank=True):
    i = next(k for k, r in enumerate(rows) if r[0] and str(r[0]).startswith(start_label))
    hdr = [c for c in rows[i + 1] if c is not None]
    out = []
    for r in rows[i + 2:]:
        vals = r[:len(hdr)]
        if vals[0] is None or str(vals[0]).startswith("Note"):
            break
        out.append(vals)
    return pd.DataFrame(out, columns=hdr)


overall = block("Diagnostic-marker concordance between")
sample = block("Sample-level diagnostic-marker concordan")
robust = block("Robustness of ASE classification after")
consist = block("ASE classification consistency metrics")
trans = block("Overall ASE class transitions among gene")
for df in (overall, sample, robust, consist, trans):
    for c in df.columns:
        df[c] = pd.to_numeric(df[c], errors="ignore")

# export panel source tables
ups.to_csv(DATA / "ED3a_upset_intersections.tsv", sep="\t", index=False)
sizes.to_csv(DATA / "ED3a_upset_set_sizes.tsv", sep="\t", index=False)
pd.concat(chrom.values()).to_csv(DATA / "ED3bc_DBA_DEBA_1Mb_windows.tsv", sep="\t", index=False)
atlas.to_csv(DATA / "ED3d_lipid_gene_ASE_atlas.tsv", sep="\t", index=False)
sample.to_csv(DATA / "ED3e_sample_level_marker_concordance.tsv", sep="\t", index=False)
overall.to_csv(DATA / "ED3e_overall_marker_concordance.tsv", sep="\t", index=False)
robust.to_csv(DATA / "ED3f_class_composition_before_after_filtering.tsv", sep="\t", index=False)
consist.to_csv(DATA / "ED3f_classification_consistency.tsv", sep="\t", index=False)
trans.to_csv(DATA / "ED3f_class_transitions.tsv", sep="\t", index=False)

fig = plt.figure(figsize=(180 * MM, 170 * MM))

# ------------------------------------------------------------------ a UpSet
N_SHOW = 18
p = ups.head(N_SHOW).reset_index(drop=True)
set_order = ["FL_Early", "FL_Mid", "FL_Late", "FL_Postharvest",
             "TN_Early", "TN_Mid", "TN_Late", "TN_Postharvest"]
set_lab = ["FL Early", "FL Mid", "FL Late", "FL Post", "TN Early", "TN Mid", "TN Late", "TN Post"]
B, T = 0.655, 0.965
ax_bar = fig.add_axes([0.145, B + 0.16, 0.2, T - B - 0.16])
ax_dot = fig.add_axes([0.145, B, 0.2, 0.15], sharex=ax_bar)
ax_set = fig.add_axes([0.03, B, 0.05, 0.15], sharey=ax_dot)
x = np.arange(len(p))
ax_bar.bar(x, p.orthogroups, color=DARK, width=0.62)
ax_bar.set_ylabel("Intersection size\n(orthogroups)")
ax_bar.tick_params(axis="x", bottom=False, labelbottom=False)
clean(ax_bar)
ypos = {s: 7 - i for i, s in enumerate(set_order)}
for i, pat in enumerate(p.pattern):
    pat = pat.zfill(8)
    on = [ypos[s] for s, bit in zip(set_order, pat) if bit == "1"]
    ax_dot.scatter([i] * 8, list(ypos.values()), s=7, color="#E3E5E7", zorder=1, lw=0)
    ax_dot.scatter([i] * len(on), on, s=8, color=DARK, zorder=3, lw=0)
    if len(on) > 1:
        ax_dot.plot([i, i], [min(on), max(on)], color=DARK, lw=0.8, zorder=2)
for k in range(4):
    ax_dot.axhspan(7 - k - 0.5, 7 - k + 0.5, color="#EEF7F6", zorder=0, lw=0)
ax_dot.set_xlim(-0.6, len(p) - 0.4)
ax_dot.set_ylim(-0.6, 7.6)
ax_dot.axis("off")
sz = sizes.set_index("set").orthogroups
ax_set.barh([ypos[s] for s in set_order], [sz[s] for s in set_order], height=0.62,
            color=[FL_COL] * 4 + [TN_COL] * 4, alpha=0.75)
ax_set.invert_xaxis()
ax_set.set_yticks([])
for s_, lab_ in zip(set_order, set_lab):
    ax_dot.text(-0.9, ypos[s_], lab_, ha="right", va="center", fontsize=5.8, clip_on=False)
ax_set.set_xticks([6000, 0], ["6,000", "0"])
ax_set.set_xlabel("Set size", labelpad=1)
for s_ in ("top", "left", "right"):
    ax_set.spines[s_].set_visible(False)
letter(fig, 0.012, 0.985, "a")

# ------------------------------------------------------------------ b,c chromosomes
cm_dba = LinearSegmentedColormap.from_list("dba", ["#F4F5F4", "#9FCB7E", "#2E7D32", "#123F16"])
cm_deba = LinearSegmentedColormap.from_list("deba", ["#F7F2F0", "#F2A58E", "#D9502F", "#8C2A12"])
for k, (m, anchor) in enumerate((("FL", "FL-Hap2 (African-derived) coordinates"),
                                 ("TN", "TN-Hap1 (dura/TK-like) coordinates"))):
    d = chrom[m]
    x0 = 0.41 + k * 0.275
    ax = fig.add_axes([x0, 0.66, 0.225, 0.305])
    vmax = float(np.ceil(pd.concat(chrom.values()).DBA_genes.quantile(0.995)))
    lengths = d.groupby("chrom_num").window_mb.max() + 1
    for c in range(1, 17):
        q = d[d.chrom_num == c].set_index("window_mb")
        y = 16 - c
        for col, cmap, dy in (("DBA_genes", cm_dba, 0.2), ("DEBA_genes", cm_deba, -0.2)):
            vals = q[col].reindex(range(int(lengths[c])), fill_value=0).values
            ax.imshow(vals[None, :], aspect="auto", cmap=cmap, vmin=0, vmax=vmax,
                      extent=(0, len(vals), y + dy - 0.18, y + dy + 0.18), interpolation="nearest")
    ax.set_xlim(0, lengths.max() + 1)
    ax.set_ylim(-0.6, 15.6)
    ax.set_yticks(range(16), [f"chr{c:02d}" for c in range(16, 0, -1)])
    ax.tick_params(axis="y", length=0, pad=1.5)
    ax.set_xlabel("Position (Mb)")
    for s in ("top", "right", "left"):
        ax.spines[s].set_visible(False)
    ax.set_title(f"{m} (1-Mb windows)", pad=3)
    if k == 0:
        continue
    cax1 = fig.add_axes([0.93, 0.83, 0.008, 0.12])
    cax2 = fig.add_axes([0.93, 0.68, 0.008, 0.12])
    for cax, cmap, lab in ((cax1, cm_dba, "DBA"), (cax2, cm_deba, "DEBA")):
        cb = fig.colorbar(mpl.cm.ScalarMappable(mpl.colors.Normalize(0, vmax), cmap), cax=cax)
        cb.outline.set_linewidth(0.4)
        cb.ax.tick_params(labelsize=5, length=1.5, width=0.4)
        cb.set_label(f"{lab} genes", fontsize=5.2, labelpad=1)
for k, lab in enumerate("bc"):
    letter(fig, 0.41 + k * 0.275 - 0.045, 0.985, lab)

# ------------------------------------------------------------------ d gene atlas
FAMILY_GROUPS = [
    ("Fatty-acid synthesis", ["ACCase", "ACP", "FabD (MCAT)", "FabG (KAR)", "FabI (ENR)",
                              "KASIII", "KAS I/II", "FATA/B", "LACS"]),
    ("Desaturation and elongation", ["SAD", "FAD2", "FAD3", "FAD6", "FAD7/FAD8", "KCS"]),
    ("TAG assembly", ["GPAT", "LPAT", "PAP", "DGAT", "PDAT", "PDCT"]),
    ("Lipid oxidation and antioxidants", ["LOX9", "VTE1", "APX7"]),
]
fams = [f for _, fs in FAMILY_GROUPS for f in fs]
order = {"A recurrent": 0, "B recurrent": 1, "Switching recurrent": 2, "No/stage-limited ASE": 3}
ax = fig.add_axes([0.13, 0.045, 0.47, 0.545])
maxg = int(atlas.groupby(["analysis", "family"]).size().max())
nrow = len(fams) * 2
yy = 0
ytick, ylab = [], []
for gi, (gname, fs) in enumerate(FAMILY_GROUPS):
    g0 = yy
    for f in fs:
        for m in ("FL", "TN"):
            sub = atlas[(atlas.analysis == m) & (atlas.family == f)].copy()
            sub["o"] = sub.ASE_status.map(order)
            sub = sub.sort_values(["o", "robust_stages"], ascending=[True, False])
            y = nrow - 1 - yy
            for j, st in enumerate(sub.ASE_status):
                ax.add_patch(Rectangle((j + 0.06, y - 0.36), 0.88, 0.72, facecolor=ATLAS_COLOR[st],
                                       edgecolor="none"))
            ax.text(-0.35, y, m, ha="right", va="center", fontsize=4.8, color="#555")
            ax.text(maxg + 0.3, y, f"{len(sub)}", ha="left", va="center", fontsize=4.8, color="#555")
            yy += 1
        ytick.append(nrow - 1 - (yy - 1.5))
        ylab.append(f)
        ax.axhline(nrow - 1 - yy + 0.5, color="#E6E8EA", lw=0.4)
    xg = maxg + 1.6
    ax.plot([xg, xg], [nrow - 1 - g0 + 0.45, nrow - 1 - (yy - 1) - 0.45], color=DARK, lw=0.6,
            clip_on=False)
    ax.text(xg + 0.35, nrow - 1 - (g0 + yy - 1) / 2, gname.replace(" and ", "\nand "), rotation=-90,
            ha="center", va="center", fontsize=5.2, clip_on=False, linespacing=1.0)
ax.set_yticks(ytick, ylab, fontsize=5.4)
ax.tick_params(axis="y", length=0, pad=9)
ax.set_xlim(-0.2, maxg + 1.2)
ax.set_ylim(-0.6, nrow - 0.4)
ax.set_xticks([])
for s in ("top", "right", "bottom", "left"):
    ax.spines[s].set_visible(False)
ax.text(maxg + 0.3, nrow - 0.2, "n", ha="left", va="bottom", fontsize=5, color="#555")
ax.set_xlabel("Genes in family (one tile per gene)", labelpad=2)
ax.legend(handles=[Patch(facecolor=ATLAS_COLOR[k], label=l) for k, l in (
    ("A recurrent", "Recurrent A bias"), ("B recurrent", "Recurrent B bias"),
    ("Switching recurrent", "Recurrent, direction switching"),
    ("No/stage-limited ASE", "No or stage-limited ASE"))],
          ncol=4, frameon=False, loc="lower center", bbox_to_anchor=(0.40, 0.997), handlelength=1.0,
          handleheight=0.8, columnspacing=0.8, handletextpad=0.3)
letter(fig, 0.012, 0.625, "d")

# ------------------------------------------------------------------ e concordance
ax = fig.add_axes([0.735, 0.395, 0.25, 0.19])
rng = np.random.default_rng(7)
col = "Expected-allele concordance (%)"
ov = overall.set_index("Material")
for i, (m, c) in enumerate((("FL", FL_COL), ("TN", TN_COL))):
    v = sample.loc[sample.Material == m, col].astype(float).values
    ax.scatter(i + rng.uniform(-0.17, 0.17, len(v)), v, s=4, color=c, alpha=0.6, lw=0, zorder=2)
    ax.hlines(np.median(v), i - 0.25, i + 0.25, color=DARK, lw=0.8, zorder=3)
    ax.scatter(i + 0.38, float(ov.loc[m, "Expected-allele concordance (%)"]), marker="D", s=11,
               facecolor="white", edgecolor=DARK, lw=0.6, zorder=4)
    ax.scatter(i + 0.38, float(ov.loc[m, "Filtered expected-allele concordance (%)"]), marker="D", s=11,
               color=DARK, lw=0, zorder=4)
ax.set_xticks([0, 1], [f"FL\n({(sample.Material == 'FL').sum()} libraries)",
                       f"TN\n({(sample.Material == 'TN').sum()} libraries)"])
ax.set_xlim(-0.55, 1.65)
ax.set_ylabel("Expected-allele concordance (%)")
ax.set_ylim(99.5, 100.0)
ax.legend(handles=[
    mpl.lines.Line2D([], [], marker="o", ls="", color=GREY, ms=2.5, label="Per library"),
    mpl.lines.Line2D([], [], marker="D", ls="", mfc="white", mec=DARK, ms=3, label="Overall, all markers"),
    mpl.lines.Line2D([], [], marker="D", ls="", color=DARK, ms=3, label="Overall, filtered markers")],
    frameon=False, loc="lower left", handletextpad=0.2, borderaxespad=0.1, labelspacing=0.3)
clean(ax)
ax.set_title("Diagnostic-marker concordance", pad=3)
letter(fig, 0.665, 0.625, "e")

# ------------------------------------------------------------------ f robustness
ax = fig.add_axes([0.735, 0.085, 0.17, 0.2])
rb = robust.set_index("Metric")
cls_rows = [("HapDom", "HapDom (%)"), ("Sub", "Sub (%)"),
            ("NoDiff", "Stage-limited / weak-switching (%)"), ("NoASE", "NoASE (%)")]
cols = ["FL before filtering", "FL after filtering", "TN before filtering", "TN after filtering"]
xpos = [0, 0.8, 2.0, 2.8]
for xp, c in zip(xpos, cols):
    bottom = 0
    for cls, key in cls_rows:
        v = float(rb.loc[key, c])
        ax.bar(xp, v, bottom=bottom, width=0.65, color=CLASS_COLOR[cls], edgecolor="white", lw=0.4)
        bottom += v
    n = int(rb.loc["Testable genes (n)", c])
    ax.text(xp, 101.5, f"{n:,}", ha="center", va="bottom", fontsize=4.6)
ax.set_xticks(xpos, ["Before", "After", "Before", "After"])
for xm, m in ((0.4, "FL"), (2.4, "TN")):
    ax.text(xm, -12.5, m, ha="center", va="top", fontsize=6.3, transform=ax.transData)
cs = consist[consist.Level == "gene_overall"].pivot(index="Material", columns="Metric", values="Value")
agree_txt = "Same class after\nfiltering:\n" + "\n".join(
    f"{m} {float(cs.loc[m, 'exact_class_agreement_pct']):.1f}%\n($\\kappa$ = {float(cs.loc[m, 'cohen_kappa']):.2f})"
    for m in ("FL", "TN"))
ax.text(1.04, 0.36, agree_txt, transform=ax.transAxes, ha="left", va="top", fontsize=5.0, linespacing=1.15)
ax.set_ylim(0, 100)
ax.set_ylabel("Testable genes (%)")
ax.set_xlim(-0.5, 3.3)
clean(ax)
ax.text(-0.55, 101.5, "n", ha="right", va="bottom", fontsize=4.6)
ax.legend(handles=[Patch(facecolor=CLASS_COLOR[c], label=l) for c, l in (
    ("HapDom", "HapDom"), ("Sub", "Sub"), ("NoDiff", "NoDiff"), ("NoASE", "NoASE"))],
    frameon=False, loc="upper left", bbox_to_anchor=(1.0, 1.0), handlelength=0.9, handleheight=0.8,
    handletextpad=0.3, labelspacing=0.35)
letter(fig, 0.665, 0.345, "f")

for ext in ("pdf", "png"):
    fig.savefig(OUT / f"Extended_Data_Fig_03.{ext}", dpi=600, facecolor="white")

# ------------------------------------------------------------------ recomputed values for README
print("overall", overall[["Material", "Expected-allele concordance (%)",
                          "Filtered expected-allele concordance (%)"]].to_dict("records"))
for m in ("FL", "TN"):
    v = sample.loc[sample.Material == m, col].astype(float)
    print(m, "libs", len(v), "median", v.median(), "min", v.min(), "max", v.max())
print(robust.to_string())
print(cs.to_string())
t = trans.copy()
t["Genes (n)"] = t["Genes (n)"].astype(int)
for m in ("FL", "TN"):
    tm = t[t.Material == m]
    tot = tm["Genes (n)"].sum()
    same = tm.loc[tm["Class before"] == tm["Class after"], "Genes (n)"].sum()
    print(m, "transition total", tot, "same", same, f"{same / tot * 100:.2f}%")
print("atlas", atlas.groupby(["analysis", "ASE_status"]).size().to_dict())
print("upset top", ups.head(3).to_dict("records"), "shown", N_SHOW, "of", len(ups))
