#!/usr/bin/env python3
"""Supplementary Fig. 4 redraw (metabolome LC-MS quality control).

Data: Source Data SF10a/b/c (REVIEW workbook, via ../data/f8_raw.pkl) and
F1's feature-level RSD tables SF4c_feature_RSD_{pos,neg}.tsv (read only).
"""
import sys
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, "${WORK_DIR}/redraw")  # shared f8_style (read-only)
from f8_style import BLUE, DARK, GREY, MM, ORANGE, RED, clean, sheet  # noqa: E402
sys.path.insert(0, "${WORK_DIR}/fix/beautify/common")
import palA  # noqa: E402


def letter(fig, ax, s, dx=-0.055, dy=0.012):
    """f8_style.letter with the Nature 8-pt panel letter (was 9 pt)."""
    bb = ax.get_position()
    fig.text(bb.x0 + dx, bb.y1 + dy, s, fontsize=palA.LETTER_PT, fontweight="bold", va="bottom", ha="left")

HERE = Path(__file__).resolve().parent
MODES = {"pos": "Positive-ion mode", "neg": "Negative-ion mode"}
MODE_COL = {"pos": BLUE, "neg": ORANGE}
CMAP = mpl.colormaps["viridis"]

fig = plt.figure(figsize=(180 * MM, 205 * MM))
gsa = fig.add_gridspec(2, 2, left=0.085, right=0.87, top=0.955, bottom=0.60, hspace=0.18, wspace=0.22)
gsb = fig.add_gridspec(1, 2, left=0.085, right=0.985, top=0.515, bottom=0.36, wspace=0.22)
gsc = fig.add_gridspec(1, 2, left=0.085, right=0.985, top=0.25, bottom=0.06, wspace=0.3)

# ---- a: TIC / BPC traces -------------------------------------------------------
summary = []
ax_a = {}
for j, mode in enumerate(MODES):
    tr = sheet(f"SF10a_{mode}_traces")
    tr["rt_min"] = tr.retention_time_seconds.astype(float) / 60
    n_inj = tr.injection_order.max()
    for i, (col, lab, scale, exp) in enumerate([("total_ion_current", "MS1 total ion current", 1e8, 8),
                                                 ("base_peak_intensity", "MS1 base-peak intensity", 1e7, 7)]):
        ax = fig.add_subplot(gsa[i, j])
        for (sn, st, io), g in tr.groupby(["sample_name", "sample_type", "injection_order"], sort=False):
            if st == "QC":
                continue
            ax.plot(g.rt_min, g[col] / scale, lw=0.35, color=CMAP((io - 1) / (n_inj - 1)), alpha=0.8, zorder=1)
        for (sn, st), g in tr.groupby(["sample_name", "sample_type"], sort=False):
            if st == "QC":
                ax.plot(g.rt_min, g[col] / scale, lw=0.55, color=DARK, alpha=0.9, zorder=2)
        ax.set_xlim(-0.2, 10.2)
        ax.set_ylim(bottom=0)
        if j == 0:
            ax.set_ylabel(f"{lab}\n" + rf"($\times 10^{{{exp}}}$)")
        if i == 0:
            n_bio = tr[tr.sample_type != "QC"].sample_name.nunique()
            n_qc = tr[tr.sample_type == "QC"].sample_name.nunique()
            ax.set_title(f"{MODES[mode]} ({n_bio} biological + {n_qc} pooled QC)", pad=3)
            ax.set_xticklabels([])
            summary.append((mode, n_bio, n_qc))
        else:
            ax.set_xlabel("Retention time (min)")
        clean(ax)
        ax_a[(i, j)] = ax
cax = fig.add_axes([0.89, 0.60, 0.01, 0.27])
cb = mpl.colorbar.ColorbarBase(cax, cmap=CMAP, norm=mpl.colors.Normalize(1, 128))
cb.set_label("Injection order\n(biological samples)", fontsize=6)
cb.set_ticks([1, 32, 64, 96, 128])
cax.tick_params(labelsize=5.5)
fig.legend(handles=[mpl.lines.Line2D([], [], color=DARK, lw=0.9, label="Pooled QC")],
           loc="upper left", bbox_to_anchor=(0.875, 0.93), frameon=False, fontsize=6)

# ---- b: total peak-area drift ---------------------------------------------------
ax_b = []
drift_rows = []
for j, mode in enumerate(MODES):
    q = sheet(f"SF10b_{mode}_QC")
    q["final_total_area"] = q.final_total_area.astype(float)
    med = q.loc[q.sample_type == "biological", "final_total_area"].median()
    q["rel"] = q.final_total_area / med
    ax = fig.add_subplot(gsb[0, j])
    bio = q[q.sample_type == "biological"]
    aqc = q[(q.sample_type == "QC") & (q.analysis_qc == "yes")]
    cqc = q[(q.sample_type == "QC") & (q.conditioning_qc == "yes")]
    ax.scatter(bio.injection_order, bio.rel, s=6, color=GREY, lw=0, label="Biological", zorder=2)
    ax.scatter(aqc.injection_order, aqc.rel, s=11, color=RED, lw=0, label="Pooled QC (analysis)", zorder=3)
    ax.scatter(cqc.injection_order, cqc.rel, s=11, facecolor="white", edgecolor=RED, lw=0.7,
               label="Pooled QC (conditioning)", zorder=3)
    ax.axhline(1, color=DARK, ls=":", lw=0.6)
    ax.set_xlabel("Injection order")
    if j == 0:
        ax.set_ylabel("Total peak area\n(relative to biological median)")
    ax.set_title(MODES[mode], pad=3)
    if j == 1:
        h, l = ax.get_legend_handles_labels()
        fig.legend(h, l, loc="center", bbox_to_anchor=(0.53, 0.305), ncol=3, frameon=False,
                   handletextpad=0.2, columnspacing=1.5)
    clean(ax)
    ax_b.append(ax)
    qc = q[q.sample_type == "QC"]
    drift_rows.append((mode, len(bio), len(qc), round(qc.rel.min(), 3), round(qc.rel.max(), 3),
                       round(qc.rel.std() / qc.rel.mean() * 100, 2)))
    q[["sample_name", "injection_order", "sample_type", "conditioning_qc", "analysis_qc",
       "final_total_area", "rel"]].rename(columns={"rel": "relative_to_biological_median"}).to_csv(
        HERE / f"SF4b_{mode}_total_area_relative.tsv", sep="\t", index=False)

# ---- c: RSD distribution + feature-filter counts -------------------------------
ax = fig.add_subplot(gsc[0, 0])
rsd_med = {}
bins = np.arange(0, 31, 2)
for mode in MODES:
    r = pd.read_csv(HERE / f"SF4c_feature_RSD_{mode}.tsv", sep="\t")
    keep = r[r.retained_in_final_matrix.astype(str) == "True"]
    v = keep.pooled_QC_RSD_after_correction_percent.astype(float)
    rsd_med[mode] = (len(v), v.median(), keep.pooled_QC_RSD_before_correction_percent.astype(float).median())
    ax.hist(v, bins=bins, histtype="step", lw=1.1, color=MODE_COL[mode],
            label=f"{MODES[mode].split('-')[0]} (n = {len(v):,}; median {v.median():.1f}%)")
ax.axvline(30, color=RED, ls="--", lw=0.6)
ax.set_xlim(0, 31)
ax.set_xlabel("Analysis-QC RSD after drift correction (%)")
ax.set_ylabel("Retained features")
ax.set_title("Retained features", pad=3)
ax.legend(frameon=False, loc="upper right", fontsize=5.8)
clean(ax)
ax_c = ax

ax = fig.add_subplot(gsc[0, 1])
steps = ["Aligned\nfeatures", "QC detected\n(≥9 of 11)", "Presence\nfilter", "RSD ≤ 30%\n(final)"]
flow_rows = []
w = 0.38
for k, mode in enumerate(MODES):
    f = sheet(f"SF10c_{mode}_flow")
    counts = f.set_index("Stage").feature_count.astype(int)
    vals = [counts["aligned_features"], counts["QC detected in at least 9 of 11"],
            counts["combined presence filter"], counts["post-correction QC RSD <= 30%"]]
    x = np.arange(4) + (k - 0.5) * w
    ax.bar(x, np.array(vals) / 1e3, width=w, color=MODE_COL[mode], label=MODES[mode].split("-")[0])
    for xi, v in zip(x, vals):
        ax.text(xi, v / 1e3 + 1.2, f"{v:,}", ha="center", va="bottom", fontsize=5, rotation=90)
    flow_rows.append((mode, *vals))
ax.set_xticks(np.arange(4), steps, fontsize=5.8)
ax.set_ylabel(r"Features ($\times 10^{3}$)")
ax.set_ylim(0, 100)
ax.set_title("Feature filtering", pad=3)
ax.legend(frameon=False, loc="upper right")
clean(ax)

letter(fig, ax_a[(0, 0)], "a", dx=-0.075)
letter(fig, ax_b[0], "b", dx=-0.075)
letter(fig, ax_c, "c", dx=-0.075)
for ext in ("pdf", "png"):
    fig.savefig(HERE / f"Supplementary_Fig_04.{ext}", dpi=600, facecolor="white")

pd.DataFrame(flow_rows, columns=["mode", "aligned", "QC_detected_9of11", "combined_presence", "RSD_le30_final"]
             ).to_csv(HERE / "SF4c_filter_flow.tsv", sep="\t", index=False)
print("acquisitions", summary)
print("drift (mode, n_bio, n_qc, qc_rel_min, qc_rel_max, qc_CV%)", drift_rows)
print("RSD retained (n, median after, median before)", rsd_med)
print("flow", flow_rows)
