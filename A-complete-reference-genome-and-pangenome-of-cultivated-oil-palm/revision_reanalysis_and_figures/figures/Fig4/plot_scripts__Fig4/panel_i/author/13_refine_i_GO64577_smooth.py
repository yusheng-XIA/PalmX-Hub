#!/usr/bin/env python3
"""Smooth panel-i display within biological regions without changing raw data."""

import json
import os
from datetime import datetime, timezone
from pathlib import Path

os.environ.setdefault("MPLBACKEND", "Agg")
os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib-fig4i-go64577-smooth")

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.signal import savgol_filter

RUN = Path("${ANALYSIS_DIR}/22_answer_reviews/00_ms/03_V3/04_figure4/Fig4_d_i_pan39_material33_singletons_20260811")
TABLE = RUN / "tables/TE_density_class_position_curve_GO64577_material33.tsv"
CLASSES = ("Core", "Soft-core", "Shell", "Cloud")
COLORS = {"Core": "#82C7B8", "Soft-core": "#E65D6D", "Shell": "#29AFD4", "Cloud": "#7DCDF8"}
REGIONS = ((0, 20, 7), (20, 120, 11), (120, 140, 7))


def region_smooth(values):
    result = np.asarray(values, dtype=float).copy()
    for start, end, window in REGIONS:
        result[start:end] = savgol_filter(result[start:end], window_length=window, polyorder=2, mode="interp")
    return result


def plot(frame, labelled):
    plt.rcParams.update({
        "font.family": "sans-serif", "font.sans-serif": ["Liberation Sans", "Arial", "DejaVu Sans"],
        "font.size": 8, "axes.labelsize": 9, "axes.linewidth": 0.8,
        "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
    })
    fig, ax = plt.subplots(figsize=(5.2, 3.8), facecolor="white")
    for category in CLASSES:
        data = frame.loc[frame["Pangenome_class"].eq(category)].sort_values("Profile_bin")
        if len(data) != 140:
            raise AssertionError(f"{category} has {len(data)} bins")
        x = data["Profile_bin"].to_numpy(float) + 0.5
        y = region_smooth(data["Mean_material_TE_percent"].to_numpy(float))
        lo = np.clip(region_smooth(data["CI95_low"].to_numpy(float)), 0, 100)
        hi = np.clip(region_smooth(data["CI95_high"].to_numpy(float)), 0, 100)
        lo = np.minimum(lo, y)
        hi = np.maximum(hi, y)
        ax.plot(x, y, color=COLORS[category], lw=1.9, label=category, zorder=3,
                solid_capstyle="round", solid_joinstyle="round")
        ax.fill_between(x, lo, hi, color=COLORS[category], alpha=0.16, linewidth=0)
    ax.axvline(20, color="#333333", ls=(0, (4, 3)), lw=0.8)
    ax.axvline(120, color="#333333", ls=(0, (4, 3)), lw=0.8)
    ax.set_xlim(0, 140)
    ax.set_ylim(bottom=0)
    ax.set_xticks([0, 20, 120, 140])
    ax.set_xticklabels(["−2 kb", "TSS", "TES", "+2 kb"])
    ax.set_ylabel("TE density (%)")
    ax.legend(frameon=False, ncol=1, loc="upper center", bbox_to_anchor=(0.56, 1.0), handlelength=1.8)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(width=0.8, length=3)
    fig.subplots_adjust(left=0.14, right=0.98, bottom=0.15, top=0.96)
    if labelled:
        fig.text(0.012, 0.985, "i", ha="left", va="top", fontsize=13, fontweight="bold")
    stem = f"Fig4i_TE_density_pan39_material33_{'labelled' if labelled else 'no_label'}"
    fig.savefig(RUN / "figures" / f"{stem}.pdf", facecolor="white")
    fig.savefig(RUN / "figures" / f"{stem}.svg", facecolor="white")
    fig.savefig(RUN / "figures" / f"{stem}_600dpi.png", dpi=600, facecolor="white")
    plt.close(fig)


frame = pd.read_csv(TABLE, sep="\t")
if len(frame) != 560:
    raise AssertionError("Panel-i source table must contain 560 rows")
for labelled in (False, True):
    plot(frame, labelled)
report = {
    "status": "PASS", "completed_utc": datetime.now(timezone.utc).isoformat(),
    "raw_source_rows_unchanged": len(frame),
    "display_only_smoothing": "Savitzky-Golay, polynomial order 2; windows 7/11/7 bins for upstream/gene-body/downstream",
    "boundary_rule": "regions smoothed independently; no smoothing across TSS or TES",
}
(RUN / "provenance/GO64577_full_redraw/panel_i_smoothing_record.json").write_text(json.dumps(report, indent=2) + "\n")
print(json.dumps(report, indent=2))
