#!/usr/bin/env python3
"""Aggregate current WGD outputs and render Fig. 4h in the accepted style."""

from __future__ import annotations

import csv
import json
import os
from pathlib import Path

os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("MPLBACKEND", "Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Patch
from scipy.stats import gaussian_kde, kruskal, mannwhitneyu


RUN = Path("${ANALYSIS_DIR}/22_answer_reviews/00_ms/03_V3/04_figure4/Fig4_d_i_pan39_material33_redraw_20260808")
ATTEMPT = os.environ.get("ATTEMPT_ID", "attempt1")
WORK = RUN / "work/h_wgd" / ATTEMPT
CLASSES = ("Core", "Soft-core", "Shell", "Cloud")
SLUG = {"Core": "core", "Soft-core": "softcore", "Shell": "shell", "Cloud": "cloud"}
FILLS = {"Core": "#E7A8B8", "Soft-core": "#C3B6E6", "Shell": "#98B0D1", "Cloud": "#A7D9DD"}
LINES = {"Core": "#C8788C", "Soft-core": "#9B8CC3", "Shell": "#6E8CAF", "Cloud": "#78B9BE"}


def setup_style() -> None:
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Liberation Sans", "Arial", "DejaVu Sans"],
        "font.size": 8.0,
        "axes.labelsize": 9.0,
        "axes.linewidth": 0.65,
        "xtick.labelsize": 7.5,
        "ytick.labelsize": 7.5,
        "xtick.major.width": 0.6,
        "ytick.major.width": 0.6,
        "xtick.major.size": 2.8,
        "ytick.major.size": 2.8,
        "legend.fontsize": 7.0,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "svg.fonttype": "none",
    })


def clean_axis(ax) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_color("#777777")
    ax.spines["bottom"].set_color("#777777")
    ax.tick_params(colors="#202020")


def read_outputs() -> tuple[dict[str, pd.DataFrame], list[dict[str, object]]]:
    data: dict[str, pd.DataFrame] = {}
    summary: list[dict[str, object]] = []
    combined = RUN / "tables/Fig4h_WGD_pairwise_current.tsv"
    first = True
    for category in CLASSES:
        path = WORK / "wgd_output" / SLUG[category] / "all_selected_cds.fa.ks.tsv"
        if not path.is_file() or path.stat().st_size == 0:
            raise FileNotFoundError(f"Missing WGD output: {path}")
        frame = pd.read_csv(path, sep="\t", low_memory=False)
        required = {"Ka", "Ks", "Omega", "Family", "Paralog1", "Paralog2"}
        missing = required - set(frame.columns)
        if missing:
            raise ValueError(f"{path} lacks columns: {sorted(missing)}")
        for column in ("Ka", "Ks", "Omega"):
            frame[column] = pd.to_numeric(frame[column], errors="coerce")
        frame.insert(0, "Category", category)
        frame.to_csv(combined, sep="\t", index=False, mode="w" if first else "a", header=first)
        first = False
        valid = frame.dropna(subset=["Ka", "Ks", "Omega"]).copy()
        ks = valid.loc[(valid["Ks"] > 0) & (valid["Ks"] <= 3), "Ks"].to_numpy(float)
        omega = valid.loc[(valid["Ks"] > 0.01) & (valid["Ks"] <= 3) &
                          (valid["Omega"] > 0) & (valid["Omega"] < 5), "Omega"].to_numpy(float)
        if len(ks) < 10 or len(omega) < 10:
            raise RuntimeError(f"Too few valid observations for {category}: Ks={len(ks)}, Ka/Ks={len(omega)}")
        data[category] = pd.DataFrame({"Ks": pd.Series(ks), "Omega": pd.Series(omega)})
        summary.append({
            "Category": category,
            "all_pairwise_rows": int(len(frame)),
            "complete_numeric_rows": int(len(valid)),
            "Ks_0_to_3_n": int(len(ks)),
            "Ks_median": float(np.median(ks)),
            "Ks_q05": float(np.quantile(ks, 0.05)),
            "Ks_q25": float(np.quantile(ks, 0.25)),
            "Ks_q75": float(np.quantile(ks, 0.75)),
            "Ks_q95": float(np.quantile(ks, 0.95)),
            "KaKs_filtered_n": int(len(omega)),
            "KaKs_median": float(np.median(omega)),
            "KaKs_mean": float(np.mean(omega)),
            "KaKs_q05": float(np.quantile(omega, 0.05)),
            "KaKs_q25": float(np.quantile(omega, 0.25)),
            "KaKs_q75": float(np.quantile(omega, 0.75)),
            "KaKs_q95": float(np.quantile(omega, 0.95)),
        })
    return data, summary


def fdr_bh(values: list[float]) -> list[float]:
    p = np.asarray(values, dtype=float)
    order = np.argsort(p)
    ranked = p[order]
    adjusted = np.minimum.accumulate((ranked * len(p) / np.arange(1, len(p) + 1))[::-1])[::-1]
    out = np.empty_like(adjusted)
    out[order] = np.minimum(adjusted, 1.0)
    return out.tolist()


def write_statistics(data: dict[str, pd.DataFrame], summary: list[dict[str, object]]) -> dict[str, object]:
    summary_path = RUN / "tables/Fig4h_WGD_summary.tsv"
    with summary_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(summary[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(summary)

    groups = [data[c]["Omega"].dropna().to_numpy() for c in CLASSES]
    h_stat, kw_p = kruskal(*groups)
    rows = []
    raw_p = []
    for other in CLASSES[1:]:
        stat, p = mannwhitneyu(groups[0], data[other]["Omega"].dropna().to_numpy(), alternative="two-sided")
        raw_p.append(float(p))
        rows.append({"comparison": f"Core vs {other}", "U": float(stat), "P_raw": float(p)})
    for row, adj in zip(rows, fdr_bh(raw_p)):
        row["P_BH"] = adj
    stats_path = RUN / "tables/Fig4h_WGD_statistics.tsv"
    with stats_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["comparison", "U", "P_raw", "P_BH"], delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    return {"Kruskal_Wallis_H": float(h_stat), "Kruskal_Wallis_P": float(kw_p), "pairwise": rows}


def significance(p: float) -> str:
    if p < 0.001:
        return "***"
    if p < 0.01:
        return "**"
    if p < 0.05:
        return "*"
    return "ns"


def save_plot(data: dict[str, pd.DataFrame], stats: dict[str, object], labelled: bool) -> None:
    setup_style()
    fig, ax = plt.subplots(figsize=(5.4, 3.95), facecolor="white")

    # Display the complete accepted Ks interval on a log10 coordinate so that
    # the biologically informative near-zero region is not compressed against
    # the y-axis. Tick labels remain in the original Ks units.
    ks_ticks = np.asarray([0.001, 0.003, 0.01, 0.03, 0.1, 0.3, 1.0, 3.0])
    log_grid = np.linspace(np.log10(ks_ticks[0]), np.log10(ks_ticks[-1]), 501)
    for category in CLASSES:
        values = data[category]["Ks"].dropna().to_numpy()
        display_values = values[values >= ks_ticks[0]]
        kde = gaussian_kde(np.log10(display_values), bw_method="scott")
        density = kde(log_grid)
        ax.plot(log_grid, density, color=LINES[category], lw=1.35, label=category)
    ax.set_xlim(log_grid[0], log_grid[-1])
    ax.set_ylim(bottom=0)
    ax.set_xticks(np.log10(ks_ticks))
    ax.set_xticklabels(["0.001", "0.003", "0.01", "0.03", "0.1", "0.3", "1", "3"], rotation=35, ha="right")
    ax.set_xlabel(r"$K_s$")
    ax.set_ylabel(r"Density (log$_{10}$ $K_s$)")
    clean_axis(ax)
    legend_handles = [Patch(facecolor="white", edgecolor=LINES[c], linewidth=1.2, label=c) for c in CLASSES]
    ax.legend(handles=legend_handles, frameon=False, loc="lower center",
              bbox_to_anchor=(0.50, 1.015), ncol=4, handlelength=1.1, handleheight=0.8,
              borderaxespad=0, columnspacing=0.9, handletextpad=0.35)

    # A large inset keeps h as one assembly-ready panel while leaving clear
    # space from the main x/y axes. It occupies the low-information upper-right
    # portion of the Ks density field.
    omega_ax = ax.inset_axes([0.585, 0.48, 0.39, 0.43])
    arrays = [data[c]["Omega"].dropna().to_numpy() for c in CLASSES]
    boxes = omega_ax.boxplot(arrays, positions=np.arange(1, 5), widths=0.58, patch_artist=True,
                             showfliers=False, whis=(5, 95),
                             medianprops={"color": "#333333", "linewidth": 0.9},
                             whiskerprops={"color": "#666666", "linewidth": 0.65},
                             capprops={"color": "#666666", "linewidth": 0.65})
    for box, category in zip(boxes["boxes"], CLASSES):
        box.set_facecolor(FILLS[category])
        box.set_edgecolor(LINES[category])
        box.set_alpha(0.72)
        box.set_linewidth(0.7)
    omega_ax.axhline(1, color="#777777", lw=0.7, ls=(0, (3, 2)), zorder=0)
    q95_max = max(float(np.quantile(values, 0.95)) for values in arrays)
    omega_top = max(1.25, q95_max * 1.24)
    omega_ax.set_ylim(0, omega_top)
    omega_ax.set_ylabel(r"$K_a/K_s$", fontsize=7.0, labelpad=1.0)
    omega_ax.set_xticks([1, 2, 3, 4])
    omega_ax.set_xticklabels(["Core", "Soft-core", "Shell", "Cloud"], rotation=28, ha="right", fontsize=5.8)
    omega_ax.tick_params(axis="y", labelsize=6.0, width=0.5, length=2.2)
    omega_ax.tick_params(axis="x", width=0.5, length=2.2, pad=1.0)
    clean_axis(omega_ax)
    annotation_y = omega_top * 0.965
    for xpos, row in enumerate(stats["pairwise"], start=2):
        omega_ax.text(xpos, annotation_y, significance(float(row["P_BH"])),
                      ha="center", va="top", fontsize=7.2, color="#333333", fontweight="bold")
    omega_ax.text(0.02, 0.995, "vs Core (BH-adjusted)", transform=omega_ax.transAxes,
                  ha="left", va="top", fontsize=5.6, color="#555555")

    fig.subplots_adjust(left=0.14, right=0.98, bottom=0.19, top=0.86)
    if labelled:
        fig.text(0.012, 0.99, "h", ha="left", va="top", fontsize=12.5, fontweight="bold")
    stem = f"Fig4h_WGD_Ks_KaKs_material33_{'labelled' if labelled else 'no_label'}"
    fig.savefig(RUN / "figures" / f"{stem}.pdf", facecolor="white")
    fig.savefig(RUN / "figures" / f"{stem}.svg", facecolor="white")
    fig.savefig(RUN / "figures" / f"{stem}_600dpi.png", dpi=600, facecolor="white")
    plt.close(fig)


def main() -> None:
    data, summary = read_outputs()
    stats = write_statistics(data, summary)
    save_plot(data, stats, labelled=True)
    save_plot(data, stats, labelled=False)
    audit = {
        "attempt": ATTEMPT,
        "status": "PASS",
        "classes": list(CLASSES),
        "class_order": list(CLASSES),
        "layout": "single Ks main panel with a large Ka/Ks inset positioned away from the main axes",
        "filters": {"Ks_analysis": "0 < Ks <= 3", "Ks_display": "0.001 <= Ks <= 3, plotted on a log10 coordinate with ticks labelled in original Ks units", "KaKs_boxplot": "0.01 < Ks <= 3 and 0 < Omega < 5; whiskers show the 5th and 95th percentiles"},
        "colors": {c: {"fill": FILLS[c], "line": LINES[c]} for c in CLASSES},
        "summary": summary,
        "statistics": stats,
        "font": "Liberation Sans (Arial metric-compatible)",
    }
    (RUN / "provenance" / f"Fig4h_audit.{ATTEMPT}.json").write_text(json.dumps(audit, indent=2) + "\n")
    (RUN / "provenance" / f"Fig4h_plot.{ATTEMPT}.SUCCESS").touch()
    print(json.dumps(audit, indent=2))


if __name__ == "__main__":
    main()
