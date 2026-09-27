#!/usr/bin/env python3
"""Replot K=4 LD decay from 2-kb-thinned PopLDdecay stats as smooth curves.

The 2-kb thinning is useful for finishing PopLDdecay quickly, but plotting the
native distance rows directly creates periodic oscillations.  This script keeps
the underlying PopLDdecay results unchanged, then produces publication-style
curves by distance binning and pair-count weighted averaging.
"""

from __future__ import annotations

import gzip
from pathlib import Path

import matplotlib as mpl

mpl.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path("${CLUSTER_WORK}/snp_repair_rerun/ld")
STATS_DIR = ROOT / "stats"
FIG_DIR = ROOT / "figures_smooth"
TABLE_DIR = ROOT / "tables_smooth"

GROUPS = [
    ("K4_Pop1", "K4-Pop1", 89, "#70B4E7"),
    ("K4_Pop2", "K4-Pop2", 23, "#F3766A"),
    ("K4_Pop3", "K4-Pop3", 153, "#82C7B8"),
    ("K4_Pop4", "K4-Pop4", 43, "#BBA8D8"),
]


def read_stat(stat_path: Path, group_id: str, label: str, sample_n: int) -> pd.DataFrame:
    rows: list[dict[str, float | str | int]] = []
    with gzip.open(stat_path, "rt") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 6:
                continue
            if parts[1] == "NA" or parts[5] == "NA":
                continue
            dist_bp = int(parts[0])
            mean_r2 = float(parts[1])
            pairs = int(parts[5])
            if pairs <= 0:
                continue
            rows.append(
                {
                    "group_id": group_id,
                    "label": label,
                    "sample_n": sample_n,
                    "dist_kb": dist_bp / 1000.0,
                    "mean_r2": mean_r2,
                    "pairs": pairs,
                }
            )
    if not rows:
        raise RuntimeError(f"No LD rows read from {stat_path}")
    return pd.DataFrame(rows)


def distance_bin(dist_kb: float) -> float:
    """Distance bins tuned for 2-kb-thinned SNPs.

    Fine bins on thinned input create artificial waves.  The early decay is still
    kept at 2-kb resolution, while longer distances use wider bins.
    """
    if dist_kb <= 100:
        return round(dist_kb / 2.0) * 2.0
    if dist_kb <= 250:
        return round(dist_kb / 5.0) * 5.0
    return round(dist_kb / 10.0) * 10.0


def build_smooth_table(raw: pd.DataFrame) -> pd.DataFrame:
    df = raw[(raw["dist_kb"] > 0) & (raw["dist_kb"] <= 500)].copy()
    df["bin_kb"] = df["dist_kb"].map(distance_bin)
    df = df[df["bin_kb"] > 0].copy()
    df["weighted_r2"] = df["mean_r2"] * df["pairs"]

    binned = (
        df.groupby(["group_id", "label", "sample_n", "bin_kb"], as_index=False)
        .agg(weighted_r2_sum=("weighted_r2", "sum"), pairs=("pairs", "sum"))
        .rename(columns={"bin_kb": "dist_kb"})
    )
    binned["r2_binned"] = binned["weighted_r2_sum"] / binned["pairs"]

    smooth_frames = []
    for group_id, group_df in binned.groupby("group_id", sort=False):
        group_df = group_df.sort_values("dist_kb").reset_index(drop=True)
        # Pair-count weighted distance bins are used for plotting.  Keeping the
        # binned value as the plotted value preserves the initial LD level used
        # for LD50 estimation while removing the 2-kb thinning oscillation.
        group_df["r2_smooth"] = group_df["r2_binned"].clip(lower=0)
        smooth_frames.append(group_df)

    out = pd.concat(smooth_frames, ignore_index=True)
    out = out[
        [
            "group_id",
            "label",
            "sample_n",
            "dist_kb",
            "pairs",
            "r2_binned",
            "r2_smooth",
        ]
    ]
    return out


def estimate_ld50(curve: pd.DataFrame) -> tuple[float, float, str]:
    curve = curve.sort_values("dist_kb").reset_index(drop=True)
    initial = float(curve["r2_smooth"].iloc[0])
    half = initial / 2.0

    above = curve[curve["r2_smooth"] >= half]
    below = curve[curve["r2_smooth"] < half]
    if below.empty:
        return initial, half, ">500"

    first_below_idx = int(below.index[0])
    if first_below_idx == 0:
        return initial, half, f"{curve.loc[first_below_idx, 'dist_kb']:.1f}"

    left = curve.loc[first_below_idx - 1]
    right = curve.loc[first_below_idx]
    x1, y1 = float(left["dist_kb"]), float(left["r2_smooth"])
    x2, y2 = float(right["dist_kb"]), float(right["r2_smooth"])
    if np.isclose(y1, y2):
        ld50 = x2
    else:
        ld50 = x1 + (half - y1) * (x2 - x1) / (y2 - y1)
    return initial, half, f"{ld50:.1f}"


def write_ld50_table(smooth: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for group_id, label, sample_n, _color in GROUPS:
        curve = smooth[smooth["group_id"] == group_id]
        initial, half, ld50 = estimate_ld50(curve)
        rows.append(
            {
                "group_id": group_id,
                "label": label,
                "sample_n": sample_n,
                "initial_smooth_r2": initial,
                "half_initial_smooth_r2": half,
                "LD50_kb_smooth": ld50,
            }
        )
    summary = pd.DataFrame(rows)
    summary.to_csv(TABLE_DIR / "K4_LD50_summary_smooth.tsv", sep="\t", index=False)
    return summary


def draw_panel(smooth: pd.DataFrame, xmax: int, summary: pd.DataFrame) -> None:
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "axes.linewidth": 0.8,
            "xtick.major.width": 0.8,
            "ytick.major.width": 0.8,
        }
    )

    fig, ax = plt.subplots(figsize=(6.6, 4.4))
    y_max = 0.0
    for group_id, label, sample_n, color in GROUPS:
        curve = smooth[(smooth["group_id"] == group_id) & (smooth["dist_kb"] <= xmax)]
        ld50 = summary.loc[summary["group_id"] == group_id, "LD50_kb_smooth"].iloc[0]
        ax.plot(
            curve["dist_kb"],
            curve["r2_smooth"],
            color=color,
            lw=2.0,
            alpha=0.98,
            solid_capstyle="round",
            label=f"{label} (n={sample_n}, LD50={ld50} kb)",
        )
        y_max = max(y_max, float(curve["r2_smooth"].max()))

    ax.set_xlim(0, xmax)
    ax.set_ylim(0, min(0.55, y_max * 1.12))
    ax.set_xlabel("Distance (kb)", fontsize=11)
    ax.set_ylabel("Mean r$^2$", fontsize=11)
    ax.tick_params(axis="both", labelsize=9, length=3.5)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(axis="y", color="#D9D9D9", linewidth=0.5, alpha=0.55)
    ax.legend(frameon=False, fontsize=8, loc="upper right", handlelength=2.7)
    fig.tight_layout()

    for suffix in ("pdf", "png", "svg"):
        fig.savefig(
            FIG_DIR / f"K4_LD_decay_smooth_0_{xmax}kb.{suffix}",
            dpi=600,
            bbox_inches="tight",
        )
    plt.close(fig)


def main() -> None:
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    TABLE_DIR.mkdir(parents=True, exist_ok=True)

    frames = []
    for group_id, label, sample_n, _color in GROUPS:
        stat_path = STATS_DIR / f"{group_id}.stat.gz"
        frames.append(read_stat(stat_path, group_id, label, sample_n))

    raw = pd.concat(frames, ignore_index=True)
    smooth = build_smooth_table(raw)
    smooth.to_csv(TABLE_DIR / "K4_LD_decay_binned_smooth_points.tsv", sep="\t", index=False)
    summary = write_ld50_table(smooth)

    draw_panel(smooth, xmax=100, summary=summary)
    draw_panel(smooth, xmax=500, summary=summary)

    print(f"Wrote smooth figures to: {FIG_DIR}")
    print(f"Wrote smooth tables to: {TABLE_DIR}")


if __name__ == "__main__":
    main()
