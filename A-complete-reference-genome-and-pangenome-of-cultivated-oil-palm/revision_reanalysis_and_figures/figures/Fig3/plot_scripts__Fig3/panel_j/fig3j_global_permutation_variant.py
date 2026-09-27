#!/usr/bin/env python3
"""Render the reviewed Figure 3j with the corrected global permutation P value."""

from pathlib import Path
import math

import matplotlib as mpl

mpl.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm


ROOT = Path(__file__).resolve().parent.parent
OUT = (
    ROOT.parent
    / "FINAL_ALL_FIGURES_ONE_FOLDER_20260901"
)
PLOTDATA = ROOT / "Figure3j_plotdata_gene_clustered.tsv"
TESTS = ROOT / "Figure3j_gene_trajectory_permutation_tests.tsv"

INHERITANCE_ORDER = ["PDO", "DO", "ODO"]
REG_ORDER = [
    "I.Cis_only",
    "II.Trans_only",
    "III.Cis_trans_enhancing",
    "IV.Cis_trans_compensating",
    "V.Compensatory",
    "VI.Conserved",
    "VII.Ambiguous",
]
REG_SHORT = [
    "I  Cis",
    "II  Trans",
    "III  Enh.",
    "IV  Cis+trans comp.",
    "V  Compensatory",
    "VI  Conserved",
    "VII  Ambiguous",
]
CMAP = LinearSegmentedColormap.from_list(
    "dpm_blue_gold_red",
    ["#3767A6", "#8FAED0", "#F7F5EE", "#D8B64C", "#E56565"],
)


def matrix(data: pd.DataFrame, value: str) -> np.ndarray:
    return (
        data.pivot(
            index="inheritance_class",
            columns="regulatory_class",
            values=value,
        )
        .reindex(index=INHERITANCE_ORDER, columns=REG_ORDER)
        .to_numpy(float)
    )


def render(
    percentage: np.ndarray,
    residual: np.ndarray,
    global_p: float,
    *,
    show_cell_numbers: bool,
    stem: str,
) -> None:
    vmax = max(75.0, 25.0 * math.ceil(float(np.nanmax(np.abs(residual))) / 25.0))
    norm = TwoSlopeNorm(vmin=-vmax, vcenter=0.0, vmax=vmax)

    mpl.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Liberation Sans", "DejaVu Sans"],
            "font.size": 8.2,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "savefig.facecolor": "white",
        }
    )
    fig = plt.figure(figsize=(3.30, 2.90), facecolor="white")
    ax = fig.add_axes([0.13, 0.27, 0.80, 0.53])
    image = ax.pcolormesh(
        np.arange(len(REG_ORDER) + 1) - 0.5,
        np.arange(len(INHERITANCE_ORDER) + 1) - 0.5,
        residual,
        cmap=CMAP,
        norm=norm,
        shading="flat",
        edgecolors="none",
        rasterized=False,
    )
    ax.set_xlim(-0.5, len(REG_ORDER) - 0.5)
    ax.set_ylim(len(INHERITANCE_ORDER) - 0.5, -0.5)
    if show_cell_numbers:
        for row in range(len(INHERITANCE_ORDER)):
            for column in range(len(REG_ORDER)):
                ax.text(
                    column,
                    row,
                    f"{percentage[row, column]:.1f}%",
                    ha="center",
                    va="center",
                    fontsize=5.1,
                    color=(
                        "white"
                        if abs(residual[row, column]) > 0.55 * vmax
                        else "#222222"
                    ),
                )
    ax.set_xticks(range(len(REG_ORDER)))
    ax.set_xticklabels(
        REG_SHORT,
        rotation=37,
        ha="right",
        rotation_mode="anchor",
        fontsize=5.3,
    )
    ax.set_yticks(range(len(INHERITANCE_ORDER)), INHERITANCE_ORDER)
    ax.tick_params(axis="y", labelsize=6.8)
    for spine in ax.spines.values():
        spine.set_visible(False)

    color_axis = fig.add_axes([0.13, 0.87, 0.80, 0.025])
    colorbar = fig.colorbar(image, cax=color_axis, orientation="horizontal")
    ticks = np.linspace(-vmax, vmax, 7)
    colorbar.set_ticks(ticks)
    colorbar.set_ticklabels([f"{tick:.0f}" for tick in ticks])
    colorbar.ax.xaxis.set_ticks_position("top")
    colorbar.ax.tick_params(labelsize=5.4, pad=1.2, length=2.0)
    colorbar.outline.set_linewidth(0.75)
    fig.text(0.93, 0.985, "Pearson residual", ha="right", va="top", fontsize=5.8)
    fig.text(
        0.53,
        0.835,
        f"Global gene-trajectory permutation P = {global_p:.4f}",
        ha="center",
        va="center",
        fontsize=4.8,
        color="#27323A",
    )
    fig.text(0.012, 0.985, "j", ha="left", va="top", fontsize=12.5, fontweight="bold")

    fig.savefig(OUT / f"{stem}.pdf", bbox_inches="tight", pad_inches=0.045)
    fig.savefig(
        OUT / f"{stem}.png",
        dpi=600,
        bbox_inches="tight",
        pad_inches=0.045,
    )
    plt.close(fig)


def main() -> None:
    data = pd.read_csv(PLOTDATA, sep="\t")
    tests = pd.read_csv(TESTS, sep="\t")
    global_rows = tests.loc[tests.scope.eq("global_3x7")]
    if len(global_rows) != 1:
        raise ValueError("Expected exactly one global 3x7 permutation result")
    global_p = float(global_rows.iloc[0].empirical_p)
    if not np.isclose(global_p, 0.0005):
        raise AssertionError(f"Unexpected global permutation P: {global_p}")

    cell_tests = tests.loc[tests.scope.eq("cell_two_sided_residual")]
    if len(cell_tests) != 21 or not np.allclose(cell_tests.holm_p, 0.0105):
        raise AssertionError("The 21 corrected cell-level Holm P values changed")

    percentage = matrix(data, "row_percentage")
    residual = matrix(data, "trajectory_permutation_residual")
    render(
        percentage,
        residual,
        global_p,
        show_cell_numbers=True,
        stem="Figure3j_revised_gene_clustered_WITH_GLOBAL_PERMUTATION_P",
    )
    render(
        percentage,
        residual,
        global_p,
        show_cell_numbers=False,
        stem="Figure3j_revised_gene_clustered_NO_CELL_NUMBERS_WITH_GLOBAL_PERMUTATION_P",
    )


if __name__ == "__main__":
    main()
