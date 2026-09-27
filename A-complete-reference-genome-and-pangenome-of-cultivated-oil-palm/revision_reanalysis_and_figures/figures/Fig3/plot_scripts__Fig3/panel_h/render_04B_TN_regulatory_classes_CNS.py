#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
CNS-style redraw of:
04B_TN_regulatory_classes_19stages

Encoding:
    bubble color = genes within stage (%)
    bubble area  = genes within stage (%)

Outputs:
    04B_TN_regulatory_classes_19stages.pdf
    04B_TN_regulatory_classes_19stages.svg
    04B_TN_regulatory_classes_19stages.png
    source_04B_TN_regulatory_classes_19stages.tsv

The script first tries to rebuild the summary table from:
    ../01_ASE/00_shared/runs/RUN-ASE-HETEROSIS-DOWNSTREAM-V2-001/
        output/TN_cis_trans_classified.tsv

If that file is unavailable, it falls back to:
    source_04B_TN_regulatory_classes_19stages.tsv
"""

from pathlib import Path
import math

import matplotlib as mpl

# Safe for SSH / headless server
mpl.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


# ============================================================
# Paths
# ============================================================

OUT = Path(__file__).resolve().parent

SHARED = OUT.parent / "01_ASE" / "00_shared"

ASE_RUN = (
    SHARED
    / "runs"
    / "RUN-ASE-HETEROSIS-DOWNSTREAM-V2-001"
)

RAW_INPUT = (
    ASE_RUN
    / "output"
    / "TN_cis_trans_classified.tsv"
)

SOURCE_TSV = (
    OUT
    / "source_04B_TN_regulatory_classes_19stages.tsv"
)

STEM = "04B_TN_regulatory_classes_19stages"


# ============================================================
# Biological ordering
# ============================================================

STAGES = [
    "0d",
    "15d",
    "35d",
    "50d",
    "65d",
    "80d",
    "95d",
    "110d",
    "125d",
    "140d",
    "155d",
    "170d",
    "185d",
    "12h",
    "24h",
    "36h",
    "48h",
    "60h",
    "72h",
]

REG_ORDER = [
    "I.Cis_only",
    "II.Trans_only",
    "III.Cis_trans_enhancing",
    "IV.Cis_trans_compensating",
    "V.Compensatory",
    "VI.Conserved",
    "VII.Ambiguous",
]

REG_LABELS = [
    "I  Cis only",
    "II  Trans only",
    "III  Cis + trans\nenhancing",
    "IV  Cis + trans\ncompensating",
    "V  Compensatory",
    "VI  Conserved",
    "VII  Ambiguous",
]


# ============================================================
# Global publication style
# ============================================================

def set_style():
    mpl.rcParams.update(
        {
            # Prefer Arial/Helvetica if installed; otherwise DejaVu Sans.
            "font.family": "sans-serif",
            "font.sans-serif": [
                "Arial",
                "Helvetica",
                "Liberation Sans",
                "DejaVu Sans",
            ],

            # CNS-like final-figure typography
            "font.size": 8.0,
            "axes.labelsize": 8.8,
            "axes.titlesize": 10.0,
            "xtick.labelsize": 7.4,
            "ytick.labelsize": 7.7,

            "axes.linewidth": 0.85,

            "xtick.major.width": 0.8,
            "ytick.major.width": 0.8,
            "xtick.major.size": 3.0,
            "ytick.major.size": 3.0,

            # Editable text in PDF/SVG
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",

            "figure.facecolor": "white",
            "savefig.facecolor": "white",

            # Publication raster export
            "savefig.dpi": 600,
        }
    )


def clean_axes(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    ax.spines["left"].set_linewidth(0.85)
    ax.spines["bottom"].set_linewidth(0.85)


def panel_letter(fig, letter="B"):
    fig.text(
        0.012,
        0.985,
        letter,
        ha="left",
        va="top",
        fontsize=16,
        fontweight="bold",
    )


# ============================================================
# Load / reconstruct source table
# ============================================================

def load_regulatory_table():

    if RAW_INPUT.exists():

        print(f"[INFO] Reading raw classification:")
        print(f"       {RAW_INPUT}")

        df = pd.read_csv(
            RAW_INPUT,
            sep="\t",
        )

        required = {
            "stage_index",
            "stage",
            "regulatory_class",
        }

        missing = required - set(df.columns)

        if missing:
            raise ValueError(
                f"Raw input is missing columns: {sorted(missing)}"
            )

        tab = (
            df.groupby(
                [
                    "stage_index",
                    "stage",
                    "regulatory_class",
                ],
                observed=True,
            )
            .size()
            .rename("genes")
            .reset_index()
        )

        totals = (
            df.groupby(
                [
                    "stage_index",
                    "stage",
                ],
                observed=True,
            )
            .size()
            .rename("stage_total")
            .reset_index()
        )

        tab = tab.merge(
            totals,
            on=[
                "stage_index",
                "stage",
            ],
            how="left",
        )

        tab["percentage"] = (
            100.0
            * tab["genes"]
            / tab["stage_total"]
        )

        tab["inference_level"] = (
            "exploratory_parent_n1"
        )

        tab.to_csv(
            SOURCE_TSV,
            sep="\t",
            index=False,
        )

        return tab, RAW_INPUT

    # --------------------------------------------------------
    # Fallback to already generated summary TSV
    # --------------------------------------------------------

    if SOURCE_TSV.exists():

        print(
            "[INFO] Raw classification file not found; "
            "using existing source TSV:"
        )
        print(f"       {SOURCE_TSV}")

        tab = pd.read_csv(
            SOURCE_TSV,
            sep="\t",
        )

        required = {
            "stage",
            "regulatory_class",
            "percentage",
        }

        missing = required - set(tab.columns)

        if missing:
            raise ValueError(
                f"Source TSV missing columns: {sorted(missing)}"
            )

        return tab, SOURCE_TSV

    raise FileNotFoundError(
        "\nNeither input file was found:\n"
        f"  {RAW_INPUT}\n"
        f"  {SOURCE_TSV}\n"
    )


# ============================================================
# Main plot
# ============================================================

def draw_regulatory_bubble():

    tab, data_source = load_regulatory_table()

    # --------------------------------------------------------
    # Matrix: regulatory classes × stages
    # --------------------------------------------------------

    piv = (
        tab.pivot(
            index="regulatory_class",
            columns="stage",
            values="percentage",
        )
        .reindex(
            index=REG_ORDER,
            columns=STAGES,
        )
        .fillna(0.0)
    )

    matrix = piv.to_numpy(dtype=float)

    maxp = float(
        np.nanmax(matrix)
    )

    # Use a clean upper bound rather than ending colorbar at
    # an awkward value such as 53.1193%.
    #
    # Your current dataset max ~53%, therefore display vmax=55.
    display_max = max(
        55.0,
        5.0 * math.ceil(maxp / 5.0),
    )

    print(
        f"[INFO] Maximum observed percentage = "
        f"{maxp:.3f}%"
    )

    print(
        f"[INFO] Plotting scale maximum      = "
        f"{display_max:.1f}%"
    )

    # --------------------------------------------------------
    # Bubble area transformation
    #
    # IMPORTANT:
    # matplotlib scatter(s=...) uses AREA in points^2.
    #
    # Both the main bubbles and bubble-size legend call
    # this exact same function.
    # --------------------------------------------------------

    def bubble_size(pct):

        arr = np.asarray(
            pct,
            dtype=float,
        )

        arr = np.clip(
            arr,
            0.0,
            display_max,
        )

        # Minimum visible area + proportional component.
        #
        # Approximate range:
        # 0%  -> 16 pt²
        # 10% -> ~91 pt²
        # 30% -> ~240 pt²
        # 50% -> ~389 pt²
        #
        # This preserves the original figure's visual hierarchy
        # while leaving adequate separation between bubbles.

        sizes = (
            16.0
            + 410.0
            * arr
            / display_max
        )

        if np.ndim(sizes) == 0:
            return float(sizes)

        return sizes

    # ========================================================
    # Canvas
    # ========================================================

    fig, ax = plt.subplots(
        figsize=(8.65, 3.85)
    )

    # Leave dedicated room at right for legends.
    fig.subplots_adjust(
        left=0.205,
        right=0.795,
        bottom=0.255,
        top=0.835,
    )

    # ========================================================
    # Build bubble coordinates
    # ========================================================

    xs = []
    ys = []
    vals = []

    for yi, reg in enumerate(
        REG_ORDER[::-1]
    ):
        for xi, stage in enumerate(
            STAGES
        ):

            xs.append(xi)
            ys.append(yi)

            vals.append(
                float(
                    piv.loc[
                        reg,
                        stage,
                    ]
                )
            )

    xs = np.asarray(xs)
    ys = np.asarray(ys)
    vals = np.asarray(vals)

    # ========================================================
    # Bubble plot
    # ========================================================

    sc = ax.scatter(
        xs,
        ys,

        # Bubble AREA
        s=bubble_size(vals),

        # Bubble COLOR
        c=vals,

        cmap="Reds",

        vmin=0.0,
        vmax=display_max,

        # Thin white separation makes dense bubbles cleaner.
        edgecolors="white",
        linewidths=0.45,

        alpha=0.98,

        zorder=3,
    )

    # ========================================================
    # X axis
    # ========================================================

    ax.set_xticks(
        np.arange(
            len(STAGES)
        )
    )

    ax.set_xticklabels(
        STAGES,
        rotation=45,
        ha="right",
        rotation_mode="anchor",
        fontsize=7.4,
    )

    ax.tick_params(
        axis="x",
        pad=2.5,
    )

    ax.set_xlabel(
        "Stage",
        fontsize=8.8,
        labelpad=5,
    )

    # ========================================================
    # Y axis
    # ========================================================

    ax.set_yticks(
        np.arange(
            len(REG_ORDER)
        )
    )

    ax.set_yticklabels(
        REG_LABELS[::-1],
        fontsize=7.8,
    )

    ax.tick_params(
        axis="y",
        pad=4,
    )

    ax.set_ylabel(
        "Regulatory class",
        fontsize=8.8,
        labelpad=10,
    )

    # ========================================================
    # Plot limits
    # ========================================================

    ax.set_xlim(
        -0.7,
        18.7,
    )

    ax.set_ylim(
        -0.65,
        6.65,
    )

    # ========================================================
    # Development vs postharvest separator
    # ========================================================

    # Between 185d (index 12) and 12h (index 13)
    ax.axvline(
        12.5,
        color="#777777",
        linewidth=0.80,
        linestyle="--",
        zorder=1,
    )

    # Phase headings
    ax.text(
        6.0,
        6.49,
        "Seed development",
        ha="center",
        va="bottom",
        fontsize=7.8,
        color="#3D7A69",
        fontweight="medium",
    )

    ax.text(
        15.5,
        6.49,
        "Postharvest",
        ha="center",
        va="bottom",
        fontsize=7.8,
        color="#B68A25",
        fontweight="medium",
    )

    clean_axes(ax)

    # ========================================================
    # Main title
    # ========================================================

    ax.set_title(
        "TN regulatory-class landscape across 19 stages",
        fontsize=10.0,
        fontweight="bold",
        pad=11,
    )

    # Panel B
    panel_letter(
        fig,
        "B",
    )

    # ========================================================
    # Colorbar
    # ========================================================

    # Dedicated colorbar axis prevents resize / crowding.
    cax = fig.add_axes(
        [
            0.830,  # left
            0.320,  # bottom
            0.017,  # width
            0.500,  # height
        ]
    )

    cbar = fig.colorbar(
        sc,
        cax=cax,
    )

    cbar.set_label(
        "Genes within stage (%)",
        fontsize=8.0,
        labelpad=7,
    )

    # Publication-friendly ticks.
    tick_values = [
        x
        for x in
        [0, 10, 20, 30, 40, 50]
        if x <= display_max
    ]

    cbar.set_ticks(
        tick_values
    )

    cbar.ax.tick_params(
        labelsize=7.1,
        width=0.75,
        length=3.0,
        pad=2.5,
    )

    cbar.outline.set_linewidth(
        0.8
    )

    # ========================================================
    # Bubble-size legend
    # ========================================================

    # Create a dedicated miniature legend axis.
    # This gives substantially more stable publication layout
    # than an automatically placed Matplotlib legend.

    lax = fig.add_axes(
        [
            0.805,  # left
            0.075,  # bottom
            0.175,  # width
            0.175,  # height
        ]
    )

    lax.set_xlim(
        0,
        1,
    )

    lax.set_ylim(
        0,
        1,
    )

    lax.axis(
        "off"
    )

    lax.text(
        0.5,
        0.94,
        "Bubble size (%)",
        ha="center",
        va="top",
        fontsize=7.7,
        fontweight="medium",
        color="#222222",
    )

    size_levels = [
        10,
        30,
        50,
    ]

    legend_x = [
        0.18,
        0.50,
        0.82,
    ]

    legend_y = 0.55

    # IMPORTANT:
    # Same bubble_size() mapping as the main plot.
    lax.scatter(
        legend_x,
        [legend_y] * 3,
        s=[
            bubble_size(v)
            for v in size_levels
        ],
        facecolor="#D9D9D9",
        edgecolor="#666666",
        linewidth=0.55,
        zorder=3,
    )

    for x, value in zip(
        legend_x,
        size_levels,
    ):

        lax.text(
            x,
            0.16,
            f"{value}%",
            ha="center",
            va="center",
            fontsize=7.0,
            color="#333333",
        )

    # ========================================================
    # Exploratory-analysis note
    # ========================================================

    fig.text(
        0.795,
        0.060,
        "Exploratory: TK/NS parental expression n = 1 per stage",
        ha="right",
        va="bottom",
        fontsize=6.35,
        color="#555555",
    )

    # ========================================================
    # Save
    # ========================================================

    pdf_path = OUT / f"{STEM}.pdf"
    svg_path = OUT / f"{STEM}.svg"
    png_path = OUT / f"{STEM}.png"

    fig.savefig(
        pdf_path,
        bbox_inches="tight",
    )

    fig.savefig(
        svg_path,
        bbox_inches="tight",
    )

    fig.savefig(
        png_path,
        dpi=600,
        bbox_inches="tight",
    )

    plt.close(
        fig
    )

    print()
    print("[DONE] 04B regenerated")
    print(f"       Data source: {data_source}")
    print(f"       PDF: {pdf_path}")
    print(f"       SVG: {svg_path}")
    print(f"       PNG: {png_path}")
    print(f"       TSV: {SOURCE_TSV}")
    print()


# ============================================================
# Run
# ============================================================

if __name__ == "__main__":

    set_style()

    draw_regulatory_bubble()

