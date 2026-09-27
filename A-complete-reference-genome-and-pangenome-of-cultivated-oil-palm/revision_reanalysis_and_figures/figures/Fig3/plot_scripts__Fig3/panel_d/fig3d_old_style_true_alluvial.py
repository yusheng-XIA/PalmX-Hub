#!/usr/bin/env python3
"""Render the reviewed shared-orthogroup Figure 3d in the legacy alluvial style."""

from pathlib import Path

import matplotlib as mpl

mpl.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.patches import PathPatch, Patch, Rectangle
from matplotlib.path import Path as MplPath


ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT.parent / "FINAL_ALL_FIGURES_ONE_FOLDER_20260901"
SOURCE = (
    Path(
        "${ANALYSIS_DIR}/"
        "22_answer_reviews/00_ms/03_V3/03_figure3/"
        "04_FL_TN_ASE_reference_panels/source_E_FL_TN_alluvial.tsv"
    )
)

RED = "#F8766D"
BLUE = "#4EA5DF"
TEAL = "#8DD3C7"
YELLOW = "#FFF59D"
DARK = "#333333"
ORDER = ["NoDiff", "HapDom", "Sub", "NoASE"]
CLASS_COLOR = {
    "NoDiff": RED,
    "HapDom": YELLOW,
    "Sub": TEAL,
    "NoASE": BLUE,
}


def intervals(values: pd.Series, gap: float) -> dict[str, tuple[float, float]]:
    usable = 1.0 - gap * (len(ORDER) - 1)
    result: dict[str, tuple[float, float]] = {}
    top = 1.0
    for category in ORDER:
        height = usable * float(values[category]) / float(values.sum())
        result[category] = (top - height, top)
        top -= height + gap
    return result


def main() -> None:
    data = pd.read_csv(SOURCE, sep="\t")
    total = int(data.orthogroups.sum())
    if total != 7059:
        raise AssertionError(f"Expected 7,059 shared orthogroups, found {total:,}")
    expected_pairs = {(left, right) for left in ORDER for right in ORDER}
    observed_pairs = set(zip(data.overall_class_FL, data.overall_class_TN))
    if observed_pairs != expected_pairs:
        raise AssertionError("The complete 4 x 4 transition table is required")

    left = (
        data.groupby("overall_class_FL").orthogroups.sum().reindex(ORDER).fillna(0)
    )
    right = (
        data.groupby("overall_class_TN").orthogroups.sum().reindex(ORDER).fillna(0)
    )
    gap = 0.012
    usable = 1.0 - gap * (len(ORDER) - 1)
    left_intervals = intervals(left, gap)
    right_intervals = intervals(right, gap)
    left_cursor = {category: left_intervals[category][0] for category in ORDER}
    right_cursor = {category: right_intervals[category][0] for category in ORDER}

    mpl.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 7.5,
            "legend.fontsize": 6.4,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )
    fig, ax = plt.subplots(figsize=(4.0, 3.35), facecolor="white")
    fig.subplots_adjust(left=0.13, right=0.94, bottom=0.12, top=0.84)
    fig.text(0.015, 0.985, "d", fontsize=12.5, fontweight="bold", ha="left", va="top")
    ax.set_xlim(0, 1)
    ax.set_ylim(-0.04, 1.03)
    ax.axis("off")

    for row in data.sort_values(["overall_class_FL", "overall_class_TN"]).itertuples():
        source_class = row.overall_class_FL
        target_class = row.overall_class_TN
        height = usable * int(row.orthogroups) / total
        left_bottom = left_cursor[source_class]
        right_bottom = right_cursor[target_class]
        vertices = [
            (0.20, left_bottom),
            (0.45, left_bottom),
            (0.55, right_bottom),
            (0.80, right_bottom),
            (0.80, right_bottom + height),
            (0.55, right_bottom + height),
            (0.45, left_bottom + height),
            (0.20, left_bottom + height),
            (0.20, left_bottom),
        ]
        codes = [
            MplPath.MOVETO,
            MplPath.CURVE4,
            MplPath.CURVE4,
            MplPath.CURVE4,
            MplPath.LINETO,
            MplPath.CURVE4,
            MplPath.CURVE4,
            MplPath.CURVE4,
            MplPath.CLOSEPOLY,
        ]
        ax.add_patch(
            PathPatch(
                MplPath(vertices, codes),
                facecolor=CLASS_COLOR[source_class],
                edgecolor="none",
                alpha=0.40,
            )
        )
        left_cursor[source_class] += height
        right_cursor[target_class] += height

    for x, category_intervals, label in [
        (0.14, left_intervals, "FL"),
        (0.80, right_intervals, "TN"),
    ]:
        for category in ORDER:
            lower, upper = category_intervals[category]
            ax.add_patch(
                Rectangle(
                    (x, lower),
                    0.07,
                    upper - lower,
                    facecolor=CLASS_COLOR[category],
                    edgecolor=DARK,
                    lw=0.65,
                )
            )
        ax.text(x + 0.035, -0.025, label, ha="center", va="top", fontsize=8)

    ax.text(
        0.5,
        -0.025,
        f"Shared ASE-eligible orthogroups, n = {total:,}",
        ha="center",
        va="top",
        fontsize=6.3,
        color="#666666",
    )
    handles = [
        Patch(
            facecolor=CLASS_COLOR[category],
            edgecolor=DARK,
            lw=0.4,
            label=category,
        )
        for category in ORDER
    ]
    fig.legend(
        handles=handles,
        title="Gene type",
        frameon=False,
        ncol=4,
        loc="upper center",
        bbox_to_anchor=(0.57, 0.99),
        columnspacing=0.9,
    )

    stems = [
        "Figure3d_RECOMMENDED_SHARED7059_OLD_STYLE_TRUE_ALLUVIAL",
        "Figure3d_ALTERNATIVE_SHARED7059_OLD_STYLE_TRUE_ALLUVIAL",
    ]
    for stem in stems:
        fig.savefig(
            OUT / f"{stem}.pdf",
            bbox_inches="tight",
            pad_inches=0.045,
            facecolor="white",
        )
        fig.savefig(
            OUT / f"{stem}.png",
            dpi=600,
            bbox_inches="tight",
            pad_inches=0.045,
            facecolor="white",
        )
    plt.close(fig)


if __name__ == "__main__":
    main()
