#!/usr/bin/env python3
"""Draw the Figure 4 pan-BGC panel for 33 de-redundant varieties."""

from __future__ import annotations

import importlib.util
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.colors import ListedColormap
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Patch


RUN = Path(
    "${ANALYSIS_DIR}/22_answer_reviews/00_ms/03_V3/04_figure4/"
    "RGA_BGC_meizhou4_33varieties_20260805"
)
SOURCE = Path(
    "${ANALYSIS_DIR}/22_answer_reviews/00_ms/03_V3/04_figure4/"
    "panBGC_visualization_revised_diversity_20260708/scripts/build_revised_panbgc_diversity_figures.py"
)

# Match the submitted Figure 4k encoding: absent blue, one copy yellow, >=2 orange.
COPY_COLORS = ["#3B7FB6", "#E7B65A", "#F28E2B"]


def load_original():
    spec = importlib.util.spec_from_file_location("bgc_original", SOURCE)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot import {SOURCE}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    module.OUT_ROOT = RUN
    module.FIG_DIR = RUN / "figures"
    module.TAB_DIR = RUN / "tables/bgc_plot"
    return module


def label_variety(name: str) -> str:
    return "Meizhou4" if name == "meizhou4" else name


def draw_panel(module, summary, families, presence, row_order, col_order) -> None:
    row_meta = families.set_index("PanBGC_Family").loc[row_order]
    mat = presence.loc[row_order, col_order].astype(int).clip(upper=2)
    n_varieties = len(col_order)

    # Manuscript-composite variant requested on 2026-08-05: preserve width and
    # increase panel height by 2.2-fold (6.4 -> 14.08 inches).
    fig = plt.figure(figsize=(16.8, 14.08), constrained_layout=False)
    gs = GridSpec(
        1, 5, figure=fig,
        width_ratios=[0.14, 0.14, 9.35, 0.55, 2.65],
        left=0.075, right=0.985, top=0.84, bottom=0.255, wspace=0.035,
    )
    ax_class = fig.add_subplot(gs[0, 0])
    ax_type = fig.add_subplot(gs[0, 1])
    ax = fig.add_subplot(gs[0, 2])
    ax_prev = fig.add_subplot(gs[0, 3], sharey=ax)
    ax_lab = fig.add_subplot(gs[0, 4], sharey=ax)

    module.draw_category_strip(ax_class, row_meta, "Family_category", module.CATEGORY_COLORS)
    module.draw_category_strip(ax_type, row_meta, "Major_Type", module.TYPE_COLORS)
    ax_class.set_title("class", fontsize=7, pad=5)
    ax_type.set_title("type", fontsize=7, pad=5)

    ax.imshow(
        mat.values, aspect="auto", interpolation="nearest",
        cmap=ListedColormap(COPY_COLORS), vmin=0, vmax=2,
    )
    ax.set_xticks(np.arange(n_varieties))
    ax.set_xticklabels([label_variety(x) for x in col_order], rotation=58, ha="right", fontsize=5.9)
    ax.set_yticks([])
    ax.tick_params(axis="y", left=False, labelleft=False)
    ax.set_xlabel("Variety", fontsize=9, labelpad=9)
    ax.set_title(
        f"Copy-number diversity of pan-BGC families across {n_varieties} varieties",
        fontsize=12, loc="left", pad=14,
    )
    ax.set_xticks(np.arange(-0.5, n_varieties, 1), minor=True)
    ax.set_yticks(np.arange(-0.5, len(row_order), 1), minor=True)
    ax.grid(which="minor", color="white", linewidth=0.28)
    ax.tick_params(which="minor", bottom=False, left=False)
    module.add_row_separators(ax, row_meta, color="#D7DEE5", linewidth=0.45)

    prevalence = row_meta["Genome_count"].astype(int).values
    bar_colors = [module.CATEGORY_COLORS[value] for value in row_meta["Family_category"]]
    y_pos = np.arange(len(row_order))
    ax_prev.hlines(y_pos, 0, prevalence, color="#CBD3D8", linewidth=1.0, zorder=1)
    ax_prev.scatter(
        prevalence, y_pos, s=40, c=bar_colors, edgecolors="white",
        linewidths=0.55, alpha=0.98, zorder=3, clip_on=False,
    )
    ax_prev.set_xlim(-0.5, n_varieties + 1.2)
    ax_prev.set_xticks([0, n_varieties // 2, n_varieties])
    ax_prev.set_xlabel("varieties", fontsize=7, labelpad=7)
    ax_prev.set_title("prevalence", fontsize=7, pad=5)
    ax_prev.tick_params(axis="y", left=False, labelleft=False)
    ax_prev.tick_params(axis="x", labelsize=7)
    sns.despine(ax=ax_prev, left=True)
    module.add_row_separators(ax_prev, row_meta, color="#E7E7E7", linewidth=0.38)

    ax_lab.set_xlim(0, 1)
    ax_lab.set_xticks([])
    ax_lab.tick_params(axis="y", left=False, labelleft=False)
    for spine in ax_lab.spines.values():
        spine.set_visible(False)
    ax_lab.axvline(0.015, color="#BBBBBB", lw=0.8, ymin=0.03, ymax=0.97)
    private_examples = module.selected_private_examples(row_meta)
    if private_examples:
        title_y = min(row_order.index(family) for family, _ in private_examples) - 2.4
        ax_lab.text(
            0.04, title_y, "Selected private family examples", fontsize=7.8,
            fontweight="bold", ha="left", va="bottom",
        )
    for family, label in private_examples:
        y = row_order.index(family)
        ax_lab.plot([0.015, 0.06], [y, y], color="#999999", lw=0.8, clip_on=False)
        ax_lab.text(0.075, y, label, ha="left", va="center", fontsize=6.4, color="#333333")

    copy_handles = [
        Patch(facecolor=COPY_COLORS[0], edgecolor="none", label="absent"),
        Patch(facecolor=COPY_COLORS[1], edgecolor="none", label="1 BGC"),
        Patch(facecolor=COPY_COLORS[2], edgecolor="none", label=">=2 BGCs"),
    ]
    class_handles = [
        Patch(facecolor=color, edgecolor="none", label=label.replace("_", " "))
        for label, color in module.CATEGORY_COLORS.items()
    ]
    key_types = [
        "saccharide", "cyclopeptide", "fatty_acid", "polyketide",
        "putative", "alkaloid", "fatty_acid-polyketide",
    ]
    type_handles = [
        Patch(facecolor=module.TYPE_COLORS[value], edgecolor="none", label=value.replace("_", " "))
        for value in key_types
    ]
    legends = [
        fig.legend(
            handles=copy_handles, loc="lower left", bbox_to_anchor=(0.075, 0.012),
            ncol=3, frameon=False, fontsize=6.5, title="Copy number",
            title_fontsize=6.8, columnspacing=0.9, handlelength=1.8,
        ),
        fig.legend(
            handles=class_handles, loc="lower left", bbox_to_anchor=(0.265, 0.012),
            ncol=4, frameon=False, fontsize=6.5, title="Family class",
            title_fontsize=6.8, columnspacing=0.9, handlelength=1.8,
        ),
        fig.legend(
            handles=type_handles, loc="lower left", bbox_to_anchor=(0.49, 0.012),
            ncol=7, frameon=False, fontsize=6.5, title="BGC type",
            title_fontsize=6.8, columnspacing=0.9, handlelength=1.8,
        ),
    ]
    for legend in legends:
        legend._legend_box.align = "left"

    fig.suptitle(
        f"Pan-BGC family diversity in the {n_varieties}-variety oil-palm panel",
        fontsize=13, y=0.965,
    )
    fig.text(
        0.075, 0.91,
        f"{n_varieties} varieties | {int(summary['Total_Clusters'].sum())} de-redundant BGCs | "
        f"{len(row_order)} PFAM-fingerprint pan-BGC families",
        ha="left", va="center", fontsize=8, color="#444444",
    )
    module.save_figure(fig, "FigB_panBGC_family_copy_number_diversity_panel_35genomes")


def main() -> None:
    module = load_original()
    module.configure_matplotlib()
    summary = pd.read_csv(RUN / "tables/BGC_summary_33varieties.tsv", sep="\t")
    families = pd.read_csv(RUN / "tables/panbgc_family_summary_33varieties.tsv", sep="\t")
    presence = pd.read_csv(RUN / "tables/panbgc_family_presence_33varieties.tsv", sep="\t", index_col=0)
    # Submitted Figure 4k runs rare/private families at top and core families at bottom.
    row_order = list(reversed(module.order_rows(families)))
    col_order = module.order_columns(summary, presence)
    module.write_order_tables(summary, families, row_order, col_order)
    draw_panel(module, summary, families, presence, row_order, col_order)


if __name__ == "__main__":
    main()
