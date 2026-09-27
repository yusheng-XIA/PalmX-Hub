#!/usr/bin/env python3
"""Redesign RGA tree heatmaps as a tree-integrated dot matrix.

The figure combines the two earlier heatmaps:
- dot size: raw RGA count;
- dot color: column-wise z-score of RGA abundance.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from Bio import Phylo
from matplotlib.colors import LinearSegmentedColormap, Normalize, to_rgba
from matplotlib.cm import ScalarMappable
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Patch, PathPatch
from matplotlib.path import Path as MplPath


RGA_DIR = Path(
    "${ANALYSIS_DIR}/21_MS/05_result/01_RGA"
)
TREE_FILE = Path(
    "${ANALYSIS_DIR}/14_pan_genome/02_geneClusterPAV_contigs/data/OrthoFinder/Results_Jan01/Species_Tree/SpeciesTree_rooted.txt"
)
SYNTENY_COUNTS_FILE = Path(
    "${DATA_DIR3}/data_backup/for_others/user/RGA/synteny/34spe_joined_RGA_withtandem.i1.blocks_sorted.edit.edit_counts.merge"
)
OUT_ROOT = Path(
    "${ANALYSIS_DIR}/22_answer_reviews/00_ms/03_V3/04_figure4/RGA_tree_dotmatrix_redesign_20260708"
)
FIG_DIR = OUT_ROOT / "figures"
TAB_DIR = OUT_ROOT / "tables"


GROUP_COLORS = {
    "African/American haplotypes": "#62B6B7",
    "BK haplotypes": "#9C8AC7",
    "Reference genomes": "#D89C35",
    "EG assemblies": "#86B68A",
}

RGA_COLUMNS = ["NBS-class", "RLK", "RLP", "TM-CC", "Other"]
RGA_DISPLAY_LABELS = {
    "NBS-class": "NBS",
    "RLK": "RLK",
    "RLP": "RLP",
    "TM-CC": "TM-CC",
    "Other": "Other",
}

ZSCORE_CMAP = LinearSegmentedColormap.from_list(
    "soft_demo_blue_white_red",
    ["#4B5AA7", "#A7AED3", "#F7F7F7", "#E8A39E", "#C84B4B"],
)


def configure_matplotlib() -> None:
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "axes.linewidth": 0.75,
            "axes.edgecolor": "#333333",
            "xtick.major.width": 0.7,
            "ytick.major.width": 0.7,
        }
    )


def tree_tip_to_summary_name(tip: str) -> str:
    return tip.replace("_pa_genome", "")


def display_genome_name(name: str) -> str:
    if name.isdigit():
        return f"EG_{int(name):03d}"
    display_map = {
        "dura": "EG_dura",
        "pisifera": "EG_pisifera",
        "houke": "EG_houke",
        "niriliya": "EG_niriliya",
    }
    return display_map.get(name, name)


def genome_group(name: str) -> str:
    if name in {"Africa_hap2", "American_hap1"}:
        return "African/American haplotypes"
    if name in {"bk_hap1", "bk_hap2"}:
        return "BK haplotypes"
    if name in {"dura", "pisifera", "houke", "niriliya"}:
        return "Reference genomes"
    return "EG assemblies"


def load_rga_summary() -> pd.DataFrame:
    df = pd.read_csv(RGA_DIR / "RGA_summary_all_genomes.csv")
    df["NBS-class"] = df[["NBS", "CNL", "TNL", "CN", "TN", "NL", "TX"]].sum(axis=1)
    df["Other"] = df["Others"] if "Others" in df.columns else 0
    keep = ["Genome"] + RGA_COLUMNS
    return df[keep].copy()


def compute_tree_layout(tree: Phylo.BaseTree.Tree) -> tuple[dict[object, float], dict[object, float], list[object]]:
    tree.ladderize(reverse=True)
    terminals = tree.get_terminals()
    n = len(terminals)
    y_by_clade: dict[object, float] = {}
    for i, tip in enumerate(terminals):
        y_by_clade[tip] = n - 1 - i

    def assign_internal_y(clade: object) -> float:
        if clade in y_by_clade:
            return y_by_clade[clade]
        child_y = [assign_internal_y(child) for child in clade.clades]
        y_by_clade[clade] = float(np.mean(child_y))
        return y_by_clade[clade]

    assign_internal_y(tree.root)

    x_by_clade: dict[object, float] = {}

    def assign_x(clade: object, x_parent: float = 0.0) -> None:
        branch = 0.0 if clade is tree.root else 1.0
        x_here = x_parent + branch
        x_by_clade[clade] = x_here
        for child in clade.clades:
            assign_x(child, x_here)

    assign_x(tree.root, 0.0)
    terminal_x = max(x_by_clade[tip] for tip in terminals)
    for tip in terminals:
        x_by_clade[tip] = terminal_x
    return x_by_clade, y_by_clade, terminals


def draw_tree(ax: plt.Axes, tree: Phylo.BaseTree.Tree, x_by_clade: dict[object, float], y_by_clade: dict[object, float]) -> None:
    max_x = max(x_by_clade.values())
    tree_color = "#3F3F3F"
    tree_lw = 0.76

    def draw_rounded_branch(x0: float, y0: float, x1: float, y1: float) -> None:
        dx = max(x1 - x0, 0.001)
        dy = y1 - y0
        if abs(dy) < 1e-6:
            vertices = [(x0, y0), (x1, y1)]
            codes = [MplPath.MOVETO, MplPath.LINETO]
        else:
            direction = 1.0 if dy > 0 else -1.0
            radius_x = min(0.28, max(dx * 0.22, 0.12))
            radius_y = min(abs(dy) * 0.48, 0.48)
            y_before_corner = y1 - direction * radius_y
            x_after_corner = x0 + radius_x
            vertices = [
                (x0, y0),
                (x0, y_before_corner),
                (x0, y1),
                (x_after_corner, y1),
                (x1, y1),
            ]
            codes = [
                MplPath.MOVETO,
                MplPath.LINETO,
                MplPath.CURVE3,
                MplPath.CURVE3,
                MplPath.LINETO,
            ]
        patch = PathPatch(
            MplPath(vertices, codes),
            facecolor="none",
            edgecolor=tree_color,
            lw=tree_lw,
            capstyle="round",
            joinstyle="round",
            antialiased=True,
        )
        ax.add_patch(patch)

    def recurse(clade: object) -> None:
        if not clade.clades:
            return
        for child in clade.clades:
            draw_rounded_branch(
                x_by_clade[clade],
                y_by_clade[clade],
                x_by_clade[child],
                y_by_clade[child],
            )
            recurse(child)

    recurse(tree.root)
    for tip in tree.get_terminals():
        y = y_by_clade[tip]
        ax.plot(
            [x_by_clade[tip], max_x * 1.015],
            [y, y],
            color="#D7DADF",
            lw=0.55,
            ls=(0, (1.5, 2.2)),
        )

    ax.set_xlim(-max_x * 0.03, max_x * 1.04)
    ax.set_ylim(-0.5, len(tree.get_terminals()) - 0.5)
    ax.axis("off")


def dot_size(count: float, max_count: float) -> float:
    if count <= 0:
        return 0.0
    return 10.0 + (count / max_count) ** 0.72 * 250.0


def build_plot_data(summary: pd.DataFrame, terminals: list[object]) -> pd.DataFrame:
    summary_i = summary.set_index("Genome")
    rows = []
    for tip in terminals:
        raw_name = tree_tip_to_summary_name(tip.name)
        if raw_name not in summary_i.index:
            raise ValueError(f"No RGA summary row for tree tip {tip.name} -> {raw_name}")
        for col in RGA_COLUMNS:
            rows.append(
                {
                    "Tree_tip": tip.name,
                    "Genome": raw_name,
                    "Display": display_genome_name(raw_name),
                    "Genome_group": genome_group(raw_name),
                    "RGA_type": col,
                    "Count": float(summary_i.loc[raw_name, col]),
                }
            )
    plot_df = pd.DataFrame(rows)
    plot_df["Zscore"] = plot_df.groupby("RGA_type")["Count"].transform(
        lambda x: (x - x.mean()) / x.std(ddof=0) if x.std(ddof=0) > 0 else 0
    )
    plot_df["Zscore_clipped"] = plot_df["Zscore"].clip(-3, 3)
    return plot_df


def draw_group_strip(ax: plt.Axes, terminals: list[object]) -> None:
    colors = []
    for tip in terminals:
        raw_name = tree_tip_to_summary_name(tip.name)
        colors.append(to_rgba(GROUP_COLORS[genome_group(raw_name)]))
    rgba = np.array([[c] for c in colors])
    ax.imshow(rgba, aspect="auto", interpolation="nearest")
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)


def draw_labels(ax: plt.Axes, terminals: list[object], y_by_clade: dict[object, float]) -> None:
    ax.set_xlim(0, 1)
    ax.set_ylim(-0.5, len(terminals) - 0.5)
    ax.axis("off")
    for tip in terminals:
        y = y_by_clade[tip]
        raw_name = tree_tip_to_summary_name(tip.name)
        ax.text(
            0.98,
            y,
            display_genome_name(raw_name),
            ha="right",
            va="center",
            fontsize=7.2,
            color="#222222",
        )


def draw_dot_matrix(ax: plt.Axes, plot_df: pd.DataFrame, terminals: list[object], y_by_clade: dict[object, float]) -> None:
    x_step = 0.66
    x_map = {col: i * x_step for i, col in enumerate(RGA_COLUMNS)}
    y_map = {tree_tip_to_summary_name(tip.name): y_by_clade[tip] for tip in terminals}
    max_count = max(plot_df["Count"].max(), 1.0)
    cmap = ZSCORE_CMAP
    norm = Normalize(vmin=-3, vmax=3)

    for _, row in plot_df.iterrows():
        if row["Count"] <= 0:
            continue
        ax.scatter(
            x_map[row["RGA_type"]],
            y_map[row["Genome"]],
            s=dot_size(row["Count"], max_count),
            c=[cmap(norm(row["Zscore_clipped"]))],
            edgecolors="white",
            linewidths=0.45,
            alpha=0.88,
        )

    ax.set_xlim(-0.36, x_map[RGA_COLUMNS[-1]] + 0.36)
    ax.set_ylim(-0.5, len(terminals) - 0.5)
    ax.set_xticks([x_map[x] for x in RGA_COLUMNS])
    ax.set_xticklabels([RGA_DISPLAY_LABELS[x] for x in RGA_COLUMNS], fontsize=10, fontweight="bold")
    ax.xaxis.tick_top()
    ax.tick_params(axis="x", length=0, pad=7)
    ax.tick_params(axis="y", left=False, labelleft=False)
    for spine in ax.spines.values():
        spine.set_visible(False)


def compute_private_rga_counts() -> pd.DataFrame:
    if not SYNTENY_COUNTS_FILE.exists():
        return pd.DataFrame(columns=["Raw_genome", "Private_synteny_blocks", "Private_RGA_copies"])
    counts = pd.read_csv(SYNTENY_COUNTS_FILE, sep="\t")
    counts = counts.apply(pd.to_numeric, errors="coerce").fillna(0).astype(int)
    present = counts > 0
    private_rows = present.sum(axis=1) == 1
    rows = []
    for genome in counts.columns:
        mask = private_rows & present[genome]
        rows.append(
            {
                "Genome": display_genome_name(genome),
                "Raw_genome": genome,
                "Private_synteny_blocks": int(mask.sum()),
                "Private_RGA_copies": int(counts.loc[mask, genome].sum()),
            }
        )
    return pd.DataFrame(rows)


def private_dot_size(count: float, max_count: float) -> float:
    if count <= 0:
        return 0.0
    return 9.0 + (count / max_count) ** 0.68 * 150.0


def draw_private_rga_axis(
    ax: plt.Axes,
    terminals: list[object],
    y_by_clade: dict[object, float],
    private_counts: pd.DataFrame,
) -> None:
    ax.set_xlim(0, 1)
    ax.set_ylim(-0.5, len(terminals) - 0.5)
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)

    ax.text(
        0.5,
        len(terminals) - 0.08,
        "Private\nRGA",
        ha="center",
        va="bottom",
        fontsize=8.8,
        fontweight="bold",
        color="#222222",
        linespacing=0.9,
    )
    count_map = private_counts.set_index("Raw_genome")["Private_RGA_copies"].to_dict()
    max_count = max(count_map.values()) if count_map else 1
    for tip in terminals:
        raw_name = tree_tip_to_summary_name(tip.name)
        private_count = int(count_map.get(raw_name, 0))
        if private_count <= 0:
            continue
        y = y_by_clade[tip]
        ax.scatter(
            0.46,
            y,
            s=private_dot_size(private_count, max_count),
            c="#D86C59",
            edgecolors="white",
            linewidths=0.55,
            alpha=0.9,
            zorder=3,
        )
        if private_count <= 12 or private_count >= 100:
            ax.text(
                0.73,
                y,
                str(private_count),
                ha="left",
                va="center",
                fontsize=6.8,
                color="#9D3E34",
            )


def add_legends(fig: plt.Figure, plot_df: pd.DataFrame, private_counts: pd.DataFrame) -> None:
    group_handles = [
        Patch(facecolor=color, edgecolor="none", label=label)
        for label, color in GROUP_COLORS.items()
    ]
    leg1 = fig.legend(
        handles=group_handles,
        loc="lower left",
        bbox_to_anchor=(0.07, 0.04),
        ncol=2,
        frameon=False,
        title="Genome group",
        fontsize=8,
        title_fontsize=8.5,
        handlelength=1.1,
        columnspacing=1.35,
    )
    leg1._legend_box.align = "left"

    max_count = max(plot_df["Count"].max(), 1.0)
    size_values = [50, 150, 300, 600]
    size_values = [v for v in size_values if v <= max_count]
    handles = [
        plt.Line2D(
            [],
            [],
            marker="o",
            linestyle="",
            markerfacecolor="white",
            markeredgecolor="#333333",
            markersize=np.sqrt(dot_size(v, max_count)),
            label=str(v),
        )
        for v in size_values
    ]
    leg2 = fig.legend(
        handles=handles,
        loc="lower left",
        bbox_to_anchor=(0.41, 0.047),
        ncol=len(handles),
        frameon=False,
        title="RGA count",
        fontsize=8,
        title_fontsize=8.5,
        handletextpad=0.9,
        columnspacing=1.15,
    )
    leg2._legend_box.align = "left"

    cax = fig.add_axes([0.69, 0.055, 0.21, 0.018])
    cb = fig.colorbar(
        ScalarMappable(norm=Normalize(vmin=-3, vmax=3), cmap=ZSCORE_CMAP),
        cax=cax,
        orientation="horizontal",
    )
    cb.set_ticks([-3, -1.5, 0, 1.5, 3])
    cb.ax.tick_params(labelsize=7, length=2)
    cb.outline.set_linewidth(0.5)
    fig.text(0.69, 0.083, "Normalized abundance (column z-score)", fontsize=8.5, ha="left")

    max_private = max(private_counts["Private_RGA_copies"].max(), 1) if not private_counts.empty else 1
    private_values = [10, 50, 100]
    private_values = [v for v in private_values if v <= max_private]
    if private_values:
        handles = [
            plt.Line2D(
                [],
                [],
                marker="o",
                linestyle="",
                markerfacecolor="#D86C59",
                markeredgecolor="white",
                markersize=np.sqrt(private_dot_size(v, max_private)),
                label=str(v),
            )
            for v in private_values
        ]
        leg3 = fig.legend(
            handles=handles,
            loc="lower right",
            bbox_to_anchor=(0.965, 0.088),
            ncol=len(handles),
            frameon=False,
            title="Private RGA copies",
            fontsize=7.2,
            title_fontsize=7.7,
            handletextpad=0.6,
            columnspacing=0.85,
        )
        leg3._legend_box.align = "left"


def save_figure(fig: plt.Figure, stem: str) -> None:
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    for ext in ("pdf", "svg", "png"):
        fig.savefig(FIG_DIR / f"{stem}.{ext}", dpi=600 if ext == "png" else None, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    configure_matplotlib()
    summary = load_rga_summary()
    tree = Phylo.read(str(TREE_FILE), "newick")
    x_by_clade, y_by_clade, terminals = compute_tree_layout(tree)
    plot_df = build_plot_data(summary, terminals)
    private_counts = compute_private_rga_counts()

    TAB_DIR.mkdir(parents=True, exist_ok=True)
    plot_df.to_csv(TAB_DIR / "RGA_major_group_counts_and_zscores.tsv", sep="\t", index=False)
    pd.DataFrame(
        {
            "Tree_tip": [tip.name for tip in terminals],
            "Genome": [tree_tip_to_summary_name(tip.name) for tip in terminals],
            "Display": [display_genome_name(tree_tip_to_summary_name(tip.name)) for tip in terminals],
            "Genome_group": [genome_group(tree_tip_to_summary_name(tip.name)) for tip in terminals],
            "Y_order_top_to_bottom": range(1, len(terminals) + 1),
        }
    ).to_csv(TAB_DIR / "RGA_tree_dotmatrix_tip_order.tsv", sep="\t", index=False)
    private_counts.sort_values(
        ["Private_RGA_copies", "Private_synteny_blocks", "Genome"],
        ascending=[False, False, True],
    ).to_csv(TAB_DIR / "RGA_global_private_synteny_block_audit.tsv", sep="\t", index=False)

    fig = plt.figure(figsize=(10.2, 8.35), constrained_layout=False)
    gs = GridSpec(
        1,
        5,
        figure=fig,
        width_ratios=[2.05, 0.095, 0.58, 2.28, 0.28],
        left=0.055,
        right=0.965,
        top=0.84,
        bottom=0.18,
        wspace=0.006,
    )
    ax_tree = fig.add_subplot(gs[0, 0])
    ax_group = fig.add_subplot(gs[0, 1])
    ax_labels = fig.add_subplot(gs[0, 2])
    ax_matrix = fig.add_subplot(gs[0, 3])
    ax_private = fig.add_subplot(gs[0, 4], sharey=ax_matrix)

    draw_tree(ax_tree, tree, x_by_clade, y_by_clade)
    draw_group_strip(ax_group, terminals)
    draw_labels(ax_labels, terminals, y_by_clade)
    draw_dot_matrix(ax_matrix, plot_df, terminals, y_by_clade)
    draw_private_rga_axis(ax_private, terminals, y_by_clade, private_counts)

    fig.text(0.025, 0.955, "a", fontsize=18, fontweight="bold", ha="left", va="top")
    fig.text(
        0.5,
        0.955,
        "RGA repertoire across oil-palm genomes",
        fontsize=14,
        fontweight="bold",
        ha="center",
        va="top",
    )
    fig.text(
        0.5,
        0.925,
        "Species-tree order with RGA class counts and column-normalized abundance",
        fontsize=9,
        color="#555555",
        ha="center",
        va="top",
    )
    fig.text(
        0.055,
        0.155,
        "Left tree is drawn as a cladogram with aligned terminal tips; branch lengths are not used.",
        fontsize=6.6,
        color="#666666",
        ha="left",
        va="center",
    )
    add_legends(fig, plot_df, private_counts)
    save_figure(fig, "RGA_tree_dotmatrix_major_groups")


if __name__ == "__main__":
    main()
