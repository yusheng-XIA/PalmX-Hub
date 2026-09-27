#!/usr/bin/env python3
"""Draw three requested publication-style population genetics PDFs.

Outputs:
1. K=3/4/6 tree + STRUCTURE panel, with tree colored by K=4 groups.
2. K=4 PCA PC1/PC2 panel.
3. K=4 pi-FST network panel.
"""

from __future__ import annotations

from collections import Counter
import math
from pathlib import Path

from Bio import Phylo
import matplotlib as mpl

mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.lines import Line2D
from matplotlib.patches import Ellipse, Patch, PathPatch, Rectangle
from matplotlib.path import Path as MplPath
from mpl_toolkits.mplot3d.art3d import Line3DCollection
import numpy as np
import pandas as pd


ROOT = Path("${CLUSTER_WORK}/snp_repair_r308/fig4new")
FIG_DIR = ROOT / "figures"
SCRIPT_DIR = ROOT / "scripts"
METADATA_DIR = ROOT / "metadata"
SUMMARY_DIR = ROOT / "summary"

STRUCTURE_DIR = ROOT / "structure"
TREE_FILE = Path(
    "${ANALYSIS_DIR}/"
    "05_GWAS/00_analysis/02_genomeDB/06_iqtree/oilpalm_newick.txt"
)
TREE_ORDER_FILE = STRUCTURE_DIR / "integrated_rooted_tree_structure_sample_order.txt"
FAM_FILE = STRUCTURE_DIR / "all.fam"
PCA10_FILE = ROOT / "structure/PCA_10_oriented.eigenvec"

PC_VARIANCE = {
    "PC1": 6.93,   # eigenvalue / GRM trace (21.2803/307.208), as in the published Fig. 4b labels
    "PC2": 4.60,
    "PC3": 3.33,
}

K4_COLORS = {
    "K4_Pop1": "#70B4E7",
    "K4_Pop2": "#F3766A",
    "K4_Pop3": "#82C7B8",
    "K4_Pop4": "#BBA8D8",
}
K3_COLORS = {
    "K3_Pop1": "#70B4E7",
    "K3_Pop2": "#F2B76A",
    "K3_Pop3": "#82C7B8",
}
K6_COLORS = {
    "K6_Pop1": "#E7A8B8",
    "K6_Pop2": "#98B0D1",
    "K6_Pop3": "#A7D9DD",
    "K6_Pop4": "#FDE5B0",
    "K6_Pop5": "#C3B6E6",
    "K6_Pop6": "#F7C2B3",
}

STRUCTURE_COLORS = {
    3: [K3_COLORS[f"K3_Pop{i}"] for i in range(1, 4)],
    4: [K4_COLORS[f"K4_Pop{i}"] for i in range(1, 5)],
    6: [K6_COLORS[f"K6_Pop{i}"] for i in range(1, 7)],
}

EDGE_CMAP = LinearSegmentedColormap.from_list(
    "fst_reference_red", ["#F8D4C8", "#F08A78", "#D64A58", "#7F3348"]
)


def setup_style() -> None:
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "axes.linewidth": 0.65,
            "xtick.major.width": 0.55,
            "ytick.major.width": 0.55,
        }
    )


def read_fam_samples() -> list[str]:
    fam = pd.read_csv(FAM_FILE, sep=r"\s+", header=None)
    return fam.iloc[:, 1].astype(str).tolist()


def read_tree_order() -> list[str]:
    with TREE_ORDER_FILE.open() as handle:
        return [line.strip() for line in handle if line.strip()]


def read_q(k: int, samples: list[str]) -> pd.DataFrame:
    q_path = STRUCTURE_DIR / f"all.{k}.Q"
    if not q_path.exists():
        q_path = STRUCTURE_DIR / "result" / f"all.{k}.Q"
    q = pd.read_csv(q_path, sep=r"\s+", header=None)
    q.index = samples
    q.columns = [f"K{k}_Pop{i}" for i in range(1, k + 1)]
    return q


def read_k4_assignments() -> pd.DataFrame:
    df = pd.read_csv(METADATA_DIR / "K4_dominant_assignments.tsv", sep="\t")
    for col in ["PC1", "PC2", "PC3", "max_Q", "second_Q"]:
        df[col] = pd.to_numeric(df[col], errors="coerce")
    return df


def group_counts(assignments: pd.DataFrame) -> dict[str, int]:
    return assignments["dominant_group"].value_counts().to_dict()


def format_group_label(group: str, counts: dict[str, int]) -> str:
    return f"{group.replace('_', '-')}" + f" (n={counts.get(group, 0)})"


def lighten_color(color: str, amount: float = 0.35) -> str:
    rgb = np.array(mpl.colors.to_rgb(color))
    rgb = rgb + (1.0 - rgb) * amount
    return mpl.colors.to_hex(rgb)


def contiguous_runs(labels: list[str]) -> list[tuple[int, int, str]]:
    if not labels:
        return []
    runs = []
    start = 0
    current = labels[0]
    for i, label in enumerate(labels[1:], start=1):
        if label != current:
            runs.append((start, i - 1, current))
            start = i
            current = label
    runs.append((start, len(labels) - 1, current))
    return runs


def local_majority_labels(labels: list[str], half_window: int = 20) -> list[str]:
    smoothed = []
    for i, label in enumerate(labels):
        left = max(0, i - half_window)
        right = min(len(labels), i + half_window + 1)
        counts = Counter(labels[left:right])
        smoothed.append(max(counts, key=lambda group: (counts[group], group == label)))
    return smoothed


def clade_tip_names(clade) -> set[str]:
    return {tip.name for tip in clade.get_terminals()}


def draw_tree_structure() -> Path:
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    samples = read_fam_samples()
    order = [sample for sample in read_tree_order() if sample in samples]
    order.extend([sample for sample in samples if sample not in set(order)])
    q_by_k = {k: read_q(k, samples).loc[order] for k in (3, 4, 6)}
    assignments_df = read_k4_assignments()
    assignments_df = assignments_df[assignments_df["sample"].isin(samples)].copy()
    assignments = assignments_df.set_index("sample")
    k4_group = assignments["dominant_group"].to_dict()
    counts = group_counts(assignments.reset_index())

    tree = Phylo.read(TREE_FILE, "newick")
    valid_tips = set(order)
    for tip in list(tree.get_terminals()):
        if tip.name not in valid_tips:
            tree.prune(tip)

    x_pos = {sample: i for i, sample in enumerate(order)}
    node_x: dict[object, float] = {}
    node_y: dict[object, float] = {}
    node_groups: dict[object, set[str]] = {}

    def recurse(clade):
        if clade.is_terminal():
            node_x[clade] = x_pos.get(clade.name, 0)
            node_y[clade] = 0.0
            node_groups[clade] = {k4_group.get(clade.name, "Unknown")}
        else:
            for child in clade.clades:
                recurse(child)
            child_x = [node_x[ch] for ch in clade.clades if ch in node_x]
            node_x[clade] = float(np.mean(child_x)) if child_x else 0.0
            child_y = [node_y[ch] for ch in clade.clades if ch in node_y]
            node_y[clade] = (max(child_y) + 1.0) if child_y else 1.0
            groups = set()
            for child in clade.clades:
                groups.update(node_groups.get(child, {"Unknown"}))
            node_groups[clade] = groups

    recurse(tree.root)
    max_depth = max(node_y.values()) if node_y else 1.0

    fig = plt.figure(figsize=(20.2, 6.15))
    gs = fig.add_gridspec(
        nrows=5,
        ncols=1,
        height_ratios=[1.55, 2.25, 0.30, 0.10, 0.33],
        hspace=0.012,
    )
    ax_tree = fig.add_subplot(gs[0])
    ax_struct = fig.add_subplot(gs[1], sharex=ax_tree)
    ax_labels = fig.add_subplot(gs[2], sharex=ax_tree)
    ax_strip = fig.add_subplot(gs[3], sharex=ax_tree)
    ax_legend = fig.add_subplot(gs[4])

    ordered_groups = [k4_group.get(sample, "Unknown") for sample in order]
    background_groups = local_majority_labels(ordered_groups, half_window=20)
    block_boundaries = []
    for start, end, group in contiguous_runs(background_groups):
        if group not in K4_COLORS:
            continue
        ax_tree.axvspan(
            start - 0.5,
            end + 0.5,
            color=K4_COLORS[group],
            alpha=0.18,
            lw=0,
            zorder=-5,
        )
        block_boundaries.append(end + 0.5)

    def roundrect_edge(x0: float, y0: float, x1: float, y1: float, lw: float) -> None:
        dx = x1 - x0
        dy = y0 - y1
        if abs(dx) < 1e-6 or dy <= 0:
            verts = [(x0, y0), (x1, y0), (x1, y1)]
            codes = [MplPath.MOVETO, MplPath.LINETO, MplPath.LINETO]
        else:
            sign = 1 if dx > 0 else -1
            rx = min(max(abs(dx) * 0.24, 0.42), 2.40)
            ry = min(max(dy * 0.42, 0.20), 1.35)
            verts = [
                (x0, y0),
                (x1 - sign * rx, y0),
                (x1, y0),
                (x1, y0 - ry),
                (x1, y1),
            ]
            codes = [
                MplPath.MOVETO,
                MplPath.LINETO,
                MplPath.CURVE3,
                MplPath.CURVE3,
                MplPath.LINETO,
            ]
        path = MplPath(verts, codes)
        ax_tree.add_patch(
            PathPatch(
                path,
                facecolor="none",
                edgecolor="#6E7479",
                lw=lw,
                alpha=0.68,
                capstyle="round",
                joinstyle="round",
                zorder=2,
            )
        )

    for parent in tree.find_clades(order="preorder"):
        if parent.is_terminal() or not parent.clades:
            continue
        for child in parent.clades:
            if child not in node_x:
                continue
            roundrect_edge(
                node_x[parent],
                node_y[parent],
                node_x[child],
                node_y[child],
                0.28 if child.is_terminal() else 0.36,
            )
    terminal_x = [x_pos[s] for s in order]
    terminal_colors = [K4_COLORS.get(k4_group.get(s, ""), "#BDBDBD") for s in order]
    ax_tree.scatter(
        terminal_x,
        [0.0] * len(terminal_x),
        s=3.8,
        c=terminal_colors,
        linewidths=0,
        zorder=3,
        alpha=0.88,
    )
    ax_tree.set_xlim(-1, len(order))
    ax_tree.set_ylim(-0.03 * max_depth, max_depth * 1.03)
    ax_tree.axis("off")
    for boundary in block_boundaries[:-1]:
        ax_tree.axvline(boundary, color="white", lw=0.7, alpha=0.85, zorder=-1)

    # STRUCTURE bars for K=3, K=4 and K=6.
    ax_struct.set_ylim(0, 3)
    row_height = 1.0
    row_ys = {3: 2.0, 4: 1.0, 6: 0.0}
    for k in (3, 4, 6):
        q = q_by_k[k]
        y0 = row_ys[k]
        for i, (_sample, vals) in enumerate(q.iterrows()):
            cum = 0.0
            for j, value in enumerate(vals.values):
                h = float(value) * row_height
                if h <= 0:
                    continue
                ax_struct.add_patch(
                    Rectangle(
                        (i - 0.5, y0 + cum * row_height),
                        1.0,
                        h,
                        facecolor=STRUCTURE_COLORS[k][j],
                        edgecolor="none",
                        linewidth=0,
                    )
                )
                cum += float(value)
        ax_struct.axhline(y0, color="white", lw=0.45, alpha=0.85)
        ax_struct.text(
            -3.2,
            y0 + row_height / 2,
            f"K={k}",
            ha="right",
            va="center",
            fontsize=8.5,
            fontweight="bold",
        )

    ax_struct.axhline(3, color="white", lw=0.45, alpha=0.85)
    for boundary in block_boundaries[:-1]:
        ax_struct.axvline(boundary, color="white", lw=0.55, alpha=0.72)
    ax_struct.set_xlim(-1, len(order))
    ax_struct.axis("off")

    ax_labels.set_ylim(0, 1)
    for i, sample in enumerate(order):
        ax_labels.text(
            i,
            0.98,
            sample,
            rotation=90,
            ha="center",
            va="top",
            fontsize=2.15,
            color="#343434",
            clip_on=False,
        )
    ax_labels.axis("off")

    ax_strip.set_ylim(0, 1)
    for i, sample in enumerate(order):
        ax_strip.add_patch(
            Rectangle(
                (i - 0.5, 0),
                1,
                1,
                facecolor=K4_COLORS.get(k4_group.get(sample, ""), "#D0D0D0"),
                edgecolor="none",
            )
        )
    for boundary in block_boundaries[:-1]:
        ax_strip.axvline(boundary, color="white", lw=0.65, alpha=0.85)
    ax_strip.axis("off")

    ax_legend.axis("off")
    handles = [
        Patch(facecolor=K4_COLORS[group], edgecolor="none", label=format_group_label(group, counts))
        for group in ["K4_Pop1", "K4_Pop2", "K4_Pop3", "K4_Pop4"]
    ]
    ax_legend.legend(
        handles=handles,
        loc="center",
        ncol=4,
        frameon=False,
        fontsize=7.7,
        handlelength=1.0,
        columnspacing=1.4,
    )
    out = FIG_DIR / "requested_01_tree_structure_K3_K4_K6_K4treecolor.pdf"
    fig.savefig(out, bbox_inches="tight", pad_inches=0.03)
    fig.savefig(out.with_suffix(".png"), dpi=320, bbox_inches="tight", pad_inches=0.03)
    plt.close(fig)
    return out


def axis_label(pc: str) -> str:
    return f"{pc} ({PC_VARIANCE[pc]:.2f}%)"


def add_group_ellipse(ax, sub: pd.DataFrame, color: str) -> None:
    pts = sub[["PC1", "PC2"]].dropna().values
    if len(pts) < 4:
        return
    center = np.median(pts, axis=0)
    if len(pts) > 15:
        dist = np.sqrt(((pts - center) ** 2).sum(axis=1))
        pts = pts[dist <= np.quantile(dist, 0.86)]
    if len(pts) < 4:
        return
    cov = np.cov(pts, rowvar=False)
    if not np.all(np.isfinite(cov)):
        return
    vals, vecs = np.linalg.eigh(cov)
    vals = np.maximum(vals, 1e-10)
    order = vals.argsort()[::-1]
    vals, vecs = vals[order], vecs[:, order]
    angle = math.degrees(math.atan2(vecs[1, 0], vecs[0, 0]))
    width, height = 2.4 * np.sqrt(vals)
    width = max(width, 0.018)
    height = max(height, 0.010)
    fill = Ellipse(
        xy=center,
        width=width,
        height=height,
        angle=angle,
        facecolor=color,
        edgecolor="none",
        alpha=0.18,
        zorder=0,
    )
    outline = Ellipse(
        xy=center,
        width=width,
        height=height,
        angle=angle,
        facecolor="none",
        edgecolor=color,
        linewidth=1.75,
        alpha=0.68,
        zorder=6,
    )
    ax.add_patch(fill)
    ax.add_patch(outline)


def scatter_k4_pc12(
    ax,
    df: pd.DataFrame,
    counts: dict[str, int],
    add_ellipses: bool = True,
    size_scale: float = 1.0,
    legend_labels: bool = False,
) -> None:
    draw_order = ["K4_Pop4", "K4_Pop3", "K4_Pop1", "K4_Pop2"]
    for group in draw_order:
        sub_all = df[df["dominant_group"] == group]
        if add_ellipses:
            add_group_ellipse(ax, sub_all, K4_COLORS[group])
        sub_admixed = sub_all[sub_all["max_Q"] < 0.7]
        sub_core = sub_all[sub_all["max_Q"] >= 0.7]
        if not sub_admixed.empty:
            ax.scatter(
                sub_admixed["PC1"],
                sub_admixed["PC2"],
                s=(46 if group != "K4_Pop2" else 54) * size_scale,
                facecolor=K4_COLORS[group],
                edgecolor="white",
                linewidth=0.42,
                alpha=0.38,
                zorder=2,
            )
        if not sub_core.empty:
            ax.scatter(
                sub_core["PC1"],
                sub_core["PC2"],
                s=(72 if group != "K4_Pop2" else 84) * size_scale,
                facecolor=K4_COLORS[group],
                edgecolor="white",
                linewidth=0.58,
                alpha=0.9,
                zorder=4 if group == "K4_Pop2" else 3,
                label=format_group_label(group, counts) if legend_labels else None,
            )


def style_pca_axis(ax) -> None:
    ax.axhline(0, color="#D7DCE0", lw=0.75, zorder=0)
    ax.axvline(0, color="#D7DCE0", lw=0.75, zorder=0)
    ax.grid(color="#EEF1F3", lw=0.6)
    for spine in ax.spines.values():
        spine.set_color("#3A3A3A")
        spine.set_linewidth(0.75)
    ax.tick_params(labelsize=8.3)


def draw_k4_pca() -> Path:
    df = read_k4_assignments().dropna(subset=["PC1", "PC2", "dominant_group"]).copy()
    counts = group_counts(df)
    handles = [
        Line2D(
            [0],
            [0],
            marker="o",
            linestyle="",
            markerfacecolor=K4_COLORS[group],
            markeredgecolor="white",
            markersize=8.7,
            label=format_group_label(group, counts),
        )
        for group in ["K4_Pop1", "K4_Pop2", "K4_Pop3", "K4_Pop4"]
    ]
    handles.append(
        Line2D(
            [0],
            [0],
            marker="o",
            linestyle="",
            markerfacecolor="#BBBBBB",
            markeredgecolor="white",
            markersize=7.4,
            alpha=0.45,
            label="max Q < 0.7",
        )
    )

    # Match the horizontal proportions of panel b in the composite reference
    # (approximately 1.52:1) and preserve that ratio in the exported page.
    fig, ax = plt.subplots(figsize=(7.0, 4.6))
    scatter_k4_pc12(ax, df, counts, add_ellipses=True, size_scale=1.0, legend_labels=True)
    ax.set_xlabel(axis_label("PC1"), fontsize=10)
    ax.set_ylabel(axis_label("PC2"), fontsize=10)
    style_pca_axis(ax)
    ax.legend(handles=handles, frameon=False, fontsize=8.1, loc="upper right", handletextpad=0.35)
    fig.subplots_adjust(left=0.105, right=0.985, bottom=0.145, top=0.975)
    out = FIG_DIR / "requested_02_K4_population_PCA_PC1_PC2.pdf"
    fig.savefig(out)
    fig.savefig(out.with_suffix(".png"), dpi=320)
    plt.close(fig)

    zoom_xlim = (-0.036, 0.056)
    zoom_ylim = (-0.024, 0.040)
    fig_zoom, ax_zoom = plt.subplots(figsize=(5.6, 5.3))
    scatter_k4_pc12(ax_zoom, df, counts, add_ellipses=False, size_scale=1.12, legend_labels=False)
    ax_zoom.set_xlim(*zoom_xlim)
    ax_zoom.set_ylim(*zoom_ylim)
    ax_zoom.set_xlabel(axis_label("PC1"), fontsize=10)
    ax_zoom.set_ylabel(axis_label("PC2"), fontsize=10)
    ax_zoom.set_title("Dense K=4 PCA cluster", loc="left", fontsize=12, fontweight="bold")
    style_pca_axis(ax_zoom)
    fig_zoom.tight_layout()
    out_zoom = FIG_DIR / "requested_02_K4_population_PCA_PC1_PC2_zoom.pdf"
    fig_zoom.savefig(out_zoom, bbox_inches="tight")
    fig_zoom.savefig(out_zoom.with_suffix(".png"), dpi=320, bbox_inches="tight")
    plt.close(fig_zoom)
    return out


def read_pc10_for_embedding() -> pd.DataFrame:
    pc = pd.read_csv(PCA10_FILE, sep=r"\s+", header=None)
    pc.columns = ["fid", "sample"] + [f"PC{i}" for i in range(1, pc.shape[1] - 1)]
    return pc


def embedding_pc_columns(df: pd.DataFrame, n_pcs: int = 10) -> list[str]:
    cols = []
    for i in range(1, n_pcs + 1):
        embed_col = f"PC{i}_embed"
        raw_col = f"PC{i}"
        if embed_col in df.columns:
            cols.append(embed_col)
        elif raw_col in df.columns:
            cols.append(raw_col)
    return cols


def draw_k4_tsne_from_pcs() -> Path:
    from sklearn.manifold import TSNE
    from sklearn.preprocessing import StandardScaler

    assignments = read_k4_assignments()
    pc = read_pc10_for_embedding()
    df = assignments.merge(pc, on="sample", suffixes=("", "_embed")).dropna(subset=["dominant_group"]).copy()
    pc_cols = embedding_pc_columns(df, n_pcs=10)
    x = StandardScaler().fit_transform(df[pc_cols].astype(float).values)
    emb = TSNE(
        n_components=2,
        perplexity=30,
        init="pca",
        learning_rate="auto",
        max_iter=1500,
        random_state=42,
        metric="euclidean",
    ).fit_transform(x)
    df["tSNE1"] = emb[:, 0]
    df["tSNE2"] = emb[:, 1]
    counts = group_counts(df)

    fig, ax = plt.subplots(figsize=(7.8, 5.6))
    for group in ["K4_Pop4", "K4_Pop3", "K4_Pop1", "K4_Pop2"]:
        sub_all = df[df["dominant_group"] == group]
        sub_admixed = sub_all[sub_all["max_Q"] < 0.7]
        sub_core = sub_all[sub_all["max_Q"] >= 0.7]
        if not sub_admixed.empty:
            ax.scatter(
                sub_admixed["tSNE1"],
                sub_admixed["tSNE2"],
                s=42,
                facecolor=K4_COLORS[group],
                edgecolor="white",
                linewidth=0.4,
                alpha=0.34,
                zorder=2,
            )
        if not sub_core.empty:
            ax.scatter(
                sub_core["tSNE1"],
                sub_core["tSNE2"],
                s=68 if group != "K4_Pop2" else 80,
                facecolor=K4_COLORS[group],
                edgecolor="white",
                linewidth=0.55,
                alpha=0.9,
                zorder=4 if group == "K4_Pop2" else 3,
            )
    ax.set_xlabel("t-SNE 1", fontsize=10)
    ax.set_ylabel("t-SNE 2", fontsize=10)
    ax.set_title("K=4 groups on t-SNE of the first 10 PCs", loc="left", fontsize=13, fontweight="bold")
    ax.grid(color="#EEF1F3", lw=0.6)
    for spine in ax.spines.values():
        spine.set_color("#3A3A3A")
        spine.set_linewidth(0.75)
    handles = [
        Line2D(
            [0],
            [0],
            marker="o",
            linestyle="",
            markerfacecolor=K4_COLORS[group],
            markeredgecolor="white",
            markersize=8.2,
            label=format_group_label(group, counts),
        )
        for group in ["K4_Pop1", "K4_Pop2", "K4_Pop3", "K4_Pop4"]
    ]
    ax.legend(
        handles=handles,
        frameon=False,
        fontsize=8,
        loc="upper left",
        bbox_to_anchor=(1.01, 1.0),
        borderaxespad=0,
        handletextpad=0.35,
    )
    fig.subplots_adjust(right=0.78, bottom=0.13)
    fig.text(
        0.10,
        0.025,
        "t-SNE was computed from standardized PC1-PC10 scores and is used as an exploratory projection.",
        ha="left",
        va="bottom",
        fontsize=7.2,
        color="#666666",
    )
    out = FIG_DIR / "requested_02b_K4_population_tSNE_first10PCs.pdf"
    fig.savefig(out, bbox_inches="tight")
    fig.savefig(out.with_suffix(".png"), dpi=320, bbox_inches="tight")
    plt.close(fig)
    return out


def draw_k4_umap_from_pcs() -> Path:
    from sklearn.preprocessing import StandardScaler
    from umap import UMAP

    assignments = read_k4_assignments()
    pc = read_pc10_for_embedding()
    df = assignments.merge(pc, on="sample", suffixes=("", "_embed")).dropna(subset=["dominant_group"]).copy()
    pc_cols = embedding_pc_columns(df, n_pcs=10)
    x = StandardScaler().fit_transform(df[pc_cols].astype(float).values)
    emb = UMAP(
        n_components=2,
        n_neighbors=30,
        min_dist=0.22,
        spread=1.0,
        metric="euclidean",
        random_state=42,
        transform_seed=42,
    ).fit_transform(x)
    df["UMAP1"] = emb[:, 0]
    df["UMAP2"] = emb[:, 1]
    counts = group_counts(df)

    fig, ax = plt.subplots(figsize=(7.8, 5.6))
    for group in ["K4_Pop4", "K4_Pop3", "K4_Pop1", "K4_Pop2"]:
        sub_all = df[df["dominant_group"] == group]
        sub_admixed = sub_all[sub_all["max_Q"] < 0.7]
        sub_core = sub_all[sub_all["max_Q"] >= 0.7]
        if not sub_admixed.empty:
            ax.scatter(
                sub_admixed["UMAP1"],
                sub_admixed["UMAP2"],
                s=42,
                facecolor=K4_COLORS[group],
                edgecolor="white",
                linewidth=0.4,
                alpha=0.34,
                zorder=2,
            )
        if not sub_core.empty:
            ax.scatter(
                sub_core["UMAP1"],
                sub_core["UMAP2"],
                s=68 if group != "K4_Pop2" else 80,
                facecolor=K4_COLORS[group],
                edgecolor="white",
                linewidth=0.55,
                alpha=0.9,
                zorder=4 if group == "K4_Pop2" else 3,
            )
    ax.set_xlabel("UMAP 1", fontsize=10)
    ax.set_ylabel("UMAP 2", fontsize=10)
    ax.set_title("K=4 groups on UMAP of the first 10 PCs", loc="left", fontsize=13, fontweight="bold")
    ax.grid(color="#EEF1F3", lw=0.6)
    for spine in ax.spines.values():
        spine.set_color("#3A3A3A")
        spine.set_linewidth(0.75)
    handles = [
        Line2D(
            [0],
            [0],
            marker="o",
            linestyle="",
            markerfacecolor=K4_COLORS[group],
            markeredgecolor="white",
            markersize=8.2,
            label=format_group_label(group, counts),
        )
        for group in ["K4_Pop1", "K4_Pop2", "K4_Pop3", "K4_Pop4"]
    ]
    ax.legend(
        handles=handles,
        frameon=False,
        fontsize=8,
        loc="upper left",
        bbox_to_anchor=(1.01, 1.0),
        borderaxespad=0,
        handletextpad=0.35,
    )
    fig.subplots_adjust(right=0.78, bottom=0.13)
    fig.text(
        0.10,
        0.025,
        "UMAP was computed from standardized PC1-PC10 scores and is used as an exploratory projection.",
        ha="left",
        va="bottom",
        fontsize=7.2,
        color="#666666",
    )
    out = FIG_DIR / "requested_02c_K4_population_UMAP_first10PCs.pdf"
    fig.savefig(out, bbox_inches="tight")
    fig.savefig(out.with_suffix(".png"), dpi=320, bbox_inches="tight")
    plt.close(fig)
    return out


def robust_pc_points(sub: pd.DataFrame, cols: list[str], keep_quantile: float = 0.86) -> np.ndarray:
    pts = sub[cols].dropna().values
    if len(pts) > 15:
        center = np.median(pts, axis=0)
        dist = np.sqrt(((pts - center) ** 2).sum(axis=1))
        pts = pts[dist <= np.quantile(dist, keep_quantile)]
    return pts


def ellipsoid_segments(
    center: np.ndarray,
    cov: np.ndarray,
    min_radii: np.ndarray,
    scale: float = 1.95,
    steps: int = 48,
) -> list[np.ndarray]:
    vals, vecs = np.linalg.eigh(cov)
    vals = np.maximum(vals, 1e-10)
    order = vals.argsort()[::-1]
    vals, vecs = vals[order], vecs[:, order]
    radii = np.maximum(scale * np.sqrt(vals), min_radii)
    segments = []
    theta = np.linspace(0, 2 * np.pi, steps)
    for plane in [(0, 1), (0, 2), (1, 2)]:
        local = np.zeros((steps, 3))
        local[:, plane[0]] = radii[plane[0]] * np.cos(theta)
        local[:, plane[1]] = radii[plane[1]] * np.sin(theta)
        segments.append(local @ vecs.T + center)
    return segments


def draw_k4_pca_3d() -> Path:
    df = read_k4_assignments().dropna(subset=["PC1", "PC2", "PC3", "dominant_group"]).copy()
    counts = group_counts(df)
    fig = plt.figure(figsize=(7.3, 6.15))
    ax = fig.add_subplot(111, projection="3d")
    draw_order = ["K4_Pop4", "K4_Pop3", "K4_Pop1", "K4_Pop2"]

    for group in draw_order:
        sub_all = df[df["dominant_group"] == group]
        pts = robust_pc_points(sub_all, ["PC1", "PC2", "PC3"])
        if len(pts) >= 4:
            center = np.median(pts, axis=0)
            cov = np.cov(pts, rowvar=False)
            if np.all(np.isfinite(cov)):
                segments = ellipsoid_segments(center, cov, np.array([0.015, 0.010, 0.010]))
                collection = Line3DCollection(
                    segments,
                    colors=K4_COLORS[group],
                    linewidths=1.15,
                    alpha=0.52,
                    zorder=1,
                )
                ax.add_collection3d(collection)

        sub_admixed = sub_all[sub_all["max_Q"] < 0.7]
        sub_core = sub_all[sub_all["max_Q"] >= 0.7]
        if not sub_admixed.empty:
            ax.scatter(
                sub_admixed["PC1"],
                sub_admixed["PC2"],
                sub_admixed["PC3"],
                s=18 if group != "K4_Pop2" else 24,
                c=K4_COLORS[group],
                edgecolors="white",
                linewidths=0.25,
                alpha=0.32,
                depthshade=False,
            )
        if not sub_core.empty:
            ax.scatter(
                sub_core["PC1"],
                sub_core["PC2"],
                sub_core["PC3"],
                s=30 if group != "K4_Pop2" else 40,
                c=K4_COLORS[group],
                edgecolors="white",
                linewidths=0.35,
                alpha=0.88,
                depthshade=False,
            )

    ax.set_xlabel(axis_label("PC1"), labelpad=7, fontsize=9)
    ax.set_ylabel(axis_label("PC2"), labelpad=7, fontsize=9)
    ax.set_zlabel(axis_label("PC3"), labelpad=7, fontsize=9)
    ax.set_title("K=4 population groups in 3D PCA space", loc="left", fontsize=13.5, fontweight="bold")
    ax.view_init(elev=21, azim=-58)
    axis_limits = {
        "PC1": (float(df["PC1"].min()) - 0.015, float(df["PC1"].max()) + 0.015),
        "PC2": (float(df["PC2"].min()) - 0.030, float(df["PC2"].max()) + 0.030),
        "PC3": (float(df["PC3"].min()) - 0.040, float(df["PC3"].max()) + 0.040),
    }
    ax.set_xlim(*axis_limits["PC1"])
    ax.set_ylim(*axis_limits["PC2"])
    ax.set_zlim(*axis_limits["PC3"])
    ax.set_box_aspect((1.05, 1.00, 0.72))
    ax.tick_params(labelsize=7.4, pad=0)
    for axis in (ax.xaxis, ax.yaxis, ax.zaxis):
        axis.pane.set_facecolor((1, 1, 1, 0))
        axis.pane.set_edgecolor("#D8DDE1")
        axis._axinfo["grid"]["color"] = "#E8ECEF"
        axis._axinfo["grid"]["linewidth"] = 0.55
    handles = [
        Line2D(
            [0],
            [0],
            marker="o",
            linestyle="",
            markerfacecolor=K4_COLORS[group],
            markeredgecolor="white",
            markersize=7,
            label=format_group_label(group, counts),
        )
        for group in ["K4_Pop1", "K4_Pop2", "K4_Pop3", "K4_Pop4"]
    ]
    handles.append(
        Line2D(
            [0],
            [0],
            marker="o",
            linestyle="",
            markerfacecolor="#BBBBBB",
            markeredgecolor="white",
            markersize=6,
            alpha=0.45,
            label="max Q < 0.7",
        )
    )
    ax.legend(handles=handles, frameon=False, fontsize=8.1, loc="upper left", bbox_to_anchor=(0.02, 0.96))
    fig.tight_layout()
    out = FIG_DIR / "requested_04_K4_population_PCA_3D.pdf"
    fig.savefig(out, bbox_inches="tight")
    fig.savefig(out.with_suffix(".png"), dpi=320, bbox_inches="tight")
    plt.close(fig)
    return out


def read_pi_table(path: Path) -> pd.Series:
    df = pd.read_csv(path, sep="\t")
    return pd.Series(df["PI_weighted"].astype(float).values, index=df["group"].astype(str))


def read_fst_matrix(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", index_col=0)
    df.index = df.index.astype(str)
    df.columns = df.columns.astype(str)
    return df.astype(float)


def draw_k4_network() -> Path:
    # Balanced square layout for use as panel c in the composite.
    # The group legend is already supplied by the neighbouring PCA panel, so this
    # panel only retains the Fst scale and places the pi values inside the nodes.
    order = ["K4_Pop1", "K4_Pop4", "K4_Pop3", "K4_Pop2"]
    pi = read_pi_table(SUMMARY_DIR / "K4_PI_weighted.tsv").reindex(order)
    fst = read_fst_matrix(SUMMARY_DIR / "K4_FST_weighted_matrix.tsv").loc[order, order]

    vals = [float(fst.iloc[i, j]) for i in range(len(order)) for j in range(i + 1, len(order))]
    norm = Normalize(vmin=min(vals), vmax=max(vals))
    positions = {
        "K4_Pop1": np.array([-0.82, 0.72]),
        "K4_Pop4": np.array([0.62, 0.72]),
        "K4_Pop3": np.array([0.62, -0.72]),
        "K4_Pop2": np.array([-0.82, -0.72]),
    }

    fig, ax = plt.subplots(figsize=(6.0, 6.0))
    fig.subplots_adjust(left=0.02, right=0.98, bottom=0.02, top=0.98)
    ax.set_aspect("equal", adjustable="box")
    ax.axis("off")

    for i, left in enumerate(order):
        for j, right in enumerate(order):
            if j <= i:
                continue
            value = float(fst.loc[left, right])
            t = norm(value)
            p1, p2 = positions[left], positions[right]
            ax.plot(
                [p1[0], p2[0]],
                [p1[1], p2[1]],
                color=EDGE_CMAP(0.16 + 0.84 * t),
                lw=2.4 + 5.8 * t,
                alpha=0.94,
                solid_capstyle="round",
                zorder=1,
            )

    pi_min, pi_max = float(pi.min()), float(pi.max())
    for group in order:
        value = float(pi.loc[group])
        t = 0.5 if pi_max <= pi_min else (value - pi_min) / (pi_max - pi_min)
        # Keep the pi encoding while avoiding an undersized node whose two-line
        # label would spill outside the circle in the square panel.
        size = 1900 + 1700 * t
        x, y = positions[group]
        ax.scatter(
            x,
            y,
            s=size,
            facecolor=K4_COLORS[group],
            edgecolor="white",
            linewidth=1.5,
            zorder=4,
        )
        ax.text(
            x,
            y,
            group.replace("K4_", "") + f"\n{value * 1000:.2f} ×10" + r"$^{-3}$",
            ha="center",
            va="center",
            fontsize=8.5,
            color="#1B1B1B",
            zorder=5,
        )

    ax.set_xlim(-1.22, 1.22)
    ax.set_ylim(-1.22, 1.22)
    sm = mpl.cm.ScalarMappable(norm=norm, cmap=EDGE_CMAP)
    sm.set_array([])
    cax = ax.inset_axes([0.91, 0.32, 0.026, 0.36])
    cbar = fig.colorbar(sm, cax=cax)
    cbar.ax.set_title("Fst", fontsize=10, loc="left", pad=5)
    cbar.ax.tick_params(labelsize=7.5, width=0.5, length=2.5)
    out = FIG_DIR / "requested_03_K4_pi_fst_network.pdf"
    # Keep an exact square page; tight-bbox cropping would alter that proportion.
    fig.savefig(out)
    fig.savefig(out.with_suffix(".png"), dpi=340)
    plt.close(fig)
    return out


def main() -> None:
    import sys
    setup_style()
    if "tree" in sys.argv: print(draw_tree_structure())
    print(draw_k4_pca())
    print(draw_k4_network())


if __name__ == "__main__":
    main()
