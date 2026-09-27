#!/usr/bin/env python3
"""Rebuild corrected Figure 3 panels without changing upstream files."""

from __future__ import annotations

import hashlib
import math
from itertools import combinations
from pathlib import Path

import matplotlib as mpl

mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
from matplotlib.lines import Line2D
from matplotlib.patches import PathPatch, Rectangle
from matplotlib.path import Path as MplPath
from matplotlib.ticker import FuncFormatter
import numpy as np
import pandas as pd
from scipy import stats


OUT = Path(__file__).resolve().parent.parent
ANALYSIS = Path(
    "${ANALYSIS_DIR}/22_answer_reviews/00_ms/03_V3"
)
REF = ANALYSIS / "03_figure3/04_FL_TN_ASE_reference_panels"
CURRENT = ANALYSIS / "02_figure/_runs/RUN-ASE-FIGURES-CURRENT-V3-001/output"
HET = (
    ANALYSIS
    / "03_figure3/01_ASE/00_shared/runs/RUN-TN-EXPR-HETEROSIS-V2-001/output"
)
ASE_OUTPUT = (
    ANALYSIS
    / "03_figure3/01_ASE/00_shared/runs/"
    "RUN-ASE-HETEROSIS-DOWNSTREAM-V2-001/output"
)

CLASS_SOURCE = REF / "source_B_ASE_class_proportions.tsv"
MIRROR_SOURCE = REF / "source_D_stage_mirror.tsv"
FLOW_SOURCE = REF / "source_E_FL_TN_alluvial.tsv"
TEST_REFERENCE = REF / "source_G_SNP_density_tests.tsv"
MODE_LEGACY_SUMMARY = HET / "stage_mode_summary.tsv"
MODE_DETAIL_SOURCE = HET / "modes_classified.tsv"
REGULATORY_SOURCE = ASE_OUTPUT / "TN_cis_trans_classified.tsv"
TRAIT_SOURCE = (
    ANALYSIS
    / "03_figure3/01_ASE/00_shared/runs/"
    "RUN-ASE-HETEROSIS-DOWNSTREAM-V2-001/output/trait_haplotype_ASE.tsv.gz"
)
LEGACY_INHERITANCE_SOURCE = (
    REF / "source_04D_TN_inheritance_regulatory_association.tsv"
)
TRAIT_LEGACY_SUMMARY = REF / "source_06A_FL_TN_trait_complement.tsv"
SNP_PATHS = {
    "FL": CURRENT / "FL_diagnostic_SNP_density_current.tsv.gz",
    "TN": CURRENT / "TN_diagnostic_SNP_density_current.tsv.gz",
}

STAGES = [
    "0d", "15d", "35d", "50d", "65d", "80d", "95d", "110d", "125d",
    "140d", "155d", "170d", "185d", "12h", "24h", "36h", "48h",
    "60h", "72h",
]
CLASS_ORDER = ["HapDom", "Sub", "NoDiff", "NoASE"]
CLASS_LABEL = {
    "HapDom": "HapDom",
    "Sub": "Sub",
    "NoDiff": "NoDiff",
    "NoASE": "NoASE",
}
CLASS_LONG_LABEL = {
    "HapDom": "HapDom (stable haplotype bias)",
    "Sub": "Sub (bias switching)",
    "NoDiff": "NoDiff (stage-limited or weakly switching ASE)",
    "NoASE": "NoASE (no robust ASE)",
}
CLASS_COLOR = {
    "HapDom": "#FFF59D",
    "Sub": "#8DD3C7",
    "NoDiff": "#F8766D",
    "NoASE": "#4EA5DF",
}
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
INHERITANCE_ORDER = ["PDO", "DO", "ODO"]
PHASE_ORDER = ["Days 0–65", "Days 80–140", "Days 155–185", "Hours 12–72"]
MODULE_ORDER = [
    "Oil biosynthesis & storage",
    "TAG assembly & oil body",
    "De-novo / saturated FA",
    "Unsaturated FA",
    "Lipid oxidation / antioxidant",
    "Shell / cell wall / lignin",
]
MODULE_SHORT = {
    "Oil biosynthesis & storage": "OBS",
    "TAG assembly & oil body": "TOF",
    "De-novo / saturated FA": "DSF",
    "Unsaturated FA": "UFA",
    "Lipid oxidation / antioxidant": "LOD",
    "Shell / cell wall / lignin": "SCL",
}

TN_A = "#68ADD0"
TN_B = "#DF8796"
FL_A = "#99BFDA"
FL_B = "#ECD48D"
DARK = "#27323A"
GRID = "#D9E0E4"
DPM_HEATMAP_CMAP = LinearSegmentedColormap.from_list(
    "dpm_blue_gold_red",
    ["#3767A6", "#8FAED0", "#F7F5EE", "#D8B64C", "#E56565"],
)


def set_style() -> None:
    mpl.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Liberation Sans", "DejaVu Sans"],
            "font.size": 8.2,
            "axes.titlesize": 9.5,
            "axes.labelsize": 8.7,
            "xtick.labelsize": 7.2,
            "ytick.labelsize": 7.2,
            "legend.fontsize": 7.1,
            "axes.linewidth": 0.75,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "savefig.facecolor": "white",
            "figure.facecolor": "white",
            "axes.facecolor": "white",
        }
    )


def clean(ax: plt.Axes, grid: str | None = None) -> None:
    ax.spines[["top", "right"]].set_visible(False)
    if grid:
        ax.grid(axis=grid, color=GRID, lw=0.5, zorder=0)
        ax.set_axisbelow(True)


def panel_label(ax: plt.Axes, label: str, x: float = -0.09, y: float = 1.04) -> None:
    ax.text(
        x,
        y,
        label,
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        fontsize=12.5,
        fontweight="bold",
        color="black",
    )


def save(fig: plt.Figure, stem: str, *, tight: bool = True) -> None:
    crop = {"bbox_inches": "tight", "pad_inches": 0.045} if tight else {}
    fig.savefig(OUT / f"{stem}.pdf", facecolor="white", **crop)
    fig.savefig(OUT / f"{stem}.png", dpi=600, facecolor="white", **crop)
    plt.close(fig)


def write_tsv(df: pd.DataFrame, name: str) -> None:
    df.to_csv(OUT / name, sep="\t", index=False)


def redraw_3d_formal() -> tuple[pd.DataFrame, list[dict[str, object]]]:
    data = pd.read_csv(CLASS_SOURCE, sep="\t")
    data["overall_class"] = pd.Categorical(
        data["overall_class"], categories=CLASS_ORDER, ordered=True
    )
    data = data.sort_values(["analysis", "overall_class"]).reset_index(drop=True)
    totals = data.groupby("analysis", observed=True)["genes"].sum()
    checks: list[dict[str, object]] = []
    for analysis in ["FL", "TN"]:
        part = data[data.analysis.eq(analysis)]
        checks.append(
            {
                "panel": "3d_formal",
                "check": f"{analysis}_four_classes_sum_to_100",
                "observed": float(part.percentage.sum()),
                "expected": 100.0,
                "pass": bool(np.isclose(part.percentage.sum(), 100.0, atol=1e-9)),
            }
        )

    # FL and TN are separate full universes. Keep the paired-column visual
    # language of the final panel without inventing unsupported ribbons.
    fig, ax = plt.subplots(figsize=(2.15, 4.55), facecolor="white")
    x_positions = {"FL": 0.30, "TN": 0.70}
    bottom = {"FL": 0.0, "TN": 0.0}
    for category in CLASS_ORDER:
        for analysis in ["FL", "TN"]:
            row = data[
                data.analysis.eq(analysis) & data.overall_class.eq(category)
            ].iloc[0]
            value = float(row.percentage)
            ax.bar(
                x_positions[analysis],
                value,
                bottom=bottom[analysis],
                width=0.22,
                color=CLASS_COLOR[category],
                edgecolor=DARK,
                lw=0.55,
                zorder=3,
            )
            is_small_segment = value < 3.0
            label_y = (
                bottom[analysis] + value + 0.7
                if is_small_segment
                else bottom[analysis] + value / 2
            )
            ax.text(
                x_positions[analysis],
                label_y,
                f"{value:.1f}%",
                ha="center",
                va="bottom" if is_small_segment else "center",
                fontsize=4.0 if is_small_segment else 6.2,
                color=DARK,
                clip_on=False,
            )
            bottom[analysis] += value

    ax.set_ylim(-5, 104)
    ax.set_xlim(0.08, 0.92)
    ax.set_xticks(
        [0.30, 0.70],
        [f"FL\n(n = {int(totals['FL']):,})", f"TN\n(n = {int(totals['TN']):,})"],
    )
    ax.set_yticks([])
    ax.tick_params(axis="x", length=0, pad=2)
    for spine in ax.spines.values():
        spine.set_visible(False)
    # Matplotlib fills a two-row legend column-wise. Reorder the handles so
    # the visual reading order remains HapDom, Sub, NoDiff, NoASE.
    legend_order = ["HapDom", "NoDiff", "Sub", "NoASE"]
    handles = [
        Rectangle(
            (0, 0), 1, 1,
            facecolor=CLASS_COLOR[c], edgecolor=DARK, lw=0.4,
            label=CLASS_LABEL[c],
        )
        for c in legend_order
    ]
    fig.legend(
        handles=handles,
        frameon=False,
        ncol=2,
        loc="lower center",
        bbox_to_anchor=(0.55, 0.015),
        fontsize=6.0,
        columnspacing=0.9,
        handlelength=1.1,
        handletextpad=0.35,
    )
    panel_label(ax, "d", -0.11, 1.02)
    fig.subplots_adjust(left=0.08, right=0.98, top=0.965, bottom=0.18)
    save(fig, "Figure3d_MAIN_TEXT_RECOMMENDED_full_universe_composition")
    data["analysis_total_genes"] = data.analysis.map(totals)
    write_tsv(data, "Figure3d_formal_plotdata.tsv")
    return data, checks


def redraw_3c() -> tuple[pd.DataFrame, list[dict[str, object]]]:
    source = pd.read_csv(MIRROR_SOURCE, sep="\t")
    checks: list[dict[str, object]] = []
    for (analysis, stage), part in source.groupby(["analysis", "stage"], observed=True):
        checks.append(
            {
                "panel": "3c",
                "check": f"{analysis}_{stage}_four_calls_sum_to_100",
                "observed": float(part.percentage.sum()),
                "expected": 100.0,
                "pass": bool(np.isclose(part.percentage.sum(), 100.0, atol=1e-8)),
            }
        )

    biased = source[
        source.ase_call.isin(["Allele_A_biased", "Allele_B_biased"])
    ].copy()
    biased["stage"] = pd.Categorical(biased.stage, STAGES, ordered=True)
    biased = biased.sort_values(["stage", "analysis", "ase_call"])
    biased["signed_genes"] = np.where(
        biased.ase_call.eq("Allele_A_biased"), biased.genes, -biased.genes
    )

    def values(analysis: str, call: str) -> np.ndarray:
        rows = biased[
            biased.analysis.eq(analysis) & biased.ase_call.eq(call)
        ].set_index("stage")
        return rows.reindex(STAGES).genes.to_numpy(float)

    x = np.arange(len(STAGES))
    width = 0.32
    fig, ax = plt.subplots(figsize=(6.8, 5.6))
    series = [
        ("TN", "Allele_A_biased", -width / 2, 1, TN_A, "TN A (Dura/TK-like)"),
        ("FL", "Allele_A_biased", width / 2, 1, FL_A, "FL A (Africa hap2)"),
        ("TN", "Allele_B_biased", -width / 2, -1, TN_B, "TN B (Pisifera/NS-like)"),
        ("FL", "Allele_B_biased", width / 2, -1, FL_B, "FL B (American hap1)"),
    ]
    for analysis, call, offset, sign, color, label in series:
        ax.bar(
            x + offset,
            sign * values(analysis, call),
            width=width,
            color=color,
            edgecolor="white",
            lw=0.4,
            label=label,
            zorder=3,
        )
    limit = float(np.ceil(max(abs(biased.signed_genes)) / 500) * 500 + 250)
    ax.axhline(0, color=DARK, lw=0.8)
    ax.axvline(12.5, color="#6F7478", lw=0.8, ls=(0, (4, 3)))
    ax.axvspan(-0.5, 12.5, color="#F2F7F5", zorder=-2)
    ax.axvspan(12.5, 18.5, color="#FFF7EA", zorder=-2)
    ax.text(6.0, limit * 0.94, "Development (0-185 d)", ha="center", fontsize=7.3)
    ax.text(15.5, limit * 0.94, "Postharvest (12-72 h)", ha="center", fontsize=7.3)
    ax.text(18.45, limit * 0.78, "Allele A biased", ha="right", fontsize=7, color="#59636A")
    ax.text(18.45, -limit * 0.80, "Allele B biased", ha="right", fontsize=7, color="#59636A")
    ax.set_ylim(-limit, limit)
    ax.set_xlim(-0.65, 18.65)
    ax.set_xticks(x, STAGES, rotation=45, ha="right")
    ax.yaxis.set_major_formatter(FuncFormatter(lambda value, _: f"{abs(int(value)):,}"))
    ax.set_ylabel("Number of ASE genes")
    clean(ax, "y")
    legend_handles, legend_labels = ax.get_legend_handles_labels()
    legend_order = [0, 2, 1, 3]
    ax.legend(
        [legend_handles[index] for index in legend_order],
        [legend_labels[index] for index in legend_order],
        frameon=False,
        ncol=2,
        loc="upper center",
        bbox_to_anchor=(0.5, 1.16),
        columnspacing=1.3,
        handlelength=1.4,
        fontsize=6.5,
    )
    panel_label(ax, "c", -0.065, 1.04)
    fig.subplots_adjust(left=0.085, right=0.99, top=0.78, bottom=0.22)
    save(fig, "Figure3c_revised_stage_ASE_counts")
    write_tsv(biased, "Figure3c_plotdata.tsv")
    return biased, checks


def flow_patch(
    ax: plt.Axes,
    x0: float,
    x1: float,
    y0a: float,
    y0b: float,
    y1a: float,
    y1b: float,
    color: str,
) -> None:
    control = 0.43 * (x1 - x0)
    vertices = [
        (x0, y0a),
        (x0 + control, y0a),
        (x1 - control, y1a),
        (x1, y1a),
        (x1, y1b),
        (x1 - control, y1b),
        (x0 + control, y0b),
        (x0, y0b),
        (x0, y0a),
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
            facecolor=color,
            edgecolor="none",
            alpha=0.38,
            zorder=1,
        )
    )


def node_intervals(totals: pd.Series, total: int) -> dict[str, tuple[float, float]]:
    gap = 0.025
    scale = (1.0 - gap * (len(CLASS_ORDER) - 1)) / total
    intervals: dict[str, tuple[float, float]] = {}
    top = 1.0
    for category in CLASS_ORDER:
        height = float(totals.get(category, 0)) * scale
        intervals[category] = (top - height, top)
        top -= height + gap
    return intervals


def redraw_3d_intersection_sensitivity() -> tuple[pd.DataFrame, list[dict[str, object]]]:
    data = pd.read_csv(FLOW_SOURCE, sep="\t")
    total = int(data.orthogroups.sum())
    display_order = list(reversed(CLASS_ORDER))
    left_totals = data.groupby("overall_class_FL").orthogroups.sum().reindex(display_order, fill_value=0)
    right_totals = data.groupby("overall_class_TN").orthogroups.sum().reindex(display_order, fill_value=0)

    def packed_intervals(totals: pd.Series) -> dict[str, tuple[float, float]]:
        intervals: dict[str, tuple[float, float]] = {}
        top = 1.0
        for category in display_order:
            height = float(totals[category]) / total
            intervals[category] = (top - height, top)
            top -= height
        return intervals

    left_nodes = packed_intervals(left_totals)
    right_nodes = packed_intervals(right_totals)
    scale = 1.0 / total

    left_cursor = {c: left_nodes[c][1] for c in CLASS_ORDER}
    right_cursor = {c: right_nodes[c][1] for c in CLASS_ORDER}
    segments: list[dict[str, object]] = []
    for source_class in display_order:
        for target_class in display_order:
            matches = data[
                data.overall_class_FL.eq(source_class)
                & data.overall_class_TN.eq(target_class)
            ]
            n = int(matches.orthogroups.iloc[0]) if not matches.empty else 0
            height = n * scale
            l_top = left_cursor[source_class]
            l_bottom = l_top - height
            r_top = right_cursor[target_class]
            r_bottom = r_top - height
            left_cursor[source_class] = l_bottom
            right_cursor[target_class] = r_bottom
            segments.append(
                {
                    "overall_class_FL": source_class,
                    "overall_class_TN": target_class,
                    "orthogroups": n,
                    "percentage_of_reviewed_orthogroups": 100 * n / total,
                    "FL_class_total": int(left_totals[source_class]),
                    "TN_class_total": int(right_totals[target_class]),
                    "FL_class_percentage": 100 * left_totals[source_class] / total,
                    "TN_class_percentage": 100 * right_totals[target_class] / total,
                    "left_bottom": l_bottom,
                    "left_top": l_top,
                    "right_bottom": r_bottom,
                    "right_top": r_top,
                }
            )

    # Match the narrow/tall slot used by panel d in the final composite. The
    # four classes are packed without side labels so the asset can replace the
    # Illustrator link directly; denominator details stay in the caption.
    fig, ax = plt.subplots(figsize=(2.15, 4.55), facecolor="white")
    fig.subplots_adjust(left=0.05, right=0.98, top=0.965, bottom=0.18)
    ax.set_xlim(0, 1)
    ax.set_ylim(-0.06, 1.02)
    ax.axis("off")
    for row in segments:
        if int(row["orthogroups"]) == 0:
            continue
        flow_patch(
            ax,
            0.25,
            0.75,
            float(row["left_bottom"]),
            float(row["left_top"]),
            float(row["right_bottom"]),
            float(row["right_top"]),
            CLASS_COLOR[str(row["overall_class_FL"])],
        )

    for category in display_order:
        for x0, nodes in [
            (0.12, left_nodes),
            (0.75, right_nodes),
        ]:
            bottom, top = nodes[category]
            ax.add_patch(
                Rectangle(
                    (x0, bottom),
                    0.13,
                    top - bottom,
                    facecolor=CLASS_COLOR[category],
                    edgecolor=DARK,
                    lw=0.65,
                    zorder=4,
                )
            )
    ax.text(0.185, -0.025, "FL", ha="center", va="top", fontsize=7.4)
    ax.text(0.815, -0.025, "TN", ha="center", va="top", fontsize=7.4)
    ax.text(
        0.5,
        1.012,
        "Shared orthogroups (n = 7,059)",
        ha="center",
        va="bottom",
        fontsize=5.0,
        color=DARK,
    )
    legend_order = ["HapDom", "NoDiff", "Sub", "NoASE"]
    handles = [
        Rectangle((0, 0), 1, 1, facecolor=CLASS_COLOR[c], edgecolor=DARK, lw=0.45)
        for c in legend_order
    ]
    fig.legend(
        handles,
        [CLASS_LABEL[c] for c in legend_order],
        frameon=False,
        ncol=2,
        loc="lower center",
        bbox_to_anchor=(0.55, 0.015),
        fontsize=5.7,
        columnspacing=0.75,
        handlelength=1.05,
        handletextpad=0.35,
    )
    fig.text(0.012, 0.992, "d", fontsize=12.5, fontweight="bold", ha="left", va="top")
    save(fig, "Figure3d_ALTERNATIVE_FLOW_intersection7059_alluvial")

    plotdata = pd.DataFrame(segments)
    write_tsv(plotdata, "Figure3d_intersection_plotdata.tsv")
    checks = [
        {
            "panel": "3d_alternative_flow",
            "check": "reviewed_one_to_one_orthogroup_total",
            "observed": total,
            "expected": 7059,
            "pass": total == 7059,
        },
        {
            "panel": "3d_alternative_flow",
            "check": "all_four_FL_classes_present",
            "observed": int((left_totals > 0).sum()),
            "expected": 4,
            "pass": bool((left_totals > 0).all()),
        },
        {
            "panel": "3d_alternative_flow",
            "check": "all_four_TN_classes_present",
            "observed": int((right_totals > 0).sum()),
            "expected": 4,
            "pass": bool((right_totals > 0).all()),
        },
    ]
    return plotdata, checks


def holm_adjust(pvalues: list[float]) -> np.ndarray:
    p = np.asarray(pvalues, dtype=float)
    order = np.argsort(p)
    adjusted = np.empty(len(p), dtype=float)
    running = 0.0
    for rank, index in enumerate(order):
        running = max(running, min(1.0, p[index] * (len(p) - rank)))
        adjusted[index] = running
    return adjusted


def bh_adjust(pvalues: list[float]) -> np.ndarray:
    p = np.asarray(pvalues, dtype=float)
    order = np.argsort(p)
    ranked = p[order]
    adjusted_ranked = np.minimum.accumulate(
        (ranked * len(p) / np.arange(1, len(p) + 1))[::-1]
    )[::-1]
    adjusted = np.empty(len(p), dtype=float)
    adjusted[order] = np.minimum(adjusted_ranked, 1.0)
    return adjusted


def compact_letters(
    significant: dict[tuple[int, int], bool], medians: list[float]
) -> list[str]:
    n_significant = sum(significant.values())
    ranked = sorted(range(3), key=lambda index: medians[index], reverse=True)
    if n_significant == 0:
        return ["a", "a", "a"]
    if n_significant == 3:
        result = [""] * 3
        for label, index in zip(["a", "b", "c"], ranked):
            result[index] = label
        return result
    if n_significant == 1:
        i, j = next(pair for pair, value in significant.items() if value)
        third = ({0, 1, 2} - {i, j}).pop()
        high, low = (i, j) if medians[i] >= medians[j] else (j, i)
        result = [""] * 3
        result[high], result[low], result[third] = "a", "b", "ab"
        return result
    i, j = next(pair for pair, value in significant.items() if not value)
    isolated = ({0, 1, 2} - {i, j}).pop()
    pair_letter = "b" if medians[isolated] >= np.mean([medians[i], medians[j]]) else "a"
    isolated_letter = "a" if pair_letter == "b" else "b"
    result = [pair_letter] * 3
    result[isolated] = isolated_letter
    return result


def p_text(pvalue: float) -> str:
    if pvalue < 0.001:
        exponent = int(math.floor(math.log10(pvalue)))
        coefficient = pvalue / (10**exponent)
        return f"P = {coefficient:.2f} x 10^{exponent}"
    return f"P = {pvalue:.3f}"


def redraw_3e() -> tuple[pd.DataFrame, pd.DataFrame, list[dict[str, object]]]:
    raw = pd.concat(
        [pd.read_csv(path, sep="\t") for path in SNP_PATHS.values()],
        ignore_index=True,
    )
    classes = ["HapDom", "NoDiff", "Sub"]
    regions = [
        ("upstream2kb_per_kb", "Upstream\n2 kb"),
        ("gene_body_per_kb", "Gene body"),
        ("downstream2kb_per_kb", "Downstream\n2 kb"),
    ]
    pairs = list(combinations(range(3), 2))
    long_rows: list[pd.DataFrame] = []
    test_rows: list[dict[str, object]] = []
    plot_groups: dict[tuple[str, str, str], np.ndarray] = {}
    for analysis in ["FL", "TN"]:
        subset = raw[raw.analysis.eq(analysis)]
        for region, region_label in regions:
            values_by_class: list[np.ndarray] = []
            for category in classes:
                gene_rows = subset.loc[
                    subset.overall_class.eq(category),
                    ["gene_id", "analysis", "overall_class", region],
                ].dropna()
                gene_rows = gene_rows.rename(columns={region: "SNPs_per_1000bp"})
                gene_rows["region"] = region
                gene_rows["region_label"] = region_label.replace("\n", " ")
                gene_rows["log10_1plus_SNPs_per_1000bp"] = np.log10(
                    1.0 + gene_rows.SNPs_per_1000bp.astype(float)
                )
                long_rows.append(gene_rows)
                values = gene_rows.log10_1plus_SNPs_per_1000bp.to_numpy(float)
                values_by_class.append(values)
                plot_groups[(analysis, region, category)] = values

            kw = stats.kruskal(*values_by_class)
            raw_p = [
                stats.mannwhitneyu(
                    values_by_class[i], values_by_class[j], alternative="two-sided"
                ).pvalue
                for i, j in pairs
            ]
            adjusted = holm_adjust(raw_p)
            significant = {
                pair: bool(pvalue < 0.05)
                for pair, pvalue in zip(pairs, adjusted)
            }
            medians = [float(np.median(values)) for values in values_by_class]
            letters = compact_letters(significant, medians)
            letter_text = ";".join(
                f"{category}:{letter}" for category, letter in zip(classes, letters)
            )
            test_rows.append(
                {
                    "analysis": analysis,
                    "region": region,
                    "test": "Kruskal-Wallis",
                    "comparison": "all three ASE classes",
                    "statistic": float(kw.statistic),
                    "pvalue": float(kw.pvalue),
                    "holm_pvalue": np.nan,
                    "significant_0.05": bool(kw.pvalue < 0.05),
                    "letters": letter_text,
                }
            )
            for (i, j), pvalue, adjusted_p in zip(pairs, raw_p, adjusted):
                test_rows.append(
                    {
                        "analysis": analysis,
                        "region": region,
                        "test": "Mann-Whitney U",
                        "comparison": f"{classes[i]} vs {classes[j]}",
                        "statistic": np.nan,
                        "pvalue": float(pvalue),
                        "holm_pvalue": float(adjusted_p),
                        "significant_0.05": bool(adjusted_p < 0.05),
                        "letters": letter_text,
                    }
                )

    plotdata = pd.concat(long_rows, ignore_index=True)
    tests = pd.DataFrame(test_rows)
    write_tsv(plotdata, "Figure3e_plotdata.tsv")
    write_tsv(tests, "Figure3e_tests.tsv")

    fig, axes = plt.subplots(1, 2, figsize=(8.0, 3.9), sharey=True)
    positions = [0, 1, 2, 4, 5, 6, 8, 9, 10]
    colors = [CLASS_COLOR[c] for _ in regions for c in classes]
    for axis, analysis in zip(axes, ["FL", "TN"]):
        data_arrays = [
            plot_groups[(analysis, region, category)]
            for region, _ in regions
            for category in classes
        ]
        boxplot = axis.boxplot(
            data_arrays,
            positions=positions,
            widths=0.68,
            showfliers=False,
            patch_artist=True,
            medianprops={"color": DARK, "lw": 0.8},
            whiskerprops={"color": "#626A70", "lw": 0.6},
            capprops={"color": "#626A70", "lw": 0.6},
        )
        for box, color in zip(boxplot["boxes"], colors):
            box.set_facecolor("white")
            box.set_edgecolor(color)
            box.set_linewidth(1.1)
        for region_index, (region, _) in enumerate(regions):
            region_tests = tests[
                tests.analysis.eq(analysis) & tests.region.eq(region)
            ]
            kw = region_tests[region_tests.test.eq("Kruskal-Wallis")].iloc[0]
            axis.text(
                region_index * 4 + 1,
                2.55,
                p_text(float(kw.pvalue)),
                ha="center",
                va="bottom",
                fontsize=6.5,
            )
            letters = {
                item.split(":")[0]: item.split(":")[1]
                for item in str(kw.letters).split(";")
            }
            for class_index, category in enumerate(classes):
                axis.text(
                    region_index * 4 + class_index,
                    2.32,
                    letters[category],
                    ha="center",
                    va="bottom",
                    fontsize=7.5,
                    fontweight="bold",
                )
        axis.set_xticks([1, 5, 9], [label for _, label in regions])
        axis.set_title(analysis, fontweight="bold")
        axis.set_ylim(-0.08, 2.73)
        clean(axis, "y")
    axes[0].set_ylabel("log10(1 + SNPs per 1,000 bp)")
    handles = [
        Rectangle((0, 0), 1, 1, facecolor="white", edgecolor=CLASS_COLOR[c], lw=1.1)
        for c in classes
    ]
    fig.legend(
        handles,
        [CLASS_LABEL[c] for c in classes],
        frameon=False,
        ncol=3,
        loc="upper center",
        bbox_to_anchor=(0.55, 0.99),
    )
    axes[0].text(
        0.0,
        1.18,
        "Diagnostic SNP density by ASE class",
        transform=axes[0].transAxes,
        fontweight="bold",
        fontsize=9.5,
    )
    axes[0].text(
        0.0,
        -0.25,
        "Above regions: Kruskal-Wallis P; letters: Holm-adjusted pairwise Mann-Whitney U.",
        transform=axes[0].transAxes,
        fontsize=6.3,
        color="#59636A",
    )
    panel_label(axes[0], "e", -0.18, 1.17)
    fig.subplots_adjust(left=0.10, right=0.99, top=0.74, bottom=0.23, wspace=0.11)
    save(fig, "Figure3e_AUDIT_ONLY_SNP_density_FL_TN")

    # The final main-text panel contains FL only. Use the final composite's
    # raw-density visual language; rank-test results are unchanged by the
    # monotonic log transform used in the audit companion above.
    fig, ax = plt.subplots(figsize=(3.80, 3.52), facecolor="white")
    data_arrays = [
        plotdata.loc[
            plotdata.analysis.eq("FL")
            & plotdata.region.eq(region)
            & plotdata.overall_class.eq(category),
            "SNPs_per_1000bp",
        ].to_numpy(float)
        for region, _ in regions
        for category in classes
    ]
    boxplot = ax.boxplot(
        data_arrays,
        positions=positions,
        widths=0.68,
        showfliers=True,
        patch_artist=True,
        medianprops={"color": DARK, "lw": 0.7},
        whiskerprops={"color": DARK, "lw": 0.55},
        capprops={"color": DARK, "lw": 0.55},
        flierprops={"marker": "o", "markersize": 1.6, "markerfacecolor": "white", "markeredgewidth": 0.45},
    )
    final_colors = {"HapDom": "#7FAFDC", "NoDiff": "#EC6A6B", "Sub": "#EFD78B"}
    final_color_sequence = [final_colors[c] for _ in regions for c in classes]
    for index, (box, color) in enumerate(zip(boxplot["boxes"], final_color_sequence)):
        box.set_facecolor("white")
        box.set_edgecolor(color)
        box.set_linewidth(0.75)
        boxplot["medians"][index].set_color(color)
        for artist in boxplot["whiskers"][2 * index : 2 * index + 2]:
            artist.set_color(color)
        for artist in boxplot["caps"][2 * index : 2 * index + 2]:
            artist.set_color(color)
        boxplot["fliers"][index].set_markeredgecolor(color)
    for region_index, (region, _) in enumerate(regions):
        region_tests = tests[
            tests.analysis.eq("FL") & tests.region.eq(region)
        ]
        kw = region_tests[region_tests.test.eq("Kruskal-Wallis")].iloc[0]
        ax.text(
            region_index * 4 + 1,
            181.0,
            p_text(float(kw.pvalue)),
            ha="center",
            va="bottom",
            fontsize=5.7,
        )
        letters = {
            item.split(":")[0]: item.split(":")[1]
            for item in str(kw.letters).split(";")
        }
        for class_index, category in enumerate(classes):
            ax.text(
                region_index * 4 + class_index,
                168.0,
                letters[category],
                ha="center",
                va="bottom",
                fontsize=6.8,
                fontweight="bold",
            )
    ax.set_xticks(
        [1, 5, 9],
        ["Upstream\n2 kb", "Gene", "Downstream\n2 kb"],
    )
    ax.set_ylim(-4.0, 190.0)
    ax.set_yticks([0, 40, 80, 120, 160])
    ax.set_ylabel("SNPs per 1,000 bp")
    ax.set_title("FL", loc="left", fontsize=8.0, fontweight="bold", pad=8)
    clean(ax)
    handles = [
        Rectangle((0, 0), 1, 1, facecolor=final_colors[c], edgecolor="none")
        for c in classes
    ]
    fig.legend(
        handles,
        [CLASS_LABEL[c] for c in classes],
        frameon=False,
        ncol=3,
        loc="upper center",
        bbox_to_anchor=(0.61, 0.995),
        fontsize=6.2,
        columnspacing=0.8,
        handlelength=1.1,
    )
    fig.text(0.012, 0.992, "e", fontsize=12.5, fontweight="bold", ha="left", va="top")
    fig.subplots_adjust(left=0.20, right=0.98, top=0.79, bottom=0.16)
    save(fig, "Figure3e_revised_FL_SNP_density", tight=False)

    reference = pd.read_csv(TEST_REFERENCE, sep="\t")
    checks: list[dict[str, object]] = []
    for _, expected in reference.iterrows():
        observed_rows = tests[
            tests.analysis.eq(expected.analysis)
            & tests.region.eq(expected.region)
            & tests.test.eq(expected.test)
            & tests.comparison.eq(expected.comparison)
        ]
        if observed_rows.empty:
            checks.append(
                {
                    "panel": "3e",
                    "check": f"formal_test_reproduced:{expected.analysis}:{expected.region}:{expected.test}:{expected.comparison}",
                    "observed": "missing",
                    "expected": float(expected.pvalue),
                    "pass": False,
                }
            )
            continue
        observed = float(observed_rows.iloc[0].pvalue)
        checks.append(
            {
                "panel": "3e",
                "check": f"formal_test_reproduced:{expected.analysis}:{expected.region}:{expected.test}:{expected.comparison}",
                "observed": observed,
                "expected": float(expected.pvalue),
                "pass": bool(np.isclose(observed, float(expected.pvalue), rtol=1e-10, atol=1e-15)),
            }
        )
    formal_fl = tests[
        tests.analysis.eq("FL") & tests.test.eq("Kruskal-Wallis")
    ].set_index("region").pvalue
    for region, expected in {
        "upstream2kb_per_kb": 0.301805654312622,
        "gene_body_per_kb": 0.7758834655889496,
        "downstream2kb_per_kb": 0.24636421212683132,
    }.items():
        observed = float(formal_fl.loc[region])
        checks.append(
            {
                "panel": "3e",
                "check": f"formal_main_text_P_is_FL:{region}",
                "observed": observed,
                "expected": expected,
                "pass": bool(np.isclose(observed, expected, rtol=1e-10, atol=1e-15)),
            }
        )
    return plotdata, tests, checks


def redraw_3g() -> tuple[pd.DataFrame, list[dict[str, object]]]:
    raw = pd.read_csv(
        MODE_DETAIL_SOURCE,
        sep="\t",
        usecols=["gene_id_africa_hap2", "stage", "stage_index", "model1"],
    )
    key = ["gene_id_africa_hap2", "stage"]
    duplicate_rows = int(raw.duplicated(key).sum())
    conflicting_keys = int(
        (raw.groupby(key, observed=True).model1.nunique(dropna=False) > 1).sum()
    )
    if conflicting_keys:
        raise RuntimeError(
            f"Cannot deduplicate {MODE_DETAIL_SOURCE}: {conflicting_keys} gene-stage keys have conflicting model1 labels"
        )
    unique = raw.drop_duplicates(key).copy()
    source = (
        unique.groupby(["stage_index", "stage", "model1"], observed=True)
        .size()
        .rename("genes")
        .reset_index()
        .rename(columns={"model1": "mode"})
    )
    totals = (
        unique.groupby(["stage_index", "stage"], observed=True)
        .size()
        .rename("stage_total")
        .reset_index()
    )
    source = source.merge(totals, on=["stage_index", "stage"], how="left")
    source["percentage"] = 100.0 * source.genes / source.stage_total
    source["inference_level"] = "exploratory_parent_n1_deduplicated_gene_stage"

    modes = [f"M{i}" for i in range(1, 13)] + ["Conserved"]
    data = source[source["mode"].isin(modes)].copy()
    data["stage"] = pd.Categorical(data.stage, STAGES, ordered=True)
    data["mode"] = pd.Categorical(data["mode"], modes, ordered=True)
    data = data.sort_values(["stage", "mode"]).reset_index(drop=True)
    counts = data.pivot(index="stage", columns="mode", values="genes").reindex(
        index=STAGES, columns=modes
    )
    percentages = data.pivot(index="stage", columns="mode", values="percentage").reindex(
        index=STAGES, columns=modes
    )
    stage_totals = data.groupby("stage", observed=True).stage_total.first().reindex(STAGES)
    checks: list[dict[str, object]] = [
        {
            "panel": "3g",
            "check": "source_duplicate_gene_stage_rows_detected",
            "observed": duplicate_rows,
            "expected": 79384,
            "pass": duplicate_rows == 79384,
        },
        {
            "panel": "3g",
            "check": "duplicate_gene_stage_mode_conflicts",
            "observed": conflicting_keys,
            "expected": 0,
            "pass": conflicting_keys == 0,
        },
        {
            "panel": "3g",
            "check": "unique_gene_stage_rows_after_deduplication",
            "observed": len(unique),
            "expected": 343802,
            "pass": len(unique) == 343802,
        },
    ]
    for stage in STAGES:
        observed_count = int(counts.loc[stage].sum())
        expected_count = int(stage_totals.loc[stage])
        observed_percentage = float(percentages.loc[stage].sum())
        checks.extend(
            [
                {
                    "panel": "3g",
                    "check": f"{stage}_M1_M12_plus_Conserved_count_equals_stage_total",
                    "observed": observed_count,
                    "expected": expected_count,
                    "pass": observed_count == expected_count,
                },
                {
                    "panel": "3g",
                    "check": f"{stage}_M1_M12_plus_Conserved_percentage_equals_100",
                    "observed": observed_percentage,
                    "expected": 100.0,
                    "pass": bool(np.isclose(observed_percentage, 100.0, atol=5e-4)),
                },
            ]
        )

    array = percentages.to_numpy(float)
    vmin = float(np.nanmin(array))
    vmax = float(np.nanmax(array))
    center = float(np.nanmedian(array))
    fig, ax = plt.subplots(figsize=(6.5, 3.3))
    image = ax.imshow(
        array,
        cmap=DPM_HEATMAP_CMAP,
        norm=TwoSlopeNorm(vmin=vmin, vcenter=center, vmax=vmax),
        aspect="auto",
    )
    for row in range(len(STAGES)):
        for column in range(len(modes)):
            value = int(counts.iloc[row, column])
            percent = float(array[row, column])
            ax.text(
                column,
                row,
                f"{value:,}",
                ha="center",
                va="center",
                fontsize=4.2,
                color=(
                    "white"
                    if abs(percent - center) > 0.36 * (vmax - vmin)
                    else DARK
                ),
            )
    ax.set_xticks(range(len(modes)), modes, fontweight="bold", fontsize=5.8)
    ax.xaxis.tick_top()
    ax.tick_params(top=True, bottom=False, labeltop=True, labelbottom=False)
    ax.set_yticks(range(len(STAGES)), STAGES)
    ax.axvline(11.5, color="white", lw=2.6)
    ax.axhline(12.5, color=DARK, lw=0.8, ls="--")
    for spine in ax.spines.values():
        spine.set_linewidth(0.7)
        spine.set_color("#777777")
    colorbar = fig.colorbar(image, ax=ax, pad=0.014, fraction=0.028)
    colorbar.set_label("Within-stage genes (%)")
    panel_label(ax, "g", -0.09, 1.10)
    fig.subplots_adjust(left=0.09, right=0.94, top=0.85, bottom=0.07)
    save(fig, "Figure3g_revised_modes_with_Conserved")
    data["display_order"] = data["mode"].cat.codes + 1
    write_tsv(data, "Figure3g_plotdata.tsv")

    legacy = pd.read_csv(MODE_LEGACY_SUMMARY, sep="\t")
    legacy = legacy[legacy.model.eq("model1")][
        ["stage", "mode", "genes", "stage_total", "percentage"]
    ].rename(
        columns={
            "genes": "legacy_rows",
            "stage_total": "legacy_stage_rows",
            "percentage": "legacy_percentage",
        }
    )
    comparison = data[
        ["stage", "mode", "genes", "stage_total", "percentage"]
    ].copy()
    comparison["stage"] = comparison.stage.astype(str)
    comparison["mode"] = comparison["mode"].astype(str)
    comparison = comparison.merge(legacy, on=["stage", "mode"], how="outer")
    comparison["removed_duplicate_rows"] = (
        comparison.legacy_rows - comparison.genes
    )
    comparison["percentage_point_change"] = (
        comparison.percentage - comparison.legacy_percentage
    )
    write_tsv(comparison, "Figure3g_legacy_vs_deduplicated.tsv")

    duplicate_audit = (
        raw.assign(is_duplicate=raw.duplicated(key, keep="first"))
        .groupby(["stage_index", "stage"], observed=True)
        .agg(
            raw_rows=("gene_id_africa_hap2", "size"),
            exact_duplicate_rows=("is_duplicate", "sum"),
            unique_gene_stage_rows=("gene_id_africa_hap2", "nunique"),
        )
        .reset_index()
    )
    duplicate_audit["duplicate_row_percentage"] = (
        100.0 * duplicate_audit.exact_duplicate_rows / duplicate_audit.raw_rows
    )
    write_tsv(duplicate_audit, "Figure3g_duplicate_audit.tsv")
    return data, checks


def _contingency_residuals(counts: np.ndarray) -> tuple[float, np.ndarray]:
    counts = np.asarray(counts, dtype=float)
    expected = np.outer(counts.sum(axis=1), counts.sum(axis=0)) / counts.sum()
    residuals = np.divide(
        counts - expected,
        np.sqrt(expected),
        out=np.zeros_like(expected),
        where=expected > 0,
    )
    return float(np.square(residuals).sum()), residuals


def redraw_3j() -> tuple[pd.DataFrame, list[dict[str, object]]]:
    reg = pd.read_csv(
        REGULATORY_SOURCE,
        sep="\t",
        usecols=["gene_africa", "stage", "regulatory_class"],
    )
    inheritance_raw = pd.read_csv(
        MODE_DETAIL_SOURCE,
        sep="\t",
        usecols=["gene_id_africa_hap2", "stage", "model3"],
    )
    inheritance_key = ["gene_id_africa_hap2", "stage"]
    upstream_duplicate_rows = int(inheritance_raw.duplicated(inheritance_key).sum())
    inheritance_conflicts = int(
        (
            inheritance_raw.groupby(inheritance_key, observed=True)
            .model3.nunique(dropna=False)
            > 1
        ).sum()
    )
    if inheritance_conflicts:
        raise RuntimeError(
            f"Cannot deduplicate inheritance source: {inheritance_conflicts} conflicting gene-stage labels"
        )
    inheritance = inheritance_raw.drop_duplicates(inheritance_key)
    merged_raw = reg.merge(
        inheritance_raw,
        left_on=["gene_africa", "stage"],
        right_on=["gene_id_africa_hap2", "stage"],
        how="inner",
    )
    merged_duplicate_rows = int(
        merged_raw.duplicated(["gene_africa", "stage"]).sum()
    )
    merged = reg.merge(
        inheritance,
        left_on=["gene_africa", "stage"],
        right_on=["gene_id_africa_hap2", "stage"],
        how="inner",
        validate="one_to_one",
    )
    data = merged[
        merged.model3.isin(INHERITANCE_ORDER)
        & merged.regulatory_class.isin(REG_ORDER)
    ].copy()

    obs = (
        pd.crosstab(data.model3, data.regulatory_class)
        .reindex(index=INHERITANCE_ORDER, columns=REG_ORDER, fill_value=0)
        .astype(int)
    )
    row_totals = obs.sum(axis=1)
    all_genes = sorted(data.gene_africa.unique())
    gene_lookup = {gene: index for index, gene in enumerate(all_genes)}
    inheritance_lookup = {value: index for index, value in enumerate(INHERITANCE_ORDER)}
    reg_lookup = {value: index for index, value in enumerate(REG_ORDER)}
    gene_counts = np.zeros((len(all_genes), len(INHERITANCE_ORDER) * len(REG_ORDER)), dtype=np.int32)
    gene_codes = data.gene_africa.map(gene_lookup).to_numpy(int)
    cell_codes = (
        data.model3.map(inheritance_lookup).to_numpy(int) * len(REG_ORDER)
        + data.regulatory_class.map(reg_lookup).to_numpy(int)
    )
    np.add.at(gene_counts, (gene_codes, cell_codes), 1)

    bootstrap_replicates = 2000
    rng_boot = np.random.default_rng(20260831)
    boot_percentages = np.empty(
        (bootstrap_replicates, len(INHERITANCE_ORDER), len(REG_ORDER)),
        dtype=float,
    )
    probabilities = np.full(len(all_genes), 1.0 / len(all_genes))
    batch_size = 100
    for start in range(0, bootstrap_replicates, batch_size):
        stop = min(start + batch_size, bootstrap_replicates)
        weights = rng_boot.multinomial(
            len(all_genes), probabilities, size=stop - start
        )
        sampled = (weights @ gene_counts).reshape(
            stop - start, len(INHERITANCE_ORDER), len(REG_ORDER)
        )
        denominators = sampled.sum(axis=2, keepdims=True)
        boot_percentages[start:stop] = np.divide(
            100.0 * sampled,
            denominators,
            out=np.full(sampled.shape, np.nan, dtype=float),
            where=denominators > 0,
        )
    ci_low = np.nanpercentile(boot_percentages, 2.5, axis=0)
    ci_high = np.nanpercentile(boot_percentages, 97.5, axis=0)

    # Permute complete inheritance trajectories only among genes with identical
    # observed-stage coverage. This keeps every gene's repeated stages together
    # and preserves stage-specific missingness. Singleton coverage strata are
    # excluded from the test because they are not exchangeable.
    stage_lookup = {stage: index for index, stage in enumerate(STAGES)}
    permutation_genes = sorted(data.gene_africa.unique())
    permutation_gene_lookup = {
        gene: index for index, gene in enumerate(permutation_genes)
    }
    model_matrix = np.full((len(permutation_genes), len(STAGES)), -1, dtype=np.int8)
    reg_matrix = np.full_like(model_matrix, -1)
    row_gene = data.gene_africa.map(permutation_gene_lookup).to_numpy(int)
    row_stage = data.stage.map(stage_lookup).to_numpy(int)
    model_matrix[row_gene, row_stage] = data.model3.map(inheritance_lookup).to_numpy(np.int8)
    reg_matrix[row_gene, row_stage] = data.regulatory_class.map(reg_lookup).to_numpy(np.int8)
    coverage_groups: dict[bytes, list[int]] = {}
    for index, mask in enumerate(model_matrix >= 0):
        coverage_groups.setdefault(np.packbits(mask).tobytes(), []).append(index)
    exchangeable_groups = [
        np.asarray(indices, dtype=int)
        for indices in coverage_groups.values()
        if len(indices) >= 2
    ]
    exchangeable_indices = np.concatenate(exchangeable_groups)
    obs_model = model_matrix[exchangeable_indices].ravel()
    obs_reg = reg_matrix[exchangeable_indices].ravel()
    valid = (obs_model >= 0) & (obs_reg >= 0)
    observed_permutation_counts = np.bincount(
        obs_model[valid] * len(REG_ORDER) + obs_reg[valid],
        minlength=len(INHERITANCE_ORDER) * len(REG_ORDER),
    ).reshape(len(INHERITANCE_ORDER), len(REG_ORDER))
    observed_statistic, observed_residuals = _contingency_residuals(
        observed_permutation_counts
    )
    permutation_replicates = 1999
    rng_perm = np.random.default_rng(20260832)
    global_exceedances = 0
    cell_exceedances = np.zeros_like(observed_residuals, dtype=int)
    for _ in range(permutation_replicates):
        permuted_model = model_matrix.copy()
        for indices in exchangeable_groups:
            permuted_model[indices] = model_matrix[rng_perm.permutation(indices)]
        flat_model = permuted_model[exchangeable_indices].ravel()
        valid = (flat_model >= 0) & (obs_reg >= 0)
        counts = np.bincount(
            flat_model[valid] * len(REG_ORDER) + obs_reg[valid],
            minlength=len(INHERITANCE_ORDER) * len(REG_ORDER),
        ).reshape(len(INHERITANCE_ORDER), len(REG_ORDER))
        statistic, residuals = _contingency_residuals(counts)
        global_exceedances += int(statistic >= observed_statistic - 1e-12)
        cell_exceedances += np.abs(residuals) >= np.abs(observed_residuals) - 1e-12
    global_p = (global_exceedances + 1) / (permutation_replicates + 1)
    cell_p = (cell_exceedances.ravel() + 1) / (permutation_replicates + 1)
    cell_holm = holm_adjust(cell_p.tolist())
    cell_bh = bh_adjust(cell_p.tolist())

    rows: list[dict[str, object]] = []
    for i, inheritance_class in enumerate(INHERITANCE_ORDER):
        for j, regulatory_class in enumerate(REG_ORDER):
            cell = data[
                data.model3.eq(inheritance_class)
                & data.regulatory_class.eq(regulatory_class)
            ]
            rows.append(
                {
                    "inheritance_class": inheritance_class,
                    "regulatory_class": regulatory_class,
                    "deduplicated_gene_stage_observations": int(obs.iloc[i, j]),
                    "unique_genes_in_cell": int(cell.gene_africa.nunique()),
                    "inheritance_row_gene_stage_total": int(row_totals.iloc[i]),
                    "row_percentage": 100.0 * obs.iloc[i, j] / row_totals.iloc[i],
                    "gene_cluster_bootstrap_ci_low": float(ci_low[i, j]),
                    "gene_cluster_bootstrap_ci_high": float(ci_high[i, j]),
                    "trajectory_permutation_residual": float(observed_residuals[i, j]),
                    "trajectory_permutation_p": float(cell_p[i * len(REG_ORDER) + j]),
                    "trajectory_permutation_holm_p": float(cell_holm[i * len(REG_ORDER) + j]),
                    "trajectory_permutation_bh_q": float(cell_bh[i * len(REG_ORDER) + j]),
                    "inference_level": "exploratory_parent_n1_gene_clustered",
                }
            )
    plotdata = pd.DataFrame(rows)
    write_tsv(plotdata, "Figure3j_plotdata_gene_clustered.tsv")

    permutation_rows = [
        {
            "scope": "global_3x7",
            "inheritance_class": "ALL",
            "regulatory_class": "ALL",
            "statistic": observed_statistic,
            "empirical_p": global_p,
            "holm_p": np.nan,
            "bh_q": np.nan,
            "permutations": permutation_replicates,
            "exchangeable_genes": len(exchangeable_indices),
            "all_genes": len(permutation_genes),
            "exchangeable_stage_coverage_strata": len(exchangeable_groups),
        }
    ]
    for i, inheritance_class in enumerate(INHERITANCE_ORDER):
        for j, regulatory_class in enumerate(REG_ORDER):
            index = i * len(REG_ORDER) + j
            permutation_rows.append(
                {
                    "scope": "cell_two_sided_residual",
                    "inheritance_class": inheritance_class,
                    "regulatory_class": regulatory_class,
                    "statistic": observed_residuals[i, j],
                    "empirical_p": cell_p[index],
                    "holm_p": cell_holm[index],
                    "bh_q": cell_bh[index],
                    "permutations": permutation_replicates,
                    "exchangeable_genes": len(exchangeable_indices),
                    "all_genes": len(permutation_genes),
                    "exchangeable_stage_coverage_strata": len(exchangeable_groups),
                }
            )
    permutation_table = pd.DataFrame(permutation_rows)
    write_tsv(permutation_table, "Figure3j_gene_trajectory_permutation_tests.tsv")

    legacy = pd.read_csv(LEGACY_INHERITANCE_SOURCE, sep="\t")
    comparison = plotdata.merge(
        legacy[
            [
                "inheritance_class",
                "regulatory_class",
                "genes",
                "row_percentage",
                "fisher_p",
                "fisher_q_BH",
            ]
        ].rename(
            columns={
                "genes": "legacy_duplicated_gene_stage_rows",
                "row_percentage": "legacy_row_percentage",
                "fisher_p": "legacy_invalid_fisher_p",
                "fisher_q_BH": "legacy_invalid_fisher_q_BH",
            }
        ),
        on=["inheritance_class", "regulatory_class"],
        how="left",
    )
    comparison["removed_exact_duplicate_rows"] = (
        comparison.legacy_duplicated_gene_stage_rows
        - comparison.deduplicated_gene_stage_observations
    )
    write_tsv(comparison, "Figure3j_legacy_vs_gene_clustered.tsv")

    percentage = (
        plotdata.pivot(
            index="inheritance_class",
            columns="regulatory_class",
            values="row_percentage",
        )
        .reindex(index=INHERITANCE_ORDER, columns=REG_ORDER)
        .to_numpy(float)
    )
    low = (
        plotdata.pivot(
            index="inheritance_class",
            columns="regulatory_class",
            values="gene_cluster_bootstrap_ci_low",
        )
        .reindex(index=INHERITANCE_ORDER, columns=REG_ORDER)
        .to_numpy(float)
    )
    high = (
        plotdata.pivot(
            index="inheritance_class",
            columns="regulatory_class",
            values="gene_cluster_bootstrap_ci_high",
        )
        .reindex(index=INHERITANCE_ORDER, columns=REG_ORDER)
        .to_numpy(float)
    )
    residual = (
        plotdata.pivot(
            index="inheritance_class",
            columns="regulatory_class",
            values="trajectory_permutation_residual",
        )
        .reindex(index=INHERITANCE_ORDER, columns=REG_ORDER)
        .to_numpy(float)
    )
    vmax = max(
        75.0,
        25.0 * math.ceil(float(np.nanmax(np.abs(residual))) / 25.0),
    )
    norm = TwoSlopeNorm(vmin=-vmax, vcenter=0.0, vmax=vmax)
    fig = plt.figure(figsize=(3.30, 2.90), facecolor="white")
    ax = fig.add_axes([0.13, 0.27, 0.80, 0.53])
    x_edges = np.arange(len(REG_ORDER) + 1) - 0.5
    y_edges = np.arange(len(INHERITANCE_ORDER) + 1) - 0.5
    image = ax.pcolormesh(
        x_edges,
        y_edges,
        residual,
        cmap=DPM_HEATMAP_CMAP,
        norm=norm,
        shading="flat",
        edgecolors="none",
        rasterized=False,
    )
    ax.set_xlim(-0.5, len(REG_ORDER) - 0.5)
    ax.set_ylim(len(INHERITANCE_ORDER) - 0.5, -0.5)
    for i in range(len(INHERITANCE_ORDER)):
        for j in range(len(REG_ORDER)):
            ax.text(
                j,
                i,
                f"{percentage[i, j]:.1f}%",
                ha="center",
                va="center",
                fontsize=5.1,
                color="white" if abs(residual[i, j]) > 0.55 * vmax else "#222222",
                linespacing=1.05,
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
    fig.text(0.012, 0.985, "j", ha="left", va="top", fontsize=12.5, fontweight="bold")
    save(fig, "Figure3j_revised_gene_clustered_inheritance_regulatory")

    checks = [
        {
            "panel": "3j",
            "check": "upstream_exact_duplicate_gene_stage_rows",
            "observed": upstream_duplicate_rows,
            "expected": 79384,
            "pass": upstream_duplicate_rows == 79384,
        },
        {
            "panel": "3j",
            "check": "inheritance_duplicate_label_conflicts",
            "observed": inheritance_conflicts,
            "expected": 0,
            "pass": inheritance_conflicts == 0,
        },
        {
            "panel": "3j",
            "check": "duplicate_rows_in_original_3j_merge",
            "observed": merged_duplicate_rows,
            "expected": 25905,
            "pass": merged_duplicate_rows == 25905,
        },
        {
            "panel": "3j",
            "check": "deduplicated_gene_stage_keys_are_unique",
            "observed": int(merged.duplicated(["gene_africa", "stage"]).sum()),
            "expected": 0,
            "pass": not merged.duplicated(["gene_africa", "stage"]).any(),
        },
        {
            "panel": "3j",
            "check": "gene_clusters_in_descriptive_analysis",
            "observed": len(all_genes),
            "expected": 9009,
            "pass": len(all_genes) == 9009,
        },
        {
            "panel": "3j",
            "check": "exchangeable_gene_coverage_fraction",
            "observed": len(exchangeable_indices) / len(permutation_genes),
            "expected": ">=0.80",
            "pass": len(exchangeable_indices) / len(permutation_genes) >= 0.80,
        },
        {
            "panel": "3j",
            "check": "cluster_bootstrap_intervals_are_finite",
            "observed": int(np.isfinite(low).sum() + np.isfinite(high).sum()),
            "expected": 2 * len(INHERITANCE_ORDER) * len(REG_ORDER),
            "pass": bool(np.isfinite(low).all() and np.isfinite(high).all()),
        },
    ]
    return plotdata, checks


def redraw_3k() -> tuple[pd.DataFrame, list[dict[str, object]]]:
    raw = pd.read_csv(TRAIT_SOURCE, sep="\t")
    key = ["analysis", "trait_module", "stage_group", "gene_id", "stage"]
    duplicate_rows = int(raw.duplicated(key).sum())
    conflict_columns = ["eligible", "robust_ase", "log2_allele_ratio", "ase_call"]
    conflict_keys = 0
    for column in conflict_columns:
        conflict_keys += int(
            (raw.groupby(key, observed=True)[column].nunique(dropna=False) > 1).sum()
        )
    if conflict_keys:
        raise RuntimeError(
            f"Cannot deduplicate {TRAIT_SOURCE}: {conflict_keys} conflicting key-column combinations"
        )
    data = raw.drop_duplicates(key).copy()
    for column in ["eligible", "robust_ase"]:
        if data[column].dtype != bool:
            data[column] = data[column].astype(str).str.lower().isin(
                ["true", "1", "yes"]
            )
    eligible = data[data.eligible].copy()

    rows: list[dict[str, object]] = []
    for (analysis, module, phase), group in eligible.groupby(
        ["analysis", "trait_module", "stage_group"], observed=True
    ):
        per_gene = group.groupby("gene_id", observed=True).agg(
            eligible_stages=("stage", "nunique"),
            any_robust_ASE=("robust_ase", "any"),
        )
        robust = group[group.robust_ase]
        robust_gene_medians = robust.groupby("gene_id", observed=True)[
            "log2_allele_ratio"
        ].median()
        eligible_genes = int(len(per_gene))
        robust_genes = int(per_gene.any_robust_ASE.sum())
        rows.append(
            {
                "analysis": analysis,
                "trait_module": module,
                "stage_group": phase,
                "eligible_unique_genes": eligible_genes,
                "robust_ASE_unique_genes": robust_genes,
                "robust_ASE_unique_gene_percentage": 100.0
                * robust_genes
                / eligible_genes,
                "median_of_within_gene_robust_log2_ratios": float(
                    robust_gene_medians.median()
                ),
                "eligible_unique_gene_stage_rows": int(len(group)),
                "robust_unique_gene_stage_rows": int(group.robust_ase.sum()),
                "size_encoding": "percentage_of_eligible_unique_genes_with_at_least_one_robust_ASE_stage_in_phase",
                "colour_encoding": "median_across_genes_of_each_genes_phase_median_robust_log2_A_over_B",
                "allele_A_definition": "Africa hap2" if analysis == "FL" else "Dura/TK-like",
                "allele_B_definition": "American hap1" if analysis == "FL" else "Pisifera/NS-like",
            }
        )
    summary = pd.DataFrame(rows)
    summary["trait_module"] = pd.Categorical(
        summary.trait_module, MODULE_ORDER, ordered=True
    )
    summary["stage_group"] = pd.Categorical(
        summary.stage_group, PHASE_ORDER, ordered=True
    )
    summary = summary.sort_values(
        ["stage_group", "trait_module", "analysis"]
    ).reset_index(drop=True)
    write_tsv(summary, "Figure3k_plotdata_unique_genes.tsv")

    legacy = pd.read_csv(TRAIT_LEGACY_SUMMARY, sep="\t")
    comparison = summary.merge(
        legacy[
            [
                "analysis",
                "trait_module",
                "stage_group",
                "tested_gene_stage_rows",
                "tested_genes",
                "robust_ASE_rows",
                "robust_ASE_percentage",
                "robust_median_log2_ratio",
            ]
        ],
        on=["analysis", "trait_module", "stage_group"],
        how="left",
    )
    comparison["percentage_point_change"] = (
        comparison.robust_ASE_unique_gene_percentage
        - comparison.robust_ASE_percentage
    )
    comparison["colour_change"] = (
        comparison.median_of_within_gene_robust_log2_ratios
        - comparison.robust_median_log2_ratio
    )
    write_tsv(comparison, "Figure3k_legacy_rows_vs_unique_genes.tsv")

    duplicate_audit = pd.DataFrame(
        [
            {
                "raw_rows": len(raw),
                "unique_analysis_module_gene_stage_rows": len(data),
                "removed_exact_duplicate_rows": duplicate_rows,
                "conflicting_key_column_combinations": conflict_keys,
            }
        ]
    )
    write_tsv(duplicate_audit, "Figure3k_duplicate_audit.tsv")

    def marker_area(percentages: np.ndarray | list[float]) -> np.ndarray:
        values = np.asarray(percentages, dtype=float)
        return 34.0 + np.clip(values - 35.0, 0.0, 60.0) * 3.7

    fig, ax = plt.subplots(figsize=(6.1, 4.0), facecolor="white")
    stage_x = np.arange(len(PHASE_ORDER), dtype=float) * 2.0
    offsets = {"FL": -0.30, "TN": 0.30}
    markers = {"FL": "o", "TN": "s"}
    outlines = {"FL": "#E76F51", "TN": "#3274A1"}
    module_y = {
        module: len(MODULE_ORDER) - 1 - index
        for index, module in enumerate(MODULE_ORDER)
    }
    stage_position = {stage: index for index, stage in enumerate(PHASE_ORDER)}
    colour_limit = max(
        1.0,
        float(
            np.nanpercentile(
                np.abs(summary.median_of_within_gene_robust_log2_ratios.to_numpy(float)),
                98,
            )
        ),
    )
    norm = TwoSlopeNorm(vmin=-colour_limit, vcenter=0.0, vmax=colour_limit)
    cmap = LinearSegmentedColormap.from_list(
        "unique_gene_ASE", ["#2F6FA3", "#F7F7F3", "#C9554D"]
    )
    for row in range(len(MODULE_ORDER)):
        if row % 2:
            ax.axhspan(row - 0.48, row + 0.48, color="#F4F6F6", lw=0, zorder=0)
    for index, center in enumerate(stage_x):
        face = "#FCF4E2" if index == len(PHASE_ORDER) - 1 else "#F7F9F9"
        ax.axvspan(center - 0.78, center + 0.78, color=face, lw=0, zorder=0)
        if index < len(PHASE_ORDER) - 1:
            ax.axvline(center + 1.0, color="#D9E0E3", lw=0.7, zorder=1)
    for analysis in ["FL", "TN"]:
        part = summary[summary.analysis.eq(analysis)]
        xs = [
            stage_x[stage_position[str(stage)]] + offsets[analysis]
            for stage in part.stage_group
        ]
        ys = [module_y[str(module)] for module in part.trait_module]
        ax.scatter(
            xs,
            ys,
            s=marker_area(part.robust_ASE_unique_gene_percentage.to_numpy(float)),
            c=part.median_of_within_gene_robust_log2_ratios.to_numpy(float),
            cmap=cmap,
            norm=norm,
            marker=markers[analysis],
            edgecolor=outlines[analysis],
            linewidth=1.25,
            zorder=4,
        )
    ax.set_xlim(stage_x[0] - 0.82, stage_x[-1] + 0.82)
    ax.set_ylim(-0.72, len(MODULE_ORDER) - 0.18)
    ax.set_xticks(
        stage_x,
        ["Early\n0-65 d", "Middle\n80-140 d", "Late\n155-185 d", "Postharvest\n12-72 h"],
    )
    ax.set_yticks(
        [module_y[module] for module in MODULE_ORDER],
        [MODULE_SHORT[module] for module in MODULE_ORDER],
    )
    ax.tick_params(axis="x", length=0, pad=9, colors="#3D474B")
    ax.tick_params(axis="y", length=0, pad=10, colors="#273238")
    for spine in ax.spines.values():
        spine.set_visible(False)
    for center in stage_x:
        for analysis in ["FL", "TN"]:
            ax.text(
                center + offsets[analysis],
                len(MODULE_ORDER) - 0.39,
                analysis,
                ha="center",
                va="bottom",
                color=outlines[analysis],
                fontsize=7.1,
                fontweight="bold",
            )
    genotype_handles = [
        Line2D(
            [0],
            [0],
            marker=markers[analysis],
            color="none",
            markerfacecolor="white",
            markeredgecolor=outlines[analysis],
            markeredgewidth=1.3,
            markersize=7,
            label=analysis,
        )
        for analysis in ["FL", "TN"]
    ]
    size_handles = [
        ax.scatter(
            [],
            [],
            s=marker_area([percentage])[0],
            color="#D4D8D9",
            edgecolor="#747D81",
            linewidth=0.7,
            label=f"{percentage}%",
        )
        for percentage in [50, 70, 90]
    ]
    genotype_legend = ax.legend(
        handles=genotype_handles,
        loc="lower left",
        bbox_to_anchor=(-0.005, -0.43),
        ncol=2,
        frameon=False,
        handletextpad=0.45,
        columnspacing=1.1,
    )
    ax.add_artist(genotype_legend)
    ax.legend(
        handles=size_handles,
        title="Robust ASE (unique genes)",
        loc="lower left",
        bbox_to_anchor=(0.15, -0.45),
        ncol=3,
        frameon=False,
        handletextpad=0.35,
        columnspacing=0.85,
        fontsize=6.7,
        title_fontsize=6.7,
    )
    color_axis = fig.add_axes([0.665, 0.070, 0.245, 0.027])
    colorbar = fig.colorbar(
        mpl.cm.ScalarMappable(norm=norm, cmap=cmap),
        cax=color_axis,
        orientation="horizontal",
    )
    colorbar.set_ticks([-colour_limit, 0, colour_limit])
    colorbar.set_ticklabels(
        [f"-{colour_limit:.1f}", "0", f"+{colour_limit:.1f}"]
    )
    colorbar.outline.set_visible(False)
    colorbar.ax.tick_params(labelsize=6.6, length=2, pad=1)
    colorbar.set_label(
        "Gene-level median log2(A/B)",
        fontsize=6.7,
        labelpad=1,
    )
    panel_label(ax, "k", -0.06, 1.05)
    fig.subplots_adjust(left=0.11, right=0.98, top=0.94, bottom=0.31)
    save(fig, "Figure3k_revised_unique_gene_trait_ASE")

    checks = [
        {
            "panel": "3k",
            "check": "exact_duplicate_module_gene_stage_rows_removed",
            "observed": duplicate_rows,
            "expected": 265,
            "pass": duplicate_rows == 265,
        },
        {
            "panel": "3k",
            "check": "duplicate_key_value_conflicts",
            "observed": conflict_keys,
            "expected": 0,
            "pass": conflict_keys == 0,
        },
        {
            "panel": "3k",
            "check": "all_2x6x4_combinations_present",
            "observed": len(summary),
            "expected": 48,
            "pass": len(summary) == 48,
        },
        {
            "panel": "3k",
            "check": "robust_unique_genes_do_not_exceed_eligible_unique_genes",
            "observed": int(
                (summary.robust_ASE_unique_genes <= summary.eligible_unique_genes).sum()
            ),
            "expected": len(summary),
            "pass": bool(
                (
                    summary.robust_ASE_unique_genes
                    <= summary.eligible_unique_genes
                ).all()
            ),
        },
        {
            "panel": "3k",
            "check": "unique_gene_percentages_within_0_100",
            "observed": float(summary.robust_ASE_unique_gene_percentage.max()),
            "expected": "<=100",
            "pass": bool(
                summary.robust_ASE_unique_gene_percentage.between(0, 100).all()
            ),
        },
    ]
    return summary, checks


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_manifest() -> None:
    rows = []
    paths = sorted(OUT.glob("Figure3*")) + [Path(__file__).resolve()]
    for path in paths:
        if not path.is_file() or path.name == "Figure3_output_sha256.tsv":
            continue
        rows.append(
            {
                "file": str(path.relative_to(OUT)),
                "bytes": path.stat().st_size,
                "sha256": file_sha256(path),
            }
        )
    write_tsv(pd.DataFrame(rows), "Figure3_output_sha256.tsv")


def write_input_manifest() -> None:
    inputs = [
        CLASS_SOURCE,
        MIRROR_SOURCE,
        FLOW_SOURCE,
        TEST_REFERENCE,
        MODE_LEGACY_SUMMARY,
        MODE_DETAIL_SOURCE,
        REGULATORY_SOURCE,
        TRAIT_SOURCE,
        LEGACY_INHERITANCE_SOURCE,
        TRAIT_LEGACY_SUMMARY,
        *SNP_PATHS.values(),
    ]
    rows = []
    for path in inputs:
        rows.append(
            {
                "file": str(path),
                "bytes": path.stat().st_size,
                "sha256": file_sha256(path),
            }
        )
    write_tsv(pd.DataFrame(rows), "Figure3_input_sha256.tsv")


def main() -> None:
    set_style()
    all_checks: list[dict[str, object]] = []
    print("[1/7] Figure 3d formal full-universe composition", flush=True)
    _, checks = redraw_3d_formal()
    all_checks.extend(checks)
    print("[2/7] Figure 3c 19-stage directional ASE", flush=True)
    _, checks = redraw_3c()
    all_checks.extend(checks)
    print("[3/7] Figure 3d alternative flow (7,059 shared orthogroups)", flush=True)
    _, checks = redraw_3d_intersection_sensitivity()
    all_checks.extend(checks)
    print("[4/7] Figure 3e FL formal SNP-density tests", flush=True)
    _, _, checks = redraw_3e()
    all_checks.extend(checks)
    print("[5/7] Figure 3g deduplicated modes with Conserved", flush=True)
    _, checks = redraw_3g()
    all_checks.extend(checks)
    print("[6/7] Figure 3j gene-clustered association", flush=True)
    _, checks = redraw_3j()
    all_checks.extend(checks)
    print("[7/7] Figure 3k unique-gene trait ASE", flush=True)
    _, checks = redraw_3k()
    all_checks.extend(checks)
    validation = pd.DataFrame(all_checks)
    write_tsv(validation, "Figure3_redraw_validation.tsv")
    if not validation["pass"].all():
        failed = validation.loc[~validation["pass"], ["panel", "check"]]
        raise RuntimeError(f"Validation failed:\n{failed.to_string(index=False)}")
    write_input_manifest()
    write_manifest()
    run_summary = pd.DataFrame(
        [
            {
                "status": "PASS",
                "validation_checks": len(validation),
                "validation_failures": int((~validation["pass"]).sum()),
                "figure3b_status": "RETAINED_LEGACY_NOT_INDEPENDENTLY_REVALIDATED",
                "figure3b_recalculation": "NOT_RUN_AUTHOR_DECISION",
                "note": "Figure 3b retained by author decision; original manuscript figures were not overwritten",
            }
        ]
    )
    write_tsv(run_summary, "Figure3_redraw_run_summary.tsv")
    write_manifest()
    print(f"PASS: {len(validation)} validation checks", flush=True)


if __name__ == "__main__":
    main()
