#!/usr/bin/env python3
"""Rebuild Figure 3e and 3f with all four current ASE classes."""

from __future__ import annotations

import csv
import hashlib
import itertools
import json
import os
import platform
import sys
import tempfile
from datetime import datetime, timezone
from pathlib import Path

import matplotlib as mpl

mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import numpy as np
import pandas as pd
import scipy
from scipy import stats
from PIL import Image


HERE = Path(__file__).resolve().parent
ANALYSIS = Path("${ANALYSIS_DIR}")
CURRENT = (
    ANALYSIS
    / "22_answer_reviews/00_ms/03_V3/02_figure/_runs/"
    "RUN-ASE-FIGURES-CURRENT-V3-001/output"
)
CURRENT_WORK = CURRENT.parent / "work/tn_kaks"
PANEL_SOURCE = (
    ANALYSIS
    / "22_answer_reviews/00_ms/03_V3/03_figure3/"
    "04_FL_TN_ASE_reference_panels"
)

CLASS_FILE = CURRENT / "gene_overall_class_current.tsv"
TN_SNP_FILE = CURRENT / "TN_diagnostic_SNP_density_current.tsv.gz"
FL_KAKS_FILE = (
    ANALYSIS
    / "22_answer_reviews/00_ms/03_V3/03_figure3/03_3c/"
    "final.Fig3c_American_hap1_Africa_hap2_KaKs_all.tsv"
)
TN_KAKS_FILE = CURRENT_WORK / "TN_Dura_Pisifera_KaKs_all_current.tsv"
THREE_CLASS_CONTROL = PANEL_SOURCE / "source_F_KaKs_by_ASE_class.tsv"
FIG3E_RATIO_REFERENCE = (
    ANALYSIS
    / "22_answer_reviews/00_ms/05_MS/MS_revision_3/Final_20260823_VectorRevision/"
    "8月29日/FINAL_ALL_FIGURES_ONE_FOLDER_20260901/"
    "Figure3e_RECOMMENDED_old_style_TN_SNP_density_KW.pdf"
)
FIG3F_RATIO_REFERENCE = PANEL_SOURCE / "01F_KaKs.pdf"

# The order and colors are locked to the current Figure 3d legend.
CLASSES = ["NoDiff", "HapDom", "Sub", "NoASE"]
COLORS = {
    "NoDiff": "#F8766D",
    "HapDom": "#FFF59D",
    "Sub": "#8DD3C7",
    "NoASE": "#4EA5DF",
}
REGIONS = [
    ("upstream2kb_per_kb", "Upstream\n2 kb"),
    ("gene_body_per_kb", "Gene"),
    ("downstream2kb_per_kb", "Downstream\n2 kb"),
]

EXPECTED_CLASS_UNIVERSE = {
    ("FL", "HapDom"): 9413,
    ("FL", "NoDiff"): 3863,
    ("FL", "Sub"): 3118,
    ("FL", "NoASE"): 1222,
    ("TN", "HapDom"): 6374,
    ("TN", "NoDiff"): 1502,
    ("TN", "Sub"): 3308,
    ("TN", "NoASE"): 225,
}
EXPECTED_THREE_CLASS_KAKS = {
    ("FL", "HapDom"): 6411,
    ("FL", "NoDiff"): 2587,
    ("FL", "Sub"): 2242,
    ("TN", "HapDom"): 2423,
    ("TN", "NoDiff"): 541,
    ("TN", "Sub"): 1523,
}

FIG3E_STEM = "Figure3e_four_ASE_classes_TN_SNP_density"
FIG3F_STEM = "Figure3f_four_ASE_classes_KaKs_Ks"
COMMAND = f"${DATA_DIR}/miniconda3/bin/python {Path(__file__).resolve()}"

# Fixed canvases match the small-panel slots used in the assembled Figure 3.
FIG3E_CANVAS_IN = (3.80, 3.52)
FIG3F_REFERENCE_CANVAS_PT = (255.945, 219.856)
FIG3F_CANVAS_IN = tuple(value / 72.0 for value in FIG3F_REFERENCE_CANVAS_PT)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def atomic_write_text(path: Path, text: str) -> None:
    staging = path.with_name(f".{path.name}.staging")
    staging.write_text(text, encoding="ascii")
    os.replace(staging, path)


def atomic_write_dataframe(frame: pd.DataFrame, path: Path) -> None:
    staging = path.with_name(f".{path.name}.staging")
    frame.to_csv(staging, sep="\t", index=False, na_rep="NA")
    os.replace(staging, path)


def configure_style(base_size: float) -> None:
    mpl.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": [
                "Arial",
                "Helvetica",
                "Liberation Sans",
                "DejaVu Sans",
            ],
            "font.size": base_size,
            "axes.labelsize": base_size,
            "xtick.labelsize": base_size * 0.86,
            "ytick.labelsize": base_size * 0.86,
            "legend.fontsize": base_size * 0.78,
            "axes.linewidth": 0.8,
            "xtick.major.width": 0.8,
            "ytick.major.width": 0.8,
            "xtick.major.size": 3.5,
            "ytick.major.size": 3.5,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
            "savefig.facecolor": "white",
            "figure.facecolor": "white",
        }
    )


def save_figure(fig: plt.Figure, stem: str, title: str, subject: str) -> list[Path]:
    outputs: list[Path] = []
    with tempfile.TemporaryDirectory(prefix=f"{stem}_", dir=HERE) as temp_dir:
        temp = Path(temp_dir)
        for extension in ("pdf", "svg", "png"):
            destination = HERE / f"{stem}.{extension}"
            staged = temp / destination.name
            kwargs: dict[str, object] = {
                "format": extension,
                "facecolor": "white",
            }
            if extension == "png":
                kwargs["dpi"] = 600
            if extension == "pdf":
                kwargs["metadata"] = {
                    "Title": title,
                    "Subject": subject,
                    "Creator": Path(__file__).name,
                }
            fig.savefig(staged, **kwargs)
            os.replace(staged, destination)
            outputs.append(destination)
    plt.close(fig)
    return outputs


def normalise_tu(values: pd.Series) -> pd.Series:
    return values.astype(str).str.replace("evm.model.", "evm.TU.", regex=False)


def add_check(
    checks: list[dict[str, object]],
    name: str,
    observed: object,
    expected: object,
    passed: bool,
) -> None:
    checks.append(
        {
            "check": name,
            "observed": str(observed),
            "expected": str(expected),
            "pass": bool(passed),
        }
    )


def require_checks(checks: list[dict[str, object]], stage: str) -> None:
    failed = [row["check"] for row in checks if not bool(row["pass"])]
    if failed:
        raise RuntimeError(f"{stage} validation failed: " + ", ".join(failed))


def load_classes(checks: list[dict[str, object]]) -> pd.DataFrame:
    frame = pd.read_csv(CLASS_FILE, sep="\t")
    duplicate_n = int(frame.duplicated(["analysis", "gene_id"]).sum())
    add_check(checks, "class_map_duplicate_keys", duplicate_n, 0, duplicate_n == 0)
    counts = frame.groupby(["analysis", "overall_class"]).size().to_dict()
    for key, expected in EXPECTED_CLASS_UNIVERSE.items():
        observed = int(counts.get(key, 0))
        add_check(
            checks,
            f"class_universe:{key[0]}:{key[1]}",
            observed,
            expected,
            observed == expected,
        )
    require_checks(checks, "ASE class map")
    return frame[["analysis", "gene_id", "overall_class"]].copy()


def build_kaks_source(
    classes: pd.DataFrame, checks: list[dict[str, object]]
) -> pd.DataFrame:
    parts: list[pd.DataFrame] = []
    specifications = [
        ("FL", FL_KAKS_FILE, "africa_gene"),
        ("TN", TN_KAKS_FILE, "gene_dura"),
    ]
    for analysis, path, gene_column in specifications:
        raw = pd.read_csv(path, sep="\t", low_memory=False)
        part = pd.DataFrame(
            {
                "analysis": analysis,
                "gene_id": normalise_tu(raw[gene_column]),
                "Ka": pd.to_numeric(raw["Ka"], errors="coerce"),
                "Ks": pd.to_numeric(raw["Ks"], errors="coerce"),
                "Ka_Ks": pd.to_numeric(raw["Ka/Ks"], errors="coerce"),
                "source_pair_id": raw["pair_id"].astype(str),
            }
        )
        class_part = classes.loc[classes.analysis.eq(analysis), ["gene_id", "overall_class"]]
        part = part.merge(class_part, on="gene_id", how="left", validate="one_to_one")
        part = part.rename(columns={"overall_class": "ASE_type"})
        part = part.loc[
            part.ASE_type.isin(CLASSES)
            & part.Ka.notna()
            & part.Ks.gt(0)
            & part.Ks.le(0.10)
            & part.Ka_Ks.ge(0)
            & part.Ka_Ks.le(3)
        ].copy()
        parts.append(part)

    result = pd.concat(parts, ignore_index=True)
    result = result[
        ["analysis", "gene_id", "ASE_type", "Ka", "Ks", "Ka_Ks", "source_pair_id"]
    ].sort_values(["analysis", "ASE_type", "gene_id"])
    duplicate_n = int(result.duplicated(["analysis", "gene_id"]).sum())
    add_check(checks, "KaKs_duplicate_gene_keys", duplicate_n, 0, duplicate_n == 0)

    control = pd.read_csv(THREE_CLASS_CONTROL, sep="\t")
    for key, expected_count in EXPECTED_THREE_CLASS_KAKS.items():
        analysis, category = key
        observed_genes = set(
            result.loc[
                result.analysis.eq(analysis) & result.ASE_type.eq(category), "gene_id"
            ]
        )
        control_genes = set(
            control.loc[
                control.analysis.eq(analysis) & control.ASE_type.eq(category), "gene_id"
            ]
        )
        add_check(
            checks,
            f"KaKs_count_regression:{analysis}:{category}",
            len(observed_genes),
            expected_count,
            len(observed_genes) == expected_count,
        )
        add_check(
            checks,
            f"KaKs_gene_set_regression:{analysis}:{category}",
            len(observed_genes.symmetric_difference(control_genes)),
            0,
            observed_genes == control_genes,
        )

    for analysis in ("FL", "TN"):
        noase_n = int(
            result.loc[
                result.analysis.eq(analysis) & result.ASE_type.eq("NoASE")
            ].shape[0]
        )
        add_check(
            checks,
            f"KaKs_NoASE_nonempty:{analysis}",
            noase_n,
            ">=20",
            noase_n >= 20,
        )
    require_checks(checks, "four-class Ka/Ks reconstruction")
    return result


def holm_adjust(pvalues: list[float]) -> np.ndarray:
    values = np.asarray(pvalues, dtype=float)
    order = np.argsort(values)
    adjusted_ordered = np.maximum.accumulate(
        values[order] * (len(values) - np.arange(len(values)))
    )
    adjusted = np.empty_like(adjusted_ordered)
    adjusted[order] = np.minimum(adjusted_ordered, 1.0)
    return adjusted


def compact_letters(
    significant: dict[tuple[int, int], bool], medians: list[float]
) -> list[str]:
    nodes = tuple(range(len(CLASSES)))
    cliques: list[frozenset[int]] = []
    for size in range(1, len(nodes) + 1):
        for subset in itertools.combinations(nodes, size):
            if all(
                not significant[tuple(sorted(pair))]
                for pair in itertools.combinations(subset, 2)
            ):
                cliques.append(frozenset(subset))
    maximal = [
        clique
        for clique in cliques
        if not any(clique < other for other in cliques)
    ]
    maximal.sort(
        key=lambda clique: (
            -max(medians[index] for index in clique),
            -np.mean([medians[index] for index in clique]),
            tuple(sorted(clique)),
        )
    )
    alphabet = "abcdefghijklmnopqrstuvwxyz"
    if len(maximal) > len(alphabet):
        raise RuntimeError("Too many compact-letter groups")
    labels = ["" for _ in nodes]
    for letter, clique in zip(alphabet, maximal):
        for index in clique:
            labels[index] += letter

    for first, second in itertools.combinations(nodes, 2):
        shared = bool(set(labels[first]).intersection(labels[second]))
        if significant[(first, second)] == shared:
            raise RuntimeError("Invalid compact-letter display")
    return labels


def build_snp_sources(
    checks: list[dict[str, object]],
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, dict[str, dict[str, str]]]:
    raw = pd.read_csv(TN_SNP_FILE, sep="\t", compression="gzip")
    tn = raw.loc[raw.analysis.eq("TN") & raw.overall_class.isin(CLASSES)].copy()
    long = tn.melt(
        id_vars=["gene_id", "analysis", "overall_class"],
        value_vars=[region for region, _ in REGIONS],
        var_name="region",
        value_name="SNPs_per_1000bp",
    )
    long["SNPs_per_1000bp"] = pd.to_numeric(
        long["SNPs_per_1000bp"], errors="coerce"
    )
    invalid = int(
        long.SNPs_per_1000bp.isna().sum() + long.SNPs_per_1000bp.lt(0).sum()
    )
    add_check(checks, "Figure3e_invalid_SNP_density", invalid, 0, invalid == 0)

    for category in CLASSES:
        expected = EXPECTED_CLASS_UNIVERSE[("TN", category)]
        counts = long.loc[long.overall_class.eq(category)].groupby("region").size()
        passed = len(counts) == len(REGIONS) and bool((counts == expected).all())
        add_check(
            checks,
            f"Figure3e_count_each_region:{category}",
            ";".join(str(int(value)) for value in counts.to_numpy()),
            ";".join([str(expected)] * len(REGIONS)),
            passed,
        )

    summary_rows: list[dict[str, object]] = []
    test_rows: list[dict[str, object]] = []
    letters_by_region: dict[str, dict[str, str]] = {}
    for region, _ in REGIONS:
        groups = [
            long.loc[
                long.region.eq(region) & long.overall_class.eq(category),
                "SNPs_per_1000bp",
            ].to_numpy(float)
            for category in CLASSES
        ]
        for category, values in zip(CLASSES, groups):
            q1, median, q3 = np.quantile(values, [0.25, 0.50, 0.75])
            summary_rows.append(
                {
                    "analysis": "TN",
                    "region": region,
                    "overall_class": category,
                    "n": len(values),
                    "q1": q1,
                    "median": median,
                    "q3": q3,
                    "mean": float(np.mean(values)),
                }
            )

        kw = stats.kruskal(*groups)
        pairs = list(itertools.combinations(range(len(CLASSES)), 2))
        raw_pvalues = [
            stats.mannwhitneyu(groups[first], groups[second], alternative="two-sided").pvalue
            for first, second in pairs
        ]
        adjusted = holm_adjust(raw_pvalues)
        significant = {
            pair: bool(value < 0.05) for pair, value in zip(pairs, adjusted)
        }
        medians = [float(np.median(values)) for values in groups]
        letters = compact_letters(significant, medians)
        letters_by_region[region] = dict(zip(CLASSES, letters))
        letter_text = ";".join(
            f"{category}:{letter}" for category, letter in zip(CLASSES, letters)
        )
        test_rows.append(
            {
                "analysis": "TN",
                "region": region,
                "test": "Kruskal-Wallis",
                "comparison": "all four ASE classes",
                "group_1": "",
                "group_2": "",
                "n_1": "",
                "n_2": "",
                "pvalue": float(kw.pvalue),
                "holm_pvalue": np.nan,
                "significant_0.05": bool(kw.pvalue < 0.05),
                "letters": letter_text,
            }
        )
        for (first, second), raw_p, holm_p in zip(pairs, raw_pvalues, adjusted):
            test_rows.append(
                {
                    "analysis": "TN",
                    "region": region,
                    "test": "Mann-Whitney U",
                    "comparison": f"{CLASSES[first]} vs {CLASSES[second]}",
                    "group_1": CLASSES[first],
                    "group_2": CLASSES[second],
                    "n_1": len(groups[first]),
                    "n_2": len(groups[second]),
                    "pvalue": float(raw_p),
                    "holm_pvalue": float(holm_p),
                    "significant_0.05": bool(holm_p < 0.05),
                    "letters": letter_text,
                }
            )

    summary = pd.DataFrame(summary_rows)
    tests = pd.DataFrame(test_rows)
    add_check(
        checks,
        "Figure3e_test_row_count",
        len(tests),
        21,
        len(tests) == 21,
    )
    require_checks(checks, "four-class TN SNP-density data")
    return long, summary, tests, letters_by_region


def format_p(value: float) -> str:
    if value < 0.001:
        exponent = int(np.floor(np.log10(value)))
        coefficient = value / (10**exponent)
        return rf"$P = {coefficient:.2f} \times 10^{{{exponent}}}$"
    return f"P = {value:.3f}"


def render_figure3e(
    plotdata: pd.DataFrame,
    tests: pd.DataFrame,
    letters_by_region: dict[str, dict[str, str]],
) -> list[Path]:
    configure_style(7.8)
    fig, ax = plt.subplots(figsize=FIG3E_CANVAS_IN, facecolor="white")
    positions = [
        region_index * 5 + class_index
        for region_index in range(len(REGIONS))
        for class_index in range(len(CLASSES))
    ]
    arrays = [
        plotdata.loc[
            plotdata.region.eq(region) & plotdata.overall_class.eq(category),
            "SNPs_per_1000bp",
        ].to_numpy(float)
        for region, _ in REGIONS
        for category in CLASSES
    ]
    box = ax.boxplot(
        arrays,
        positions=positions,
        widths=0.70,
        showfliers=False,
        patch_artist=True,
        medianprops={"lw": 0.8},
        whiskerprops={"lw": 0.65},
        capprops={"lw": 0.65},
    )
    color_sequence = [COLORS[category] for _ in REGIONS for category in CLASSES]
    for index, color in enumerate(color_sequence):
        box["boxes"][index].set_facecolor(mpl.colors.to_rgba(color, 0.18))
        box["boxes"][index].set_edgecolor(color)
        box["boxes"][index].set_linewidth(0.9)
        box["medians"][index].set_color(color)
        for artist in box["whiskers"][2 * index : 2 * index + 2]:
            artist.set_color(color)
        for artist in box["caps"][2 * index : 2 * index + 2]:
            artist.set_color(color)

    whisker_max = max(float(np.max(line.get_ydata())) for line in box["whiskers"])
    letter_y = whisker_max * 1.055
    p_y = whisker_max * 1.145
    ylim_top = whisker_max * 1.27
    kw = tests.loc[tests.test.eq("Kruskal-Wallis")].set_index("region")
    for region_index, (region, _) in enumerate(REGIONS):
        center = region_index * 5 + 1.5
        ax.text(
            center,
            p_y,
            format_p(float(kw.loc[region, "pvalue"])),
            ha="center",
            va="bottom",
            fontsize=5.9,
        )
        for class_index, category in enumerate(CLASSES):
            ax.text(
                region_index * 5 + class_index,
                letter_y,
                letters_by_region[region][category],
                ha="center",
                va="bottom",
                fontsize=6.8,
                fontweight="bold",
            )

    ax.set_xlim(-0.75, 13.75)
    ax.set_ylim(-0.025 * ylim_top, ylim_top)
    ax.set_xticks([1.5, 6.5, 11.5], [label for _, label in REGIONS])
    ax.set_ylabel("SNPs per 1,000 bp")
    ax.set_title("TN", loc="left", fontsize=8.2, fontweight="bold", pad=8)
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(direction="out")
    handles = [
        Patch(
            facecolor=mpl.colors.to_rgba(COLORS[category], 0.18),
            edgecolor=COLORS[category],
            linewidth=0.9,
            label=category,
        )
        for category in CLASSES
    ]
    fig.legend(
        handles=handles,
        frameon=False,
        ncol=4,
        loc="upper center",
        bbox_to_anchor=(0.60, 0.995),
        columnspacing=0.7,
        handlelength=1.0,
    )
    fig.text(0.012, 0.992, "e", fontsize=12.5, fontweight="bold", ha="left", va="top")
    fig.subplots_adjust(left=0.18, right=0.985, top=0.79, bottom=0.16)
    return save_figure(
        fig,
        FIG3E_STEM,
        "Figure 3e: TN SNP density across four ASE classes",
        "Raw diagnostic SNP-density boxplots with four current ASE classes",
    )


def kde_curve(values: pd.Series, xmax: float) -> tuple[np.ndarray, np.ndarray]:
    vector = pd.to_numeric(values, errors="coerce").to_numpy(float)
    vector = vector[np.isfinite(vector) & (vector >= 0) & (vector <= xmax)]
    if len(vector) < 20 or np.unique(vector).size < 3:
        raise RuntimeError("Insufficient values for KDE")
    x = np.linspace(0, xmax, 420)
    return x, stats.gaussian_kde(vector)(x)


def render_figure3f(kaks: pd.DataFrame) -> list[Path]:
    configure_style(7.8)
    fig = plt.figure(figsize=FIG3F_CANVAS_IN, facecolor="white")
    ax = fig.add_axes([0.15, 0.15, 0.82, 0.79])
    inset = ax.inset_axes([0.54, 0.48, 0.43, 0.46])
    line_styles = {"FL": "-", "TN": (0, (8, 4))}

    for analysis in ("FL", "TN"):
        for category in CLASSES:
            group = kaks.loc[
                kaks.analysis.eq(analysis) & kaks.ASE_type.eq(category)
            ]
            x, y = kde_curve(group.Ka_Ks, 3.0)
            ax.plot(
                x,
                y,
                color=COLORS[category],
                lw=1.35,
                ls=line_styles[analysis],
                dash_capstyle="butt",
                zorder=3 if analysis == "TN" else 2,
            )
            x_ks, y_ks = kde_curve(group.Ks, 0.10)
            inset.plot(
                x_ks,
                y_ks,
                color=COLORS[category],
                lw=1.0,
                ls=line_styles[analysis],
                dash_capstyle="butt",
                zorder=3 if analysis == "TN" else 2,
            )

    ax.axvline(1, color="#777777", lw=0.75, ls=(0, (3, 3)), zorder=1)
    ax.set_xlim(0, 3)
    ax.set_ylim(bottom=0)
    ax.set_xticks([0, 1, 2, 3])
    ax.set_xlabel(r"$K_a/K_s$ ratio")
    ax.set_ylabel("Density")
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(direction="out")

    inset.set_xlim(0, 0.10)
    inset.set_ylim(bottom=0)
    inset.set_xticks([0, 0.025, 0.050, 0.075, 0.100])
    inset.set_xlabel(r"$K_s$ value", fontsize=6.3)
    inset.set_ylabel("Density", fontsize=6.3)
    inset.tick_params(labelsize=5.6, direction="out", length=2.4)
    for spine in inset.spines.values():
        spine.set_color("#888888")
        spine.set_linewidth(0.65)

    class_handles = [
        Line2D([0], [0], color=COLORS[category], lw=1.8, label=category)
        for category in CLASSES
    ]
    material_handles = [
        Line2D([0], [0], color="#333333", lw=1.6, ls="-", label="FL"),
        Line2D(
            [0],
            [0],
            color="#333333",
            lw=1.6,
            ls=(0, (8, 4)),
            label="TN",
        ),
    ]
    legend_classes = ax.legend(
        handles=class_handles,
        title="ASE class",
        frameon=False,
        ncol=1,
        loc="lower left",
        bbox_to_anchor=(0.57, 0.035),
        borderaxespad=0,
        fontsize=5.6,
        title_fontsize=6.0,
        labelspacing=0.25,
        handlelength=1.7,
    )
    ax.add_artist(legend_classes)
    ax.legend(
        handles=material_handles,
        title="Material",
        frameon=False,
        ncol=1,
        loc="lower left",
        bbox_to_anchor=(0.82, 0.105),
        borderaxespad=0,
        fontsize=5.6,
        title_fontsize=6.0,
        labelspacing=0.35,
        handlelength=3.0,
    )
    fig.text(0.015, 0.985, "f", fontsize=12.5, fontweight="bold", ha="left", va="top")
    return save_figure(
        fig,
        FIG3F_STEM,
        "Figure 3f: Ka/Ks distributions across four ASE classes",
        "FL solid and TN long-dashed density curves; inset shows Ks",
    )


def summarize_kaks(kaks: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for (analysis, category), group in kaks.groupby(["analysis", "ASE_type"]):
        ka_q = np.quantile(group.Ka_Ks, [0.25, 0.50, 0.75])
        ks_q = np.quantile(group.Ks, [0.25, 0.50, 0.75])
        rows.append(
            {
                "analysis": analysis,
                "ASE_type": category,
                "n": len(group),
                "Ka_Ks_q1": ka_q[0],
                "Ka_Ks_median": ka_q[1],
                "Ka_Ks_q3": ka_q[2],
                "Ks_q1": ks_q[0],
                "Ks_median": ks_q[1],
                "Ks_q3": ks_q[2],
            }
        )
    return pd.DataFrame(rows).sort_values(["analysis", "ASE_type"])


def write_input_manifest() -> Path:
    roles = {
        CLASS_FILE: "current overall ASE classes",
        TN_SNP_FILE: "current TN diagnostic SNP densities",
        FL_KAKS_FILE: "current FL allele-pair Ka/Ks",
        TN_KAKS_FILE: "current TN allele-pair Ka/Ks",
        THREE_CLASS_CONTROL: "three-class Figure 3f regression control",
        FIG3E_RATIO_REFERENCE: "assembled Figure 3e small-panel ratio reference",
        FIG3F_RATIO_REFERENCE: "assembled Figure 3f small-panel ratio reference",
    }
    rows = [
        {
            "role": role,
            "path": str(path),
            "bytes": path.stat().st_size,
            "sha256": sha256(path),
        }
        for path, role in roles.items()
    ]
    path = HERE / "input_manifest.tsv"
    atomic_write_dataframe(pd.DataFrame(rows), path)
    return path


def main() -> None:
    started = datetime.now(timezone.utc)
    checks: list[dict[str, object]] = []
    input_manifest = write_input_manifest()
    classes = load_classes(checks)
    kaks = build_kaks_source(classes, checks)
    snp_long, snp_summary, snp_tests, letters = build_snp_sources(checks)

    source_paths = [
        HERE / "Figure3e_four_ASE_classes_plotdata.tsv",
        HERE / "Figure3e_four_ASE_classes_summary.tsv",
        HERE / "Figure3e_four_ASE_classes_tests.tsv",
        HERE / "Figure3f_four_ASE_classes_plotdata.tsv",
        HERE / "Figure3f_four_ASE_classes_summary.tsv",
    ]
    atomic_write_dataframe(snp_long, source_paths[0])
    atomic_write_dataframe(snp_summary, source_paths[1])
    atomic_write_dataframe(snp_tests, source_paths[2])
    atomic_write_dataframe(kaks, source_paths[3])
    atomic_write_dataframe(summarize_kaks(kaks), source_paths[4])

    figure_paths = render_figure3e(snp_long, snp_tests, letters)
    figure_paths.extend(render_figure3f(kaks))

    f_svg = HERE / f"{FIG3F_STEM}.svg"
    dash_present = "stroke-dasharray" in f_svg.read_text(encoding="utf-8")
    add_check(
        checks,
        "Figure3f_SVG_contains_dashed_strokes",
        dash_present,
        True,
        dash_present,
    )
    expected_png_dimensions = {
        f"{FIG3E_STEM}.png": tuple(
            int(value * 600) for value in FIG3E_CANVAS_IN
        ),
        f"{FIG3F_STEM}.png": tuple(
            int(value * 600) for value in FIG3F_CANVAS_IN
        ),
    }
    for filename, expected in expected_png_dimensions.items():
        with Image.open(HERE / filename) as rendered:
            observed = rendered.size
        add_check(
            checks,
            f"fixed_small_panel_canvas:{filename}",
            f"{observed[0]}x{observed[1]}",
            f"{expected[0]}x{expected[1]}",
            observed == expected,
        )
    for path in source_paths + figure_paths:
        add_check(
            checks,
            f"nonempty_output:{path.name}",
            path.stat().st_size,
            ">100 bytes",
            path.stat().st_size > 100,
        )
    require_checks(checks, "output")

    validation_path = HERE / "validation.tsv"
    atomic_write_dataframe(pd.DataFrame(checks), validation_path)

    ended = datetime.now(timezone.utc)
    metadata = {
        "objective": "Rebuild Figure 3e and 3f with four current ASE classes",
        "command": COMMAND,
        "host": platform.node(),
        "python": sys.version,
        "platform": platform.platform(),
        "packages": {
            "matplotlib": mpl.__version__,
            "numpy": np.__version__,
            "pandas": pd.__version__,
            "scipy": scipy.__version__,
        },
        "started_utc": started.isoformat(),
        "ended_utc": ended.isoformat(),
        "elapsed_seconds": (ended - started).total_seconds(),
        "filters": {"Ka_Ks": "0 <= Ka/Ks <= 3", "Ks": "0 < Ks <= 0.10"},
        "small_panel_canvases": {
            "Figure3e_inches": FIG3E_CANVAS_IN,
            "Figure3f_points": FIG3F_REFERENCE_CANVAS_PT,
            "Figure3e_reference": str(FIG3E_RATIO_REFERENCE),
            "Figure3f_reference": str(FIG3F_RATIO_REFERENCE),
        },
        "figure3e_statistics": (
            "Kruskal-Wallis plus six two-sided Mann-Whitney U tests per region, "
            "Holm adjusted"
        ),
    }
    metadata_path = HERE / "execution_metadata.json"
    atomic_write_text(metadata_path, json.dumps(metadata, indent=2, sort_keys=True) + "\n")

    log_lines = [
        "Figure 3e/3f four-class revision completed.",
        f"Started UTC: {started.isoformat()}",
        f"Ended UTC: {ended.isoformat()}",
        f"Command: {COMMAND}",
        "Ka/Ks plotted counts:",
    ]
    for row in summarize_kaks(kaks).itertuples(index=False):
        log_lines.append(f"  {row.analysis} {row.ASE_type}: n={row.n}")
    log_lines.append("Figure 3e Kruskal-Wallis P values:")
    for row in snp_tests.loc[snp_tests.test.eq("Kruskal-Wallis")].itertuples(index=False):
        log_lines.append(f"  {row.region}: P={row.pvalue:.16g}; {row.letters}")
    log_path = HERE / "run.log"
    atomic_write_text(log_path, "\n".join(log_lines) + "\n")

    manifest_targets = (
        source_paths
        + figure_paths
        + [
            input_manifest,
            validation_path,
            metadata_path,
            log_path,
            HERE / "attempt_history.tsv",
        ]
    )
    output_manifest = pd.DataFrame(
        [
            {
                "path": str(path),
                "bytes": path.stat().st_size,
                "sha256": sha256(path),
            }
            for path in manifest_targets
        ]
    )
    output_manifest_path = HERE / "output_manifest.tsv"
    atomic_write_dataframe(output_manifest, output_manifest_path)

    checksum_targets = manifest_targets + [output_manifest_path, Path(__file__).resolve()]
    checksum_lines = [f"{sha256(path)}  {path.name}" for path in checksum_targets]
    atomic_write_text(HERE / "checksums.sha256", "\n".join(checksum_lines) + "\n")

    print(f"PASS: {len(checks)} validation checks")
    for path in figure_paths:
        print(path)


if __name__ == "__main__":
    main()
