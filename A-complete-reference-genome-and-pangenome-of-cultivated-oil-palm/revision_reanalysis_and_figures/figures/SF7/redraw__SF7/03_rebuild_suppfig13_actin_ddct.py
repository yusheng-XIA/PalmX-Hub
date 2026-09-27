#!/usr/bin/env python3
from __future__ import annotations

import csv
import hashlib
import math
import statistics
from collections import Counter
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from openpyxl import load_workbook
from scipy.stats import wilcoxon


RUN = Path(
    "${ANALYSIS_DIR}/22_answer_reviews/"
    "00_ms/05_MS/MS_revision_3/Final_20260823_VectorRevision/8月29日/"
    "Revised_Panels_20260831"
)
WORKBOOK = RUN / "inputs/qPCR_zip3/上交补充实验/数据整理.xlsx"
CSV_DIR = RUN / "inputs/qPCR_zip3/上交补充实验"
OUT = RUN / "SupplementaryFig13_final_consistent"

MATERIALS = ("TK", "TN", "FL", "NS")
PLOT_MATERIALS = ("FL", "TN", "NS", "TK")
COLORS = {
    "FL": "#159A9C",
    "TN": "#E76F51",
    "NS": "#8E6C9E",
    "TK": "#6FJOINT2",
}
PRIMER_NAMES = {1: "Primer 1", 2: "Primer 2", 3: "Primer 3"}

plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["ps.fonttype"] = 42
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Liberation Sans", "Arial", "DejaVu Sans"]
plt.rcParams["font.size"] = 8
plt.rcParams["axes.labelsize"] = 8
plt.rcParams["xtick.labelsize"] = 7
plt.rcParams["ytick.labelsize"] = 7
plt.rcParams["legend.fontsize"] = 7
plt.rcParams["axes.linewidth"] = 0.7

# Match the accepted Illustrator panel rather than introducing a new layout.
LEGACY_PAGE_SIZE_IN = (488.046 / 72.0, 462.035 / 72.0)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def cell_list(values: list[float]) -> str:
    return ",".join(f"{value:.10g}" for value in values)


def index_list(values: list[int]) -> str:
    return ",".join(str(value) for value in values)


def bh_adjust(p_values: list[float]) -> list[float]:
    n = len(p_values)
    order = sorted(range(n), key=lambda index: p_values[index])
    adjusted = [1.0] * n
    running = 1.0
    for rank_index in range(n - 1, -1, -1):
        original_index = order[rank_index]
        rank = rank_index + 1
        running = min(running, p_values[original_index] * n / rank)
        adjusted[original_index] = running
    return adjusted


def stage_blocks(sheet) -> list[tuple[int, str, int, int]]:
    markers: list[tuple[str, int]] = []
    for row in range(1, sheet.max_row + 1):
        value = sheet.cell(row, 17).value  # Q
        if value not in (None, *MATERIALS):
            markers.append((str(value), row))
    blocks = []
    for index, (stage, start) in enumerate(markers, start=1):
        end = markers[index][1] - 1 if index < len(markers) else sheet.max_row
        blocks.append((index, stage, start, end))
    return blocks


def extract_ddct() -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    sheet = load_workbook(WORKBOOK, data_only=True).active
    long_rows: list[dict[str, object]] = []

    for stage_index, stage, start, end in stage_blocks(sheet):
        material_rows = {material: [] for material in MATERIALS}
        for row in range(start, end + 1):
            material = sheet.cell(row, 1).value  # A
            if material in MATERIALS:
                material_rows[str(material)].append(row)
        for material, rows in material_rows.items():
            if len(rows) != 2:
                raise ValueError(f"{stage} {material}: expected two raw-grid rows, found {rows}")

        raw: dict[tuple[str, int], dict[str, object]] = {}
        for material, (first_row, second_row) in material_rows.items():
            actin_values = [sheet.cell(first_row, col).value for col in range(2, 6)]
            target_specs = {
                1: (second_row, range(2, 6)),
                2: (first_row, range(6, 10)),
                3: (second_row, range(6, 10)),
            }
            for primer_set, (target_row, target_cols) in target_specs.items():
                target_values = [sheet.cell(target_row, col).value for col in target_cols]
                actin_indices = [
                    index
                    for index, actin in enumerate(actin_values, start=1)
                    if isinstance(actin, (int, float))
                ]
                target_indices = [
                    index
                    for index, target in enumerate(target_values, start=1)
                    if isinstance(target, (int, float))
                ]
                if len(actin_indices) < 2 or len(target_indices) < 2:
                    raise ValueError(
                        f"{stage} {material} primer {primer_set}: "
                        f"only n_actin={len(actin_indices)}, n_target={len(target_indices)} numeric wells"
                    )
                actin_numeric = [float(actin_values[index - 1]) for index in actin_indices]
                target_numeric = [float(target_values[index - 1]) for index in target_indices]
                mean_actin_ct = statistics.mean(actin_numeric)
                mean_target_ct = statistics.mean(target_numeric)
                delta_ct = mean_target_ct - mean_actin_ct
                raw[(material, primer_set)] = {
                    "stage_index": stage_index,
                    "stage": stage,
                    "phase": "Development" if stage_index <= 13 else "Post-harvest",
                    "material": material,
                    "primer_set": primer_set,
                    "actin_row": first_row,
                    "target_row": target_row,
                    "actin_indices": actin_indices,
                    "target_indices": target_indices,
                    "n_actin": len(actin_indices),
                    "n_target": len(target_indices),
                    "actin_ct": actin_numeric,
                    "target_ct": target_numeric,
                    "mean_actin_ct": mean_actin_ct,
                    "mean_target_ct": mean_target_ct,
                    "delta_ct": delta_ct,
                }

        for primer_set in (1, 2, 3):
            for material in MATERIALS:
                record = raw[(material, primer_set)]
                delta_ct = float(record["delta_ct"])
                log2_target_over_actin = -delta_ct
                long_rows.append(
                    {
                        "stage_index": record["stage_index"],
                        "stage": record["stage"],
                        "phase": record["phase"],
                        "material": material,
                        "primer_set": primer_set,
                        "primer_label": PRIMER_NAMES[primer_set],
                        "actin_source_row": record["actin_row"],
                        "target_source_row": record["target_row"],
                        "actin_finite_indices_1to4": index_list(record["actin_indices"]),
                        "target_finite_indices_1to4": index_list(record["target_indices"]),
                        "n_actin": record["n_actin"],
                        "n_target": record["n_target"],
                        "actin_ct_values": cell_list(record["actin_ct"]),
                        "target_ct_values": cell_list(record["target_ct"]),
                        "mean_actin_ct": record["mean_actin_ct"],
                        "mean_target_ct": record["mean_target_ct"],
                        "delta_ct_mean_target_minus_mean_actin": delta_ct,
                        "log2_target_over_actin_minus_delta_ct": log2_target_over_actin,
                        "target_over_actin_2pow_minus_delta_ct": 2.0**log2_target_over_actin,
                        "calculation": "deltaCt=meanCt(target)-meanCt(Actin); log2(target/Actin)=-deltaCt",
                        "provenance_status": "workbook_raw_grid; plate_map_not_independently_available",
                    }
                )

    by_key = {
        (str(row["stage"]), str(row["material"]), int(row["primer_set"])): row
        for row in long_rows
    }
    summary_rows: list[dict[str, object]] = []
    for stage_index, stage, _, _ in stage_blocks(sheet):
        for material in MATERIALS:
            primer_rows = [by_key[(stage, material, primer)] for primer in (1, 2, 3)]
            log_values = [
                float(row["log2_target_over_actin_minus_delta_ct"])
                for row in primer_rows
            ]
            summary_rows.append(
                {
                    "stage_index": stage_index,
                    "stage": stage,
                    "phase": "Development" if stage_index <= 13 else "Post-harvest",
                    "material": material,
                    "primer1_n_actin": primer_rows[0]["n_actin"],
                    "primer1_n_target": primer_rows[0]["n_target"],
                    "primer1_log2_target_over_actin": log_values[0],
                    "primer1_target_over_actin": primer_rows[0]["target_over_actin_2pow_minus_delta_ct"],
                    "primer2_n_actin": primer_rows[1]["n_actin"],
                    "primer2_n_target": primer_rows[1]["n_target"],
                    "primer2_log2_target_over_actin": log_values[1],
                    "primer2_target_over_actin": primer_rows[1]["target_over_actin_2pow_minus_delta_ct"],
                    "primer3_n_actin": primer_rows[2]["n_actin"],
                    "primer3_n_target": primer_rows[2]["n_target"],
                    "primer3_log2_target_over_actin": log_values[2],
                    "primer3_target_over_actin": primer_rows[2]["target_over_actin_2pow_minus_delta_ct"],
                    "mean_log2_target_over_actin": statistics.mean(log_values),
                    "geometric_mean_target_over_actin": 2.0 ** statistics.mean(log_values),
                    "provenance_status": "Actin-normalized_raw-grid_reconstruction; plate_map_not_independently_available",
                }
            )
    return long_rows, summary_rows


def calculate_stats(summary_rows: list[dict[str, object]]) -> list[dict[str, object]]:
    values = {
        (str(row["phase"]), str(row["stage"]), str(row["material"])): float(
            row["mean_log2_target_over_actin"]
        )
        for row in summary_rows
    }
    output: list[dict[str, object]] = []
    for phase in ("Development", "Post-harvest"):
        stages = [
            str(row["stage"])
            for row in summary_rows
            if row["phase"] == phase and row["material"] == "TK"
        ]
        comparisons = []
        for comparator in ("TN", "NS", "TK"):
            fl = [values[(phase, stage, "FL")] for stage in stages]
            other = [values[(phase, stage, comparator)] for stage in stages]
            test = wilcoxon(fl, other, alternative="two-sided", method="exact")
            comparisons.append((comparator, float(test.statistic), float(test.pvalue)))
        q_values = bh_adjust([item[2] for item in comparisons])
        for (comparator, statistic, p_value), q_value in zip(comparisons, q_values):
            output.append(
                {
                    "phase": phase,
                    "comparison": f"FL_vs_{comparator}",
                    "n_paired_stages": len(stages),
                    "test": "two-sided exact Wilcoxon signed-rank across matched stages",
                    "statistic": statistic,
                    "p_value": p_value,
                    "multiplicity": "Benjamini-Hochberg across three FL comparisons within phase",
                    "q_value": q_value,
                    "significant_q_lt_0.05": q_value < 0.05,
                    "inference_scope": "stage-level exploratory; not biological-replicate inference",
                }
            )
    return output


def draw(
    summary_rows: list[dict[str, object]],
    png: Path,
    pdf: Path,
) -> None:
    stages = [str(row["stage"]) for row in summary_rows if row["material"] == "TK"]
    by_material = {
        material: [
            float(row["mean_log2_target_over_actin"])
            for row in summary_rows
            if row["material"] == material
        ]
        for material in PLOT_MATERIALS
    }

    fig = plt.figure(figsize=LEGACY_PAGE_SIZE_IN, constrained_layout=False)
    grid = fig.add_gridspec(
        2,
        2,
        height_ratios=(0.83, 1.0),
        left=0.082,
        right=0.992,
        bottom=0.075,
        top=0.955,
        hspace=0.36,
        wspace=0.22,
    )
    ax_top = fig.add_subplot(grid[0, :])
    ax_dev = fig.add_subplot(grid[1, 0])
    ax_post = fig.add_subplot(grid[1, 1], sharey=ax_dev)

    x = np.arange(len(stages))
    for material in PLOT_MATERIALS:
        ax_top.plot(
            x,
            by_material[material],
            color=COLORS[material],
            marker="o",
            markersize=3.6,
            linewidth=1.25,
            markeredgecolor="white",
            markeredgewidth=0.35,
            label=material,
        )
    ax_top.axhline(0, color="#A8ADB4", linewidth=0.7)
    ax_top.axvline(12.5, color="#A8ADB4", linewidth=0.8, linestyle=(0, (2, 2)))
    ax_top.set_xticks(x, stages, rotation=45, ha="right", fontsize=7)
    ax_top.set_ylabel("Mean log2 FAD2/Actin expression (-ΔCt)", fontsize=8)
    ax_top.set_xlabel("Stage", fontsize=8)
    ax_top.tick_params(axis="y", labelsize=7, length=3, width=0.6, color="#555555")
    ax_top.text(
        6,
        1.025,
        "Development",
        transform=ax_top.get_xaxis_transform(),
        ha="center",
        color="#6FJOINT2",
        fontsize=8,
    )
    ax_top.text(
        15.5,
        1.025,
        "Post-harvest",
        transform=ax_top.get_xaxis_transform(),
        ha="center",
        color="#6FJOINT2",
        fontsize=8,
    )
    ax_top.legend(
        frameon=False,
        ncol=4,
        loc="upper left",
        fontsize=7,
        handlelength=1.8,
        columnspacing=1.7,
    )
    ax_top.spines[["top", "right"]].set_visible(False)

    rng = np.random.default_rng(20260831)
    phase_slices = {
        "Development": slice(0, 13),
        "Post-harvest": slice(13, 19),
    }
    all_panel_b_values = np.concatenate(
        [np.asarray(by_material[material], dtype=float) for material in PLOT_MATERIALS]
    )
    panel_b_span = max(float(np.ptp(all_panel_b_values)), 1.0)
    panel_b_limits = (
        min(float(np.min(all_panel_b_values)) - 0.12 * panel_b_span, -0.5),
        float(np.max(all_panel_b_values)) + 0.12 * panel_b_span,
    )
    for axis, phase in ((ax_dev, "Development"), (ax_post, "Post-harvest")):
        stage_slice = phase_slices[phase]
        values = [np.asarray(by_material[material][stage_slice], dtype=float) for material in PLOT_MATERIALS]
        box = axis.boxplot(
            values,
            positions=np.arange(1, 5),
            widths=0.55,
            patch_artist=True,
            showfliers=False,
            medianprops={"color": "#222222", "linewidth": 1.0},
            whiskerprops={"color": "#6FJOINT2", "linewidth": 0.7},
            capprops={"color": "#6FJOINT2", "linewidth": 0.7},
        )
        for patch, material in zip(box["boxes"], PLOT_MATERIALS):
            patch.set(
                facecolor=COLORS[material],
                edgecolor=COLORS[material],
                linewidth=0.8,
                alpha=0.28,
            )
        for index, (material, stage_values) in enumerate(zip(PLOT_MATERIALS, values), start=1):
            jitter = rng.uniform(-0.08, 0.08, size=len(stage_values))
            axis.scatter(
                index + jitter,
                stage_values,
                s=12,
                color=COLORS[material],
                edgecolor="white",
                linewidth=0.25,
                alpha=0.75,
                zorder=3,
            )
        axis.axhline(0, color="#B0B4BA", linewidth=0.6)
        axis.set_xticks(np.arange(1, 5), PLOT_MATERIALS, fontsize=7)
        axis.set_xlabel(phase, fontsize=8)
        axis.tick_params(axis="y", labelsize=7, length=3, width=0.6, color="#555555")
        axis.spines[["top", "right"]].set_visible(False)
        axis.set_ylim(*panel_b_limits)
    ax_dev.set_ylabel("Stage-level mean log2 FAD2/Actin (-ΔCt)", fontsize=8.0)
    ax_post.tick_params(labelleft=False)

    fig.text(0.012, 0.985, "a", fontsize=12, fontweight="bold", va="top")
    fig.text(0.012, 0.505, "b", fontsize=12, fontweight="bold", va="top")
    fig.savefig(png, dpi=600)
    fig.savefig(pdf)
    plt.close(fig)


def csv_workbook_multiset_check() -> tuple[int, list[str]]:
    sheet = load_workbook(WORKBOOK, data_only=True).active
    mismatches = []
    blocks = stage_blocks(sheet)
    for _, stage, start, end in blocks:
        workbook_values = []
        for row in range(start, end + 1):
            if sheet.cell(row, 1).value in MATERIALS:
                for col in range(2, 10):
                    value = sheet.cell(row, col).value
                    if isinstance(value, (int, float)):
                        workbook_values.append(round(float(value), 8))
        csv_values = []
        with (CSV_DIR / f"lxy-{stage}.csv").open(encoding="latin1", newline="") as handle:
            rows = list(csv.reader(handle))
        for row in rows[23:87]:
            if len(row) <= 6:
                continue
            try:
                csv_values.append(round(float(row[6]), 8))
            except ValueError:
                pass
        if Counter(workbook_values) != Counter(csv_values):
            mismatches.append(stage)
    return len(blocks), mismatches


def write_source_workbook(path: Path, sheets: list[tuple[str, pd.DataFrame]]) -> None:
    with pd.ExcelWriter(path, engine="openpyxl") as writer:
        for sheet_name, frame in sheets:
            frame.to_excel(writer, sheet_name=sheet_name, index=False)
            worksheet = writer.book[sheet_name]
            worksheet.freeze_panes = "A2"
            worksheet.auto_filter.ref = worksheet.dimensions
            for column_cells in worksheet.columns:
                width = max(
                    len(str(cell.value)) if cell.value is not None else 0
                    for cell in column_cells
                )
                worksheet.column_dimensions[column_cells[0].column_letter].width = min(
                    max(width + 2, 10), 42
                )


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    long_rows, summary_rows = extract_ddct()
    stats_rows = calculate_stats(summary_rows)
    long_lookup = {
        (str(row["stage"]), str(row["material"]), int(row["primer_set"])): row
        for row in long_rows
    }
    ratio_rows = []
    for primer_set in (1, 2, 3):
        fl_row = long_lookup[("170d", "FL", primer_set)]
        tn_row = long_lookup[("170d", "TN", primer_set)]
        fl_delta_ct = float(fl_row["delta_ct_mean_target_minus_mean_actin"])
        tn_delta_ct = float(tn_row["delta_ct_mean_target_minus_mean_actin"])
        fl = 2.0 ** (-fl_delta_ct)
        tn = 2.0 ** (-tn_delta_ct)
        ratio = 2.0 ** (tn_delta_ct - fl_delta_ct)
        ratio_rows.append(
            {
                "stage": "170d",
                "primer_set": primer_set,
                "calculation": "2^(meanDeltaCt_TN - meanDeltaCt_FL); target and Actin each use four finite 170d wells",
                "FL_mean_target_ct": fl_row["mean_target_ct"],
                "FL_mean_actin_ct": fl_row["mean_actin_ct"],
                "FL_delta_ct": fl_delta_ct,
                "FL_relative_expression": fl,
                "TN_mean_target_ct": tn_row["mean_target_ct"],
                "TN_mean_actin_ct": tn_row["mean_actin_ct"],
                "TN_delta_ct": tn_delta_ct,
                "TN_relative_expression": tn,
                "FL_over_TN_ratio": ratio,
                "manuscript_rounded_ratio": f"{ratio:.2f}",
                "n_target_wells_FL": fl_row["n_target"],
                "n_actin_wells_FL": fl_row["n_actin"],
                "n_target_wells_TN": tn_row["n_target"],
                "n_actin_wells_TN": tn_row["n_actin"],
                "provenance_status": "Actin-normalized_workbook_raw_grid; plate_map_not_independently_available",
            }
        )

    png = OUT / "SupplementaryFig13_actin_delta_ct_final_consistent_600dpi.png"
    pdf = OUT / "SupplementaryFig13_actin_delta_ct_final_consistent.pdf"
    draw(summary_rows, png, pdf)

    n_csv_stages, csv_mismatches = csv_workbook_multiset_check()
    actin_counts = Counter(int(row["n_actin"]) for row in long_rows)
    target_counts = Counter(int(row["n_target"]) for row in long_rows)
    expected_ratios = (6.976483362247222, 3.951772788562064, 8.267779999193406)
    validation_rows = [
        {"check": "stage_count", "observed": len({row["stage"] for row in long_rows}), "expected": 19, "pass": len({row["stage"] for row in long_rows}) == 19, "note": "13 developmental + 6 post-harvest"},
        {"check": "stage_material_count", "observed": len(summary_rows), "expected": 76, "pass": len(summary_rows) == 76, "note": "19 stages x 4 materials"},
        {"check": "stage_material_primer_count", "observed": len(long_rows), "expected": 228, "pass": len(long_rows) == 228, "note": "19 x 4 x 3"},
        {"check": "finite_actin_well_counts", "observed": f"n4={actin_counts[4]};n3={actin_counts[3]};n2={actin_counts[2]}", "expected": "n4=228;n3=0;n2=0", "pass": actin_counts == Counter({4: 228}), "note": "Actin mean uses all finite Actin wells"},
        {"check": "finite_target_well_counts", "observed": f"n4={target_counts[4]};n3={target_counts[3]};n2={target_counts[2]}", "expected": "n4=204;n3=21;n2=3", "pass": target_counts == Counter({4: 204, 3: 21, 2: 3}), "note": "target mean uses all finite target wells; no imputation"},
        {"check": "csv_to_workbook_raw_ct_multiset", "observed": f"{n_csv_stages - len(csv_mismatches)}/{n_csv_stages}", "expected": "19/19", "pass": not csv_mismatches, "note": f"mismatch stages: {','.join(csv_mismatches) if csv_mismatches else 'none'}"},
        {"check": "manual_Q_U_values_used", "observed": False, "expected": False, "pass": True, "note": "reconstruction uses raw-grid A:I Ct values only"},
        {"check": "170d_ratio_primer1", "observed": ratio_rows[0]["FL_over_TN_ratio"], "expected": expected_ratios[0], "pass": math.isclose(float(ratio_rows[0]["FL_over_TN_ratio"]), expected_ratios[0], abs_tol=1e-12), "note": "rounds to 6.98"},
        {"check": "170d_ratio_primer2", "observed": ratio_rows[1]["FL_over_TN_ratio"], "expected": expected_ratios[1], "pass": math.isclose(float(ratio_rows[1]["FL_over_TN_ratio"]), expected_ratios[1], abs_tol=1e-12), "note": "rounds to 3.95"},
        {"check": "170d_ratio_primer3", "observed": ratio_rows[2]["FL_over_TN_ratio"], "expected": expected_ratios[2], "pass": math.isclose(float(ratio_rows[2]["FL_over_TN_ratio"]), expected_ratios[2], abs_tol=1e-12), "note": "rounds to 8.27"},
        {"check": "significance_annotations_in_figure", "observed": False, "expected": False, "pass": True, "note": "audit-only Wilcoxon/BH table is not used as a manuscript conclusion"},
        {"check": "plate_map_available", "observed": False, "expected": False, "pass": True, "note": "known limitation: sample/target/replicate identities cannot be independently reconciled"},
    ]

    long_df = pd.DataFrame(long_rows)
    summary_df = pd.DataFrame(summary_rows)
    ratio_df = pd.DataFrame(ratio_rows)
    validation_df = pd.DataFrame(validation_rows)
    stats_df = pd.DataFrame(stats_rows)

    long_path = OUT / "SupplementaryFig13_actin_delta_ct_plotdata_long.tsv"
    summary_path = OUT / "SupplementaryFig13_actin_delta_ct_stage_summary.tsv"
    stats_path = OUT / "AUDIT_ONLY_SupplementaryFig13_stage_level_wilcoxon_bh.tsv"
    ratio_path = OUT / "SupplementaryFig13_170d_FL_TN_actin_delta_ct_ratios.tsv"
    validation_path = OUT / "SupplementaryFig13_actin_delta_ct_validation.tsv"
    long_df.to_csv(long_path, sep="\t", index=False)
    summary_df.to_csv(summary_path, sep="\t", index=False)
    stats_df.to_csv(stats_path, sep="\t", index=False)
    ratio_df.to_csv(ratio_path, sep="\t", index=False)
    validation_df.to_csv(validation_path, sep="\t", index=False)

    source_xlsx = OUT / "SupplementaryFig13_source_data_Actin_normalized.xlsx"
    write_source_workbook(
        source_xlsx,
        [
            ("raw_ct_and_deltaCt", long_df),
            ("stage_summary", summary_df),
            ("170d_ratios", ratio_df),
            ("validation", validation_df),
        ],
    )

    readme = OUT / "README.md"
    readme.write_text(
        "# Supplementary Figure 13: Actin-normalized reconstruction\n\n"
        "This directory is a new, non-overwriting reconstruction from the raw-grid Ct values "
        "in ZIP 3 `数据整理.xlsx`. It does not use the workbook's manually assembled Q:U "
        "values. For each stage, material and FAD2 primer assay, `DeltaCt = mean Ct(FAD2) - "
        "mean Ct(Actin)`, and the plotted log2 target/reference value is `-DeltaCt`. Target and "
        "Actin means are calculated separately because no plate map establishes one-to-one "
        "biological pairing between their wells. There is no stage-specific calibrator, which "
        "preserves temporal variation. Non-numeric/no-Ct target values were excluded without "
        "imputation; all Actin groups retained four finite wells, while 204 target groups "
        "retained four, 21 retained three and 3 retained two.\n\n"
        "At 170 d, all FL and TN assays retain four finite target and four finite Actin wells. The ratios recalculated from the workbook raw grids "
        "FL/TN ratios are 6.976483, 3.951773 and 8.267780, which support manuscript values "
        "6.98, 3.95 and 8.27. The alternative 11.51, 5.17 and 9.98 values are ratios of the "
        "workbook W:Y arithmetic means and are not used here because W:Y averages target-only "
        "relative values without dividing by the Actin column.\n\n"
        "Panel b shows the distributions of the same 13 developmental or 6 post-harvest "
        "stage-level aggregate values. It contains no significance brackets or inferential claim. "
        "The Wilcoxon/BH TSV is retained only as an audit trace and must not be treated as "
        "biological-replicate inference or a manuscript conclusion without a verified plate map.\n\n"
        "Recoverable instrument metadata: all 19 CSV exports report a qTOWER3/G instrument "
        "(device ID 3107B-3522352), SYBR Green acquisition, an initial 95 C hold for 2 min, "
        "40 cycles of 95 C for 5 s, 60 C for 10 s and 72 C for 15 s with acquisition, followed "
        "by an initial melt acquisition at 60 C and 35 one-degree increments to 95 C "
        "(36 acquisition points). Primer sequences, master-mix identity, "
        "reaction volume and assay efficiencies are not present in the supplied archives and "
        "must be added from the wet-lab record before submission.\n\n"
        "Limitation: every instrument CSV labels the wells as unknown samples, and no plate map "
        "was supplied. The 19 CSV Ct multisets match the workbook raw grids exactly, but material, "
        "target, biological-replicate and technical-replicate identities cannot be independently "
        "reconciled from the instrument exports. This limitation must remain in the source-data "
        "description or Methods.\n\n"
        "Suggested Results sentence: `Workbook-assigned qPCR measurements with three FAD2 "
        "primer sets produced Actin-normalized FL/TN expression ratios of 6.98, 3.95 and "
        "8.27 at 170 d, respectively (Supplementary Fig. 13).`\n\n"
        "Suggested caption: `Supplementary Fig. 13 | Workbook-assigned FAD2 qPCR measurements. "
        "a, Mean log2 FAD2/Actin expression across three primer assays, "
        "calculated as -DeltaCt, where DeltaCt is mean Ct(FAD2) minus mean Ct(Actin). "
        "Target and Actin means used all finite wells separately; two to four finite target "
        "wells and four finite Actin wells were available per stage-material-assay combination, "
        "and no-Ct target wells were excluded without imputation. "
        "b, Distributions of the same stage-level means during development (13 stages) and "
        "post-harvest storage (6 stages). Boxes show the interquartile range with median lines, "
        "whiskers extend to 1.5 times the interquartile range, and points denote stages.`\n",
        encoding="utf-8",
    )

    outputs = [
        long_path,
        summary_path,
        stats_path,
        ratio_path,
        validation_path,
        source_xlsx,
        png,
        pdf,
        readme,
    ]
    manifest = OUT / "sha256.tsv"
    with manifest.open("w", encoding="utf-8", newline="") as handle:
        handle.write("sha256\tbytes\tfile\n")
        for path in outputs:
            handle.write(f"{sha256(path)}\t{path.stat().st_size}\t{path.name}\n")

    print(f"workbook={WORKBOOK}")
    print(f"long_rows={len(long_rows)}")
    print(f"summary_rows={len(summary_rows)}")
    print(f"output_dir={OUT}")


if __name__ == "__main__":
    main()
