#!/usr/bin/env python3
"""Rebuild audited Figure 2e and 2f panels without modifying source files.

Figure 2e expands the representative-metabolite heatmap from seven sentinel
stages to all 19 measured stages and states the chemical evidence tier.

Figure 2f uses the current Astral-114 directLFQ matrix for protein-detection
summaries. Its FAD2 RNA track is restricted to the confirmed chr08 locus so
the RNA and protein scopes are aligned. A separate comparison panel makes the
three incompatible FAD2 protein sources explicit; it never labels legacy
timsTOF/OAU values as Astral.
"""

from __future__ import annotations

import hashlib
import math
import platform
import re
import sys
from collections import defaultdict
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap, Normalize, TwoSlopeNorm
from matplotlib.patches import Rectangle


OUT_ROOT = Path(__file__).resolve().parents[1]
FIG_DIR = OUT_ROOT / "figures"
TAB_DIR = OUT_ROOT / "tables"
LOG_DIR = OUT_ROOT / "logs"

BASE = Path("${ANALYSIS_DIR}")
FIG2_SOURCE = (
    BASE
    / "22_answer_reviews/00_ms/03_V3/02_figure"
    / "05_Figure2_evolution_multiomics_panels_flat_20260727"
)
ATLAS = (
    BASE
    / "22_answer_reviews/00_ms/03_V3/03_figure3/05_multiomics_integration"
    / "runs/RUN-MULTIOMICS-INTEGRATION-20260721-001/outputs"
    / "stage26_full_chemical_phenotype_atlas_attempt002"
)
PROT_ROOT = (
    BASE
    / "22_answer_reviews/00_ms/03_V3/03_figure3"
    / "07_proteomics_multiomics_20260724"
)
REF_DIR = PROT_ROOT / "01_reference_database"
ASTRAL_RUN = (
    PROT_ROOT
    / "02_current_Astral_114/runs"
    / "RUN-PROT-ASTRAL-DIRECTLFQ-20260725-001"
)
ASTRAL_MATRIX = (
    ASTRAL_RUN
    / "outputs/current114_directlfq_protein_abundance_integration_keys.tsv"
)
LEGACY_PANEL_SOURCE = FIG2_SOURCE / "Fig2G_FL_TN_push_pull_package_protect_source.tsv"
LEGACY_OAU_MATRIX = (
    BASE
    / "22_answer_reviews/00_ms/03_V3/02_figure"
    / "06_Fig2_key_FA_enzyme_heatmap_panels_20260728"
    / "Protein_validated_groups_enzyme_log2_FL_TN_19.tsv"
)
FA_GENE_LIST = BASE / "21_MS/03_result/01_omic/FA_gene_list.tsv"
RNA_FA_VALUES = BASE / "21_MS/03_result/01_omic/FA_RNA_TPM_pergene.tsv"
AFRICA_ANNOTATION = (
    BASE
    / "20_results/Figure2/07_new_figure/05_omic/1.final_counts"
    / "GO_annotation/Africa_hap2/Africa_hap2.emapper.annotations"
)
OAU_METRICS = (
    BASE
    / "22_answer_reviews/00_ms/03_V3/03_figure3/05_multiomics_integration"
    / "runs/RUN-MULTIOMICS-ALLELE-CNS-V4-20260723-001/outputs"
    / "stage28_proteome_interface_attempt001/six_genome_protein_to_OAU_metrics.tsv.gz"
)
UNIFIED_MAP = REF_DIR / "FL_TN_unified_exact_sequence_nr_map.tsv"

SCORES = ATLAS / "molecular_phenotype_scores_by_sample.tsv"
COVERAGE = ATLAS / "phenotype_axis_chemical_coverage_and_status.tsv"
COMPOUND_MATRIX = (
    ATLAS.parent / "stage13_compound_abundance_attempt003/compound_consensus_z_matrix.tsv"
)
SELECTED_METABOLITES = FIG2_SOURCE / "Fig2E_FL_TN_marker_metabolite_heatmap_source.tsv"

STAGES = [
    "0d", "15d", "35d", "50d", "65d", "80d", "95d", "110d", "125d",
    "140d", "155d", "170d", "185d", "12h", "24h", "36h", "48h", "60h", "72h",
]
GENOTYPES = ["FL", "TN"]
# Match the final assembled Figure 2: FL is teal and TN is coral red.
GENOTYPE_COLOURS = {"FL": "#168A8A", "TN": "#D35F4B"}

FAD2_CHR08_UNIFIED_IDS = {
    "UFTN046834",  # FL chr08B: evm.TU.chr08B.792
    "UFTN046835",  # FL chr08A: evm.TU.chr08A.741
    "UFTN046836",  # TN hap1: evm.TU.bk_hap1_chr2.783
    "UFTN046837",  # TN hap2: evm.TU.bk_hap2_chr2.1260 (not quantified)
}
FAD2_CHR08_RNA_GENE_ID = "evm.TU.chr08B.792"

AXES = [
    ("P02", "Storage lipids", "#C8902F"),
    ("P01", "Oleic balance", "#4A8F78"),
    ("P03", "Hydrolytic rancidity", "#C65F4B"),
    ("P04", "Oxidative rancidity", "#7B67A4"),
]

MARKERS = [
    ("Push", "ACCase", "precursor input", "ACCase"),
    ("Push", "KASIII", "FA initiation", "KASIII"),
    ("Push", "ENR", "FA synthesis", "FabI (ENR)"),
    ("Push", "FATA/B", "FA export", "FATA/B"),
    ("Pull", "LACS", "acyl activation", "LACS"),
    ("Pull", "GPAT", "glycerolipid initiation", "GPAT"),
    ("Pull", "FAD2-like", "omega-6 desaturation", "FAD2"),
    ("Package / protect", "DGAT", "TAG deposition", "DGAT"),
    ("Package / protect", "OLE16", "lipid-droplet surface", "OLE16"),
    ("Package / protect", "LOX (loss risk)", "lipid oxidation", "LOX (loss risk)"),
]
PROCESS_COLOURS = {
    "Push": "#D68B3B",
    "Pull": "#3C9180",
    "Package / protect": "#517CA8",
}


def configure_plotting() -> None:
    mpl.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Liberation Sans", "DejaVu Sans"],
            "font.size": 8.0,
            "axes.titlesize": 9.2,
            "axes.labelsize": 8.0,
            "xtick.labelsize": 6.8,
            "ytick.labelsize": 7.0,
            "axes.linewidth": 0.7,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
            "savefig.facecolor": "white",
        }
    )


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def save_figure(fig: plt.Figure, stem: str) -> None:
    fig.savefig(FIG_DIR / f"{stem}.pdf", bbox_inches="tight", pad_inches=0.05)
    fig.savefig(
        FIG_DIR / f"{stem}.png",
        dpi=600,
        bbox_inches="tight",
        pad_inches=0.05,
    )
    plt.close(fig)


def stage_positions() -> dict[str, float]:
    positions = {}
    for i, stage in enumerate(STAGES):
        positions[stage] = float(i + (0.75 if i >= 13 else 0.0))
    return positions


def draw_figure2e() -> dict[str, int]:
    scores = pd.read_csv(SCORES, sep="\t")
    coverage = pd.read_csv(COVERAGE, sep="\t")
    axis_ids = [axis_id for axis_id, _, _ in AXES]
    selected_scores = scores[
        scores["axis_id"].isin(axis_ids) & scores["genotype"].isin(GENOTYPES)
    ].copy()

    expected_cells = set(
        (axis_id, genotype, stage)
        for axis_id in axis_ids
        for genotype in GENOTYPES
        for stage in STAGES
    )
    observed_cells = set(
        selected_scores[["axis_id", "genotype", "stage"]].itertuples(index=False, name=None)
    )
    if expected_cells != observed_cells:
        missing = sorted(expected_cells - observed_cells)
        raise ValueError(f"Figure 2e score cells are incomplete: {missing[:10]}")

    replicate_counts = selected_scores.groupby(
        ["axis_id", "genotype", "stage"], observed=True
    ).size()
    if not replicate_counts.eq(3).all():
        raise ValueError("Figure 2e requires exactly three biological replicates per cell")

    summary = (
        selected_scores.groupby(
            ["axis_id", "genotype", "stage", "stage_index", "phase"],
            as_index=False,
            observed=True,
        )["molecular_phenotype_score"]
        .agg(mean="mean", sd="std", n="size")
    )
    summary["se"] = summary["sd"] / np.sqrt(summary["n"])

    coverage_small = coverage[coverage["axis_id"].isin(axis_ids)][
        [
            "axis_id",
            "axis_name_en",
            "level2_compound_n",
            "level3_compound_n",
            "ms1_mass_proxy_compound_n",
            "signed_compound_n",
            "final_axis_status",
        ]
    ].copy()
    observed_coverage = {
        row.axis_id: (int(row.level2_compound_n), int(row.ms1_mass_proxy_compound_n))
        for row in coverage_small.itertuples(index=False)
    }
    expected_coverage = {
        "P01": (1, 32),
        "P02": (0, 11),
        "P03": (3, 66),
        "P04": (7, 10),
    }
    if observed_coverage != expected_coverage:
        raise AssertionError(
            f"Figure 2e chemical-evidence coverage changed: {observed_coverage}"
        )
    coverage_small.to_csv(TAB_DIR / "Figure2e_axis_coverage.tsv", sep="\t", index=False)
    summary.merge(coverage_small, on="axis_id", how="left").to_csv(
        TAB_DIR / "Figure2e_axis_summary_19stage.tsv", sep="\t", index=False
    )

    old_selected = pd.read_csv(SELECTED_METABOLITES, sep="\t")
    selected_meta = (
        old_selected[
            [
                "axis_id",
                "compound_group_id",
                "display_name",
                "component_role",
                "evidence_level",
            ]
        ]
        .drop_duplicates()
        .copy()
    )
    # Preserve the row order used in the final assembled Figure 2.
    selected_meta = selected_meta.reset_index(drop=True)
    if len(selected_meta) != 8:
        raise ValueError(f"Expected eight frozen representative metabolites, found {len(selected_meta)}")

    matrix = pd.read_csv(COMPOUND_MATRIX, sep="\t").set_index("compound_group_id")
    rows = []
    for meta in selected_meta.itertuples(index=False):
        if meta.compound_group_id not in matrix.index:
            raise KeyError(f"Missing compound in full matrix: {meta.compound_group_id}")
        for genotype in GENOTYPES:
            for stage_index, stage in enumerate(STAGES, start=1):
                columns = [f"{genotype}|{stage}|R{i}" for i in (1, 2, 3)]
                values = pd.to_numeric(matrix.loc[meta.compound_group_id, columns], errors="coerce")
                rows.append(
                    {
                        "axis_id": meta.axis_id,
                        "compound_group_id": meta.compound_group_id,
                        "display_name": str(meta.display_name).replace("†", ""),
                        "component_role": meta.component_role,
                        "evidence_level": meta.evidence_level,
                        "genotype": genotype,
                        "stage_index": stage_index,
                        "stage": stage,
                        "replicate_n": int(values.notna().sum()),
                        "mean_consensus_z": float(values.mean()),
                    }
                )
    metabolite_long = pd.DataFrame(rows)
    metabolite_long.to_csv(
        TAB_DIR / "Figure2e_representative_metabolites_19stage.tsv", sep="\t", index=False
    )
    if len(metabolite_long) != 8 * 2 * 19:
        raise AssertionError("Figure 2e full-stage metabolite table has the wrong size")
    low_coverage = metabolite_long[metabolite_long["replicate_n"].lt(3)].copy()
    observed_low_coverage = {
        (row.display_name, row.genotype, row.stage): int(row.replicate_n)
        for row in low_coverage.itertuples(index=False)
    }
    expected_low_coverage = {
        ("TG 42:0", "FL", "95d"): 1,
        ("TG 42:0", "TN", "36h"): 0,
    }
    if observed_low_coverage != expected_low_coverage:
        raise AssertionError(
            f"Figure 2e low-coverage cells changed: {observed_low_coverage}"
        )

    configure_plotting()
    # The final Figure 2 assembled the original 10.8 x 2.65 trajectory strip
    # above the 8.8 x 3.72 metabolite table.  This combined canvas retains that
    # visual ratio while accommodating all 19 measured stages.
    fig = plt.figure(figsize=(10.8, 7.15), facecolor="white")
    outer = fig.add_gridspec(
        2,
        1,
        height_ratios=[2.55, 3.72],
        left=0.055,
        right=0.995,
        top=0.94,
        bottom=0.145,
        hspace=0.22,
    )
    trajectory_grid = outer[0].subgridspec(1, 4, wspace=0.14)
    xpos = stage_positions()
    harvest_x = (xpos["185d"] + xpos["12h"]) / 2
    coverage_index = coverage_small.set_index("axis_id")
    rng = np.random.default_rng(20260728)

    for col, (axis_id, title, axis_colour) in enumerate(AXES):
        ax = fig.add_subplot(trajectory_grid[0, col])
        ax.axvspan(harvest_x, max(xpos.values()) + 0.55, color="#F3F5F5", zorder=0)
        ax.axvline(harvest_x, color="#A8B1B4", lw=0.65, ls=(0, (2.5, 2.5)), zorder=1)
        ax.axhline(0, color="#CED5D7", lw=0.55, zorder=1)
        ax.grid(axis="y", color="#E8ECEC", lw=0.5, zorder=0)
        q = selected_scores[selected_scores["axis_id"].eq(axis_id)].copy()
        for offset, genotype in zip((-0.08, 0.08), GENOTYPES):
            raw = q[q["genotype"].eq(genotype)].sort_values(
                ["stage_index", "replicate"]
            )
            jitter = rng.normal(offset, 0.024, len(raw))
            raw_x = raw["stage"].map(xpos).to_numpy(float) + jitter
            ax.scatter(
                raw_x,
                raw["molecular_phenotype_score"],
                s=6,
                color=GENOTYPE_COLOURS[genotype],
                alpha=0.24,
                edgecolor="none",
                zorder=2,
            )
            qq = summary[
                summary["axis_id"].eq(axis_id) & summary["genotype"].eq(genotype)
            ].sort_values("stage_index")
            xx = qq["stage"].map(xpos).to_numpy(float)
            yy = qq["mean"].to_numpy(float)
            se = qq["se"].to_numpy(float)
            ax.fill_between(
                xx,
                yy - se,
                yy + se,
                color=GENOTYPE_COLOURS[genotype],
                alpha=0.10,
                linewidth=0,
            )
            ax.plot(
                xx,
                yy,
                color=GENOTYPE_COLOURS[genotype],
                lw=1.35,
                marker="o",
                ms=2.2,
                markeredgecolor="white",
                markeredgewidth=0.25,
                label=genotype,
                zorder=3,
            )
        cov = coverage_index.loc[axis_id]
        subtitle = (
            f"putative axis: {int(cov.level2_compound_n)} L2 + "
            f"{int(cov.ms1_mass_proxy_compound_n)} MS1 proxies"
        )
        ax.set_title(title, loc="left", color="#263136", fontweight="semibold", fontsize=8.3, pad=13)
        ax.text(0.0, 1.015, subtitle, transform=ax.transAxes, ha="left", va="bottom", fontsize=5.4, color=axis_colour)
        ax.set_xlim(-0.55, max(xpos.values()) + 0.55)
        ax.set_ylim(-0.82, 0.86)
        shown = ["0d", "125d", "185d", "72h"]
        ax.set_xticks([xpos[x] for x in shown], shown, rotation=40, ha="right")
        ax.tick_params(length=2.0, width=0.55, color="#778185", labelsize=6.2)
        ax.spines[["top", "right"]].set_visible(False)
        ax.spines[["left", "bottom"]].set_color("#8A9497")
        ax.spines[["left", "bottom"]].set_linewidth(0.6)
        if col == 0:
            ax.set_ylabel("Score", fontsize=7.0)
        else:
            ax.tick_params(labelleft=False)
        if col == 3:
            ax.legend(frameon=False, loc="upper right", fontsize=6.2, ncol=2,
                      handlelength=1.6, columnspacing=0.8, borderaxespad=0.1)

    heat_grid = outer[1].subgridspec(1, 2, width_ratios=[0.245, 0.755], wspace=0.018)
    names_ax = fig.add_subplot(heat_grid[0, 0])
    heat_ax = fig.add_subplot(heat_grid[0, 1])
    ordered_ids = selected_meta["compound_group_id"].tolist()
    blocks = []
    for genotype in GENOTYPES:
        block = (
            metabolite_long[metabolite_long["genotype"].eq(genotype)]
            .pivot(index="compound_group_id", columns="stage", values="mean_consensus_z")
            .reindex(index=ordered_ids, columns=STAGES)
            .to_numpy(float)
        )
        blocks.append(block)
    values = np.concatenate([blocks[0], np.full((len(ordered_ids), 1), np.nan), blocks[1]], axis=1)
    cmap = LinearSegmentedColormap.from_list("metabolite", ["#177E89", "#F7F4EA", "#CC5A4A"])
    norm = TwoSlopeNorm(vmin=-2.2, vcenter=0, vmax=2.2)
    group_gap = 0.62
    y_positions = np.asarray([i + group_gap * (i // 2) for i in range(len(ordered_ids))], dtype=float)
    ymax = y_positions[-1] + 0.75
    names_ax.set_xlim(0, 1)
    names_ax.set_ylim(ymax, -0.80)
    names_ax.axis("off")
    heat_ax.set_xlim(-0.5, values.shape[1] - 0.5)
    heat_ax.set_ylim(ymax, -0.80)
    heat_ax.spines[:].set_visible(False)
    heat_ax.set_yticks([])
    for i, y in enumerate(y_positions):
        for j in range(values.shape[1]):
            if j == 19:
                continue
            face = cmap(norm(values[i, j])) if np.isfinite(values[i, j]) else "#E5E8E9"
            heat_ax.add_patch(
                Rectangle(
                    (j - 0.46, y - 0.43),
                    0.92,
                    0.86,
                    facecolor=face,
                    edgecolor="#FFFFFF",
                    linewidth=0.75,
                )
            )
    group_colours = {axis_id: colour for axis_id, _, colour in AXES}
    group_titles = {
        "P02": "Storage-lipid related",
        "P01": "C18:1 related",
        "P03": "Hydrolysis related",
        "P04": "Oxidation related",
    }
    for i, row in selected_meta.iterrows():
        y = y_positions[i]
        colour = group_colours[row["axis_id"]]
        label = str(row["display_name"]).replace("†", "")
        if row["evidence_level"] != "Level_2_candidate":
            label += "†"
        names_ax.add_patch(Rectangle((0.00, y - 0.43), 0.016, 0.86, facecolor=colour, edgecolor="none"))
        names_ax.text(0.035, y, label, ha="left", va="center", fontsize=7.2, color="#273136")
    for axis_id, *_ in AXES:
        indices = np.flatnonzero(selected_meta["axis_id"].eq(axis_id).to_numpy())
        if len(indices):
            y = y_positions[indices[0]] - 0.62
            names_ax.text(0.00, y, group_titles[axis_id], ha="left", va="bottom",
                          fontsize=7.5, color=group_colours[axis_id], fontweight="bold")

    heat_ax.set_xticks(
        list(range(19)) + list(range(20, 39)),
        STAGES + STAGES,
        rotation=48,
        ha="right",
    )
    heat_ax.tick_params(length=0, pad=3, colors="#566267", labelsize=5.1)
    for start, genotype in ((0, "FL"), (20, "TN")):
        colour = GENOTYPE_COLOURS[genotype]
        heat_ax.add_patch(Rectangle((start - 0.46, -0.72), 18.92, 0.33,
                                    facecolor=mpl.colors.to_rgba(colour, 0.14), edgecolor="none"))
        heat_ax.text(start + 9.0, -0.555, genotype, ha="center", va="center",
                     fontsize=8.2, color=colour, fontweight="bold")
        heat_ax.axvline(start + 12.5, color="#7E8B8F", lw=0.7, ls=(0, (3, 3)))
    heat_ax.axvline(19.5, color="#536166", lw=1.1)
    row_positions = {compound_id: i for i, compound_id in enumerate(ordered_ids)}
    for cell in low_coverage.itertuples(index=False):
        x = (cell.stage_index - 1) + (0 if cell.genotype == "FL" else 20)
        y = y_positions[row_positions[cell.compound_group_id]]
        label = "NA" if cell.replicate_n == 0 else f"n={cell.replicate_n}"
        heat_ax.add_patch(
            Rectangle(
                (x - 0.47, y - 0.47),
                0.94,
                0.94,
                fill=False,
                edgecolor="#22282A",
                linewidth=0.75,
            )
        )
        heat_ax.text(
            x,
            y,
            label,
            ha="center",
            va="center",
            fontsize=4.6,
            fontweight="bold",
            color="#202628",
        )

    cax = fig.add_axes([0.780, 0.045, 0.170, 0.025])
    cbar = mpl.colorbar.ColorbarBase(cax, cmap=cmap, norm=norm, orientation="horizontal", ticks=[-2, 0, 2])
    cbar.ax.tick_params(labelsize=6.2, length=2, pad=1)
    cbar.outline.set_visible(False)
    cax.set_title("relative abundance (z-score)", fontsize=6.2, pad=2, color="#5B676B")
    fig.text(
        0.055,
        0.052,
        "† putative MS1 proxy; unmarked = Level 2; boxed n=1/NA = one/no measured replicate.",
        ha="left",
        va="center",
        fontsize=5.9,
        color="#4E5A5F",
    )
    save_figure(fig, "Figure2e_full19_putative_MS1_redraw")
    return {
        "axis_rows": len(summary),
        "metabolite_rows": len(metabolite_long),
        "metabolite_low_coverage_cells": len(low_coverage),
    }


def read_africa_annotation_gene_sets() -> dict[str, set[str]]:
    annotation = pd.read_csv(
        AFRICA_ANNOTATION,
        sep="\t",
        comment="#",
        header=None,
        dtype=str,
        low_memory=False,
    )
    description = annotation.get(7, pd.Series("", index=annotation.index)).fillna("")
    preferred = annotation.get(8, pd.Series("", index=annotation.index)).fillna("")
    searchable = description.str.cat(preferred, sep=" ")
    ole = searchable.str.contains(r"oleosin 16|\bole16\b", case=False, regex=True, na=False)
    lox = searchable.str.contains(r"plant lipoxygenase|\blox9\b", case=False, regex=True, na=False)
    return {
        "OLE16": set(annotation.loc[ole, 0].astype(str)),
        "LOX (loss risk)": set(annotation.loc[lox, 0].astype(str)),
    }


def build_astral_family_assignments() -> tuple[pd.DataFrame, pd.DataFrame]:
    genes = pd.read_csv(FA_GENE_LIST, sep="\t", dtype=str)
    family_gene_ids: dict[str, set[str]] = {}
    for _, _, _, source_family in MARKERS:
        if source_family in {"OLE16", "LOX (loss risk)"}:
            continue
        family_gene_ids[source_family] = set(
            genes.loc[genes["Enzyme"].eq(source_family), "GeneID"].dropna().astype(str)
        )
    family_gene_ids.update(read_africa_annotation_gene_sets())

    oau = pd.read_csv(
        OAU_METRICS,
        sep="\t",
        usecols=["allele_unit_id", "gene_id", "genome_role"],
        dtype=str,
        low_memory=False,
    )
    gene_to_oauses = oau.groupby("gene_id", observed=True)["allele_unit_id"].agg(lambda x: set(x.dropna()))
    oau_family_labels: dict[str, set[str]] = defaultdict(set)
    for family, gene_ids in family_gene_ids.items():
        for gene_id in gene_ids:
            for oau_id in gene_to_oauses.get(gene_id, set()):
                oau_family_labels[oau_id].add(family)
    unambiguous_oau_family = {
        oau_id: next(iter(labels))
        for oau_id, labels in oau_family_labels.items()
        if len(labels) == 1
    }

    unified = pd.read_csv(
        UNIFIED_MAP,
        sep="\t",
        usecols=["unified_protein_id", "source_gene_ids"],
        dtype=str,
    )
    unified_family = {}
    mapping_rows = []
    for row in unified.itertuples(index=False):
        gene_ids = []
        for token in str(row.source_gene_ids).split(";"):
            gene_ids.append(token.split(":", 1)[1] if ":" in token else token)
        oau_ids = set()
        for gene_id in gene_ids:
            oau_ids.update(gene_to_oauses.get(gene_id, set()))
        labels = {unambiguous_oau_family[x] for x in oau_ids if x in unambiguous_oau_family}
        family = next(iter(labels)) if len(labels) == 1 else None
        if family == "FAD2" and row.unified_protein_id not in FAD2_CHR08_UNIFIED_IDS:
            family = None
        if family is not None:
            unified_family[row.unified_protein_id] = family
        mapping_rows.append(
            {
                "unified_protein_id": row.unified_protein_id,
                "source_gene_n": len(gene_ids),
                "mapped_oau_n": len(oau_ids),
                "family_label_n": len(labels),
                "assigned_family": family,
            }
        )

    matrix = pd.read_csv(ASTRAL_MATRIX, sep="\t", low_memory=False)
    group_rows = []
    assigned = []
    for group in matrix["unified_protein_group"].astype(str):
        members = group.split(";")
        labels = {unified_family[x] for x in members if x in unified_family}
        family = next(iter(labels)) if len(labels) == 1 else None
        group_rows.append(
            {
                "unified_protein_group": group,
                "member_n": len(members),
                "mapped_family_n": len(labels),
                "assigned_family": family,
            }
        )
        assigned.append(family)
    matrix.insert(1, "assigned_family", assigned)
    mapping_audit = pd.DataFrame(mapping_rows)
    mapping_audit[mapping_audit["assigned_family"].notna()].to_csv(
        TAB_DIR / "Figure2f_Astral114_unified_protein_assignment.tsv", sep="\t", index=False
    )
    group_audit = pd.DataFrame(group_rows)
    group_audit[group_audit["assigned_family"].notna()].to_csv(
        TAB_DIR / "Figure2f_Astral114_group_assignment.tsv", sep="\t", index=False
    )
    fad2_map = unified[unified["unified_protein_id"].isin(FAD2_CHR08_UNIFIED_IDS)].copy()
    quantified_members = {
        member
        for group in matrix["unified_protein_group"].astype(str)
        for member in group.split(";")
    }
    fad2_map["scope"] = "chr08 FAD2 orthogroup"
    fad2_map["quantified_in_Astral114"] = fad2_map["unified_protein_id"].isin(
        quantified_members
    )
    fad2_map.to_csv(
        TAB_DIR / "Figure2f_Astral114_chr08_FAD2_mapping.tsv", sep="\t", index=False
    )
    return matrix, group_audit


def aggregate_astral_markers(matrix: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    records = []
    replicate_records = []
    source_families = [row[3] for row in MARKERS]
    for source_family in source_families:
        subset = matrix[matrix["assigned_family"].eq(source_family)].copy()
        numeric = subset.drop(columns=["unified_protein_group", "assigned_family"]).apply(pd.to_numeric, errors="coerce")
        for genotype in GENOTYPES:
            for stage_index, stage in enumerate(STAGES, start=1):
                columns = [f"{genotype}|{stage}|R{i}" for i in (1, 2, 3)]
                replicate_values = []
                for replicate, column in enumerate(columns, start=1):
                    values = numeric[column] if column in numeric else pd.Series(dtype=float)
                    positive = values[np.isfinite(values) & values.gt(0)]
                    abundance = float(positive.sum()) if len(positive) else math.nan
                    replicate_values.append(abundance)
                    replicate_records.append(
                        {
                            "family": source_family,
                            "genotype": genotype,
                            "stage_index": stage_index,
                            "stage": stage,
                            "replicate": f"R{replicate}",
                            "detected": int(np.isfinite(abundance) and abundance > 0),
                            "summed_directLFQ_abundance": abundance,
                            "assigned_protein_group_n": len(subset),
                        }
                    )
                arr = np.asarray(replicate_values, dtype=float)
                detected_replicates = int(np.isfinite(arr).sum())
                records.append(
                    {
                        "family": source_family,
                        "genotype": genotype,
                        "stage_index": stage_index,
                        "stage": stage,
                        "detected_replicate_n": detected_replicates,
                        "detected_any_replicate": int(detected_replicates >= 1),
                        "detected_at_least_two_replicates": int(detected_replicates >= 2),
                        "mean_detected_replicate_abundance": float(np.nanmean(arr)) if detected_replicates else math.nan,
                        "assigned_protein_group_n": len(subset),
                    }
                )
    stage_table = pd.DataFrame(records)
    replicate_table = pd.DataFrame(replicate_records)
    stage_table.to_csv(TAB_DIR / "Figure2f_Astral114_stage_detection.tsv", sep="\t", index=False)
    replicate_table.to_csv(TAB_DIR / "Figure2f_Astral114_replicate_detection.tsv", sep="\t", index=False)
    return stage_table, replicate_table


def _within_track_z(values: np.ndarray) -> np.ndarray:
    """Match the legacy Figure 2 RNA-track transform for one 19-stage row."""
    values = np.asarray(values, dtype=float)
    transformed = np.log2(values + 1.0)
    valid = np.isfinite(transformed)
    out = np.full(values.shape, np.nan, dtype=float)
    if valid.sum() < 2:
        return out
    sd = float(np.nanstd(transformed[valid], ddof=0))
    if sd == 0:
        out[valid] = 0.0
    else:
        out[valid] = (
            transformed[valid] - np.nanmean(transformed[valid])
        ) / sd
    return np.clip(out, -2.5, 2.5)


def build_figure2f_rna_source() -> pd.DataFrame:
    """Retain legacy RNA tracks except for the chr08-only FAD2 replacement."""
    rna = pd.read_csv(LEGACY_PANEL_SOURCE, sep="\t")
    rna = rna[rna["omics"].eq("RNA")].copy()
    plotted_markers = {marker for _process, marker, _function, _family in MARKERS}
    rna = rna[rna["marker"].isin(plotted_markers)].copy()

    source = pd.read_csv(RNA_FA_VALUES, sep="\t")
    selected = source[
        source["Geneid"].eq(FAD2_CHR08_RNA_GENE_ID)
        & source["Enzyme"].eq("FAD2")
    ]
    if len(selected) != 1:
        raise AssertionError(
            f"Expected one chr08 FAD2 RNA row, found {len(selected)}"
        )
    selected_row = selected.iloc[0]

    for genotype in GENOTYPES:
        columns = [f"{genotype}{i:02d}" for i in range(1, 20)]
        values = pd.to_numeric(selected_row[columns], errors="coerce").to_numpy(float)
        if not np.isfinite(values).all():
            raise AssertionError(f"Non-finite chr08 FAD2 RNA value for {genotype}")
        selector = rna["marker"].eq("FAD2-like") & rna["genotype"].eq(genotype)
        target_index = rna.loc[selector].sort_values("stage_index").index
        if len(target_index) != 19:
            raise AssertionError(
                f"Expected 19 FAD2 RNA stages for {genotype}, found {len(target_index)}"
            )
        rna.loc[target_index, "raw_aggregate"] = values
        rna.loc[target_index, "within_track_zscore"] = _within_track_z(values)
        rna.loc[target_index, "matched_feature_n"] = 1
        rna.loc[target_index, "matched_features"] = FAD2_CHR08_RNA_GENE_ID

    fad2 = rna[rna["marker"].eq("FAD2-like")]
    ratio_170d = (
        fad2.query("genotype == 'FL' and stage == '170d'")["raw_aggregate"].iloc[0]
        / fad2.query("genotype == 'TN' and stage == '170d'")["raw_aggregate"].iloc[0]
    )
    if not math.isclose(ratio_170d, 4.387189342467, rel_tol=1e-12):
        raise AssertionError(f"chr08 FAD2 RNA 170d ratio changed: {ratio_170d}")

    rna.to_csv(
        TAB_DIR / "Figure2f_RNA_chr08_aligned_source.tsv", sep="\t", index=False
    )
    return rna


def build_legacy_oau_stage_table() -> pd.DataFrame:
    matrix = pd.read_csv(LEGACY_OAU_MATRIX, sep="\t", na_values=["NA"])
    records = []
    available_families = set(matrix["Enzyme"].astype(str))
    for _, _, _, source_family in MARKERS:
        available = source_family in available_families
        row = matrix[matrix["Enzyme"].eq(source_family)].iloc[0] if available else None
        for genotype in GENOTYPES:
            for stage_index, stage in enumerate(STAGES, start=1):
                value = (
                    pd.to_numeric(pd.Series([row[f"{genotype}{stage_index:02d}"]]), errors="coerce").iloc[0]
                    if available
                    else math.nan
                )
                records.append(
                    {
                        "family": source_family,
                        "genotype": genotype,
                        "stage_index": stage_index,
                        "stage": stage,
                        "detected_replicate_n": math.nan,
                        "detected_any_replicate": int(np.isfinite(value)) if available else 0,
                        "detected_at_least_two_replicates": math.nan,
                        "mean_detected_replicate_abundance": value,
                        "assigned_protein_group_n": math.nan,
                        "source_available": available,
                    }
                )
    out = pd.DataFrame(records)
    out.to_csv(
        TAB_DIR / "Figure2f_AUDIT_ONLY_legacy_timsTOF_OAU_stage_detection.tsv",
        sep="\t",
        index=False,
    )
    return out


def draw_figure2f(
    stage_table: pd.DataFrame,
    rna: pd.DataFrame,
    *,
    stem: str,
    panel_label: str,
    summary_name: str,
    footer: str,
    detection_column: str,
    detection_rule: str,
    audit_only: bool = False,
) -> dict[str, int]:
    rows = [
        (process, marker, genotype, function, source_family)
        for process, marker, function, source_family in MARKERS
        for genotype in GENOTYPES
    ]
    group_gap = 0.70
    y_positions = []
    for i, (process, *_rest) in enumerate(rows):
        group_index = 0 if process == "Push" else (1 if process == "Pull" else 2)
        y_positions.append(i + group_index * group_gap)
    y_positions = np.asarray(y_positions)

    configure_plotting()
    # Reuse the original Figure 2f/G canvas and table geometry so this panel can
    # replace the old asset without changing the assembled figure's proportions.
    fig, ax = plt.subplots(figsize=(14.0, 7.8), facecolor="white")
    ax.set_xlim(0, 61.0)
    ax.set_ylim(y_positions[-1] + 1.15, -3.10)
    ax.axis("off")
    x_marker, x_variety, x_function = 1.7, 8.3, 10.8
    matrix_x, cell_w = 20.0, 1.48
    matrix_end = matrix_x + 19 * cell_w
    summary_x = matrix_end + 1.0
    rna_peak_x = summary_x + 5.2
    cmap = LinearSegmentedColormap.from_list(
        "rna_programme", ["#FFFFFF", "#D6EDF7", "#7CC2E5", "#247CAE"]
    )
    norm = Normalize(-2, 2)

    development_end = matrix_x + 13 * cell_w
    ax.add_patch(Rectangle((matrix_x, -2.65), development_end - matrix_x - 0.18, 0.82, facecolor="#CDEAF8", edgecolor="none"))
    ax.add_patch(Rectangle((development_end + 0.15, -2.65), matrix_end - development_end - 0.15, 0.82, facecolor="#CDEAF8", edgecolor="none"))
    ax.text((matrix_x + development_end) / 2, -2.24, "Development (0–185 d)", ha="center", va="center", fontsize=8.2, color="#263136", fontweight="bold")
    ax.text((development_end + matrix_end) / 2, -2.24, "Post-harvest (12–72 h)", ha="center", va="center", fontsize=8.2, color="#263136", fontweight="bold")
    ax.text(x_marker, -1.15, "Gene family", fontsize=7.1, color="#58656A", fontweight="bold")
    ax.text(x_variety, -1.15, "Variety", fontsize=7.1, color="#58656A", fontweight="bold")
    ax.text(x_function, -1.15, "Function", fontsize=7.1, color="#58656A", fontweight="bold")
    ax.text(
        summary_x + 3.3,
        -2.24,
        "Summary",
        ha="center",
        va="center",
        fontsize=8.2,
        color="#263136",
        fontweight="bold",
    )
    ax.text(
        summary_x + 1.15,
        -1.15,
        "Protein det.",
        ha="center",
        va="center",
        fontsize=6.6,
        color="#58656A",
        fontweight="bold",
    )
    ax.text(
        rna_peak_x,
        -1.15,
        "RNA peak",
        ha="center",
        va="center",
        fontsize=6.6,
        color="#58656A",
        fontweight="bold",
    )
    for j, stage in enumerate(STAGES):
        ax.text(matrix_x + (j + 0.5) * cell_w, -0.95, stage, rotation=48, ha="left", va="bottom", fontsize=5.8, color="#465358")

    for process in PROCESS_COLOURS:
        indices = [i for i, row in enumerate(rows) if row[0] == process]
        y0 = y_positions[min(indices)] - 0.43
        y1 = y_positions[max(indices)] + 0.43
        ax.add_patch(Rectangle((0, y0), 1.15, y1 - y0, facecolor=PROCESS_COLOURS[process], edgecolor="none"))
        ax.text(0.575, (y0 + y1) / 2, process.replace(" / ", "/"), rotation=90, ha="center", va="center", color="white", fontweight="bold", fontsize=7.0)

    summary_rows = []
    for i, (process, marker, genotype, function, source_family) in enumerate(rows):
        y = y_positions[i]
        family_index = i // 2
        if family_index % 2 == 0:
            ax.add_patch(Rectangle((1.35, y - 0.46), 64.0, 0.92, facecolor="#F2F3F3", edgecolor="none", zorder=0))
        marker_label = "FAD2 (chr08)" if source_family == "FAD2" and not audit_only else marker
        ax.text(x_marker, y, marker_label, ha="left", va="center", fontsize=6.8, fontweight="semibold", color="#253035")
        ax.text(x_variety, y, genotype, ha="left", va="center", fontsize=6.8, fontweight="bold", color=GENOTYPE_COLOURS[genotype])
        ax.text(x_function, y, function, ha="left", va="center", fontsize=6.3, color="#4C585D")

        qrna = rna[
            rna["marker"].eq(marker) & rna["genotype"].eq(genotype)
        ].sort_values("stage_index")
        if len(qrna) != 19:
            raise ValueError(f"RNA source incomplete for {marker} {genotype}: {len(qrna)}")
        rna_z = qrna["within_track_zscore"].to_numpy(float)
        for j, zvalue in enumerate(rna_z):
            face = cmap(norm(zvalue)) if np.isfinite(zvalue) else "#E1E3E4"
            ax.add_patch(Rectangle((matrix_x + j * cell_w + 0.08, y - 0.39), cell_w - 0.16, 0.78, facecolor=face, edgecolor="#566267", lw=0.32))

        qpro = stage_table[
            stage_table["family"].eq(source_family) & stage_table["genotype"].eq(genotype)
        ].sort_values("stage_index")
        if len(qpro) != 19:
            raise ValueError(f"Astral source incomplete for {source_family} {genotype}: {len(qpro)}")
        available = bool(qpro["source_available"].all()) if "source_available" in qpro else True
        selected_detection = pd.to_numeric(qpro[detection_column], errors="coerce")
        detected_n = int(selected_detection.fillna(0).sum()) if available else 0
        detected_any_n = int(qpro["detected_any_replicate"].sum()) if available else 0
        strict_values = pd.to_numeric(qpro["detected_at_least_two_replicates"], errors="coerce")
        detected_strict_n = int(strict_values.sum()) if strict_values.notna().any() else math.nan
        fraction = detected_n / 19 if available else 0.0
        bar_x, bar_w = summary_x, 2.30
        ax.add_patch(Rectangle((bar_x, y - 0.22), bar_w, 0.44, facecolor="#E4E6E7", edgecolor="#AAB2B5", lw=0.35))
        if available:
            ax.add_patch(Rectangle((bar_x, y - 0.22), bar_w * fraction, 0.44, facecolor="#777F83", edgecolor="none"))
            percent = f"{fraction * 100:.0f}%" if math.isclose(fraction, 1.0) else f"{fraction * 100:.1f}%"
            detection_label = percent
        else:
            detection_label = "not OAU-mapped"
        ax.text(bar_x + bar_w + 0.18, y, detection_label, ha="left", va="center", fontsize=5.7, color="#59656A")
        raw_values = qrna["raw_aggregate"].to_numpy(float)
        peak = STAGES[int(np.nanargmax(raw_values))] if np.isfinite(raw_values).any() else "ND"
        ax.text(rna_peak_x, y, peak, ha="center", va="center", fontsize=6.1, color="#384348")
        summary_rows.append(
            {
                "process": process,
                "marker": marker,
                "source_family": source_family,
                "genotype": genotype,
                "detection_source": panel_label.replace("\n", "; "),
                "source_available": available,
                "stage_detected_n": detected_n if available else math.nan,
                "stage_detected_any_replicate_n": detected_any_n if available else math.nan,
                "stage_detected_at_least_two_replicates_n": detected_strict_n,
                "stage_total_n": 19,
                "detection_percent": fraction * 100 if available else math.nan,
                "detection_rule": detection_rule,
                "RNA_peak_stage": peak,
            }
        )

    summary = pd.DataFrame(summary_rows)
    summary.to_csv(TAB_DIR / summary_name, sep="\t", index=False)
    cax = fig.add_axes([0.355, 0.035, 0.16, 0.020])
    mpl.colorbar.ColorbarBase(cax, cmap=cmap, norm=norm, orientation="horizontal", ticks=[-2, 0, 2])
    cax.tick_params(labelsize=5.5, length=2, pad=1)
    cax.set_title("within-track z-score", fontsize=5.8, pad=2, color="#59656A")
    fig.text(
        0.535,
        0.045,
        "RNA family-level abundance",
        fontsize=5.9,
        color="#5F6B70",
        va="center",
    )
    if audit_only:
        ax.text(
            -0.012,
            1.012,
            "AUDIT ONLY - legacy source; do not use as the manuscript replacement",
            transform=ax.transAxes,
            ha="left",
            va="top",
            fontsize=7.0,
            color="#B53A3A",
            fontweight="bold",
        )
    save_figure(fig, stem)
    return {
        "marker_summary_rows": len(summary),
        "fad2_fl": int(summary.query("marker == 'FAD2-like' and genotype == 'FL'")["stage_detected_n"].iloc[0]),
        "fad2_tn": int(summary.query("marker == 'FAD2-like' and genotype == 'TN'")["stage_detected_n"].iloc[0]),
    }


def fad2_source_comparison(stage_table: pd.DataFrame) -> pd.DataFrame:
    legacy = pd.read_csv(LEGACY_PANEL_SOURCE, sep="\t")
    legacy = legacy[(legacy["marker"].eq("FAD2-like")) & (legacy["omics"].eq("Protein"))]
    old_counts = (
        legacy.assign(detected=legacy["raw_aggregate"].gt(0).astype(int))
        .groupby("genotype", observed=True)["detected"]
        .sum()
        .to_dict()
    )

    strict = pd.read_csv(LEGACY_OAU_MATRIX, sep="\t", na_values=["NA"])
    strict_row = strict[strict["Enzyme"].eq("FAD2")].iloc[0]
    strict_counts = {
        genotype: int(
            pd.to_numeric(strict_row[[f"{genotype}{i:02d}" for i in range(1, 20)]], errors="coerce").notna().sum()
        )
        for genotype in GENOTYPES
    }
    astral_any_counts = (
        stage_table[stage_table["family"].eq("FAD2")]
        .groupby("genotype", observed=True)["detected_any_replicate"]
        .sum()
        .astype(int)
        .to_dict()
    )
    astral_strict_counts = (
        stage_table[stage_table["family"].eq("FAD2")]
        .groupby("genotype", observed=True)["detected_at_least_two_replicates"]
        .sum()
        .astype(int)
        .to_dict()
    )
    rows = []
    source_specs = [
        ("Final panel source", "legacy exact-gene MaxLFQ", old_counts),
        ("Strict OAU rebuild", "legacy timsTOF directLFQ/OAU", strict_counts),
        ("Astral any replicate", "Astral-114 directLFQ >=1/3", astral_any_counts),
        (
            "Astral recommended",
            "Astral-114 directLFQ >=2/3",
            astral_strict_counts,
        ),
    ]
    for source_label, provenance, counts in source_specs:
        for genotype in GENOTYPES:
            detected_n = int(counts.get(genotype, 0))
            rows.append(
                {
                    "source": source_label,
                    "provenance": provenance,
                    "genotype": genotype,
                    "detected_stage_n": detected_n,
                    "total_stage_n": 19,
                    "detected_percent": detected_n / 19 * 100,
                }
            )
    comparison = pd.DataFrame(rows)
    comparison.to_csv(TAB_DIR / "Figure2f_FAD2_source_comparison.tsv", sep="\t", index=False)

    configure_plotting()
    fig, ax = plt.subplots(figsize=(8.6, 3.9), facecolor="white")
    source_order = [x[0] for x in source_specs]
    x = np.arange(len(source_order), dtype=float)
    width = 0.32
    for offset, genotype in zip((-width / 2, width / 2), GENOTYPES):
        q = comparison[comparison["genotype"].eq(genotype)].set_index("source").loc[source_order]
        bars = ax.bar(x + offset, q["detected_percent"], width=width, color=GENOTYPE_COLOURS[genotype], label=genotype)
        for bar, row in zip(bars, q.itertuples()):
            ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 1.3, f"{row.detected_stage_n}/19", ha="center", va="bottom", fontsize=7.0)
    ax.set_ylabel("FAD2 protein detection across stages (%)")
    ax.set_xticks(
        x,
        [
            "Final figure\nlegacy MaxLFQ",
            "Legacy timsTOF\nstrict OAU",
            "Astral-114\n>=1/3 replicates",
            "Astral-114\n>=2/3 replicates",
        ],
    )
    ax.set_ylim(0, 90)
    ax.grid(axis="y", color="#E2E7E8", lw=0.6)
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(frameon=False, ncol=2, loc="upper right")
    ax.set_title("FAD2 detection depends on the proteomics source", loc="left", fontweight="bold")
    fig.text(
        0.11,
        0.015,
        "AUDIT ONLY. Legacy strict-OAU values (36.8% FL; 42.1% TN) are not Astral-114 results; the recommended Astral rule is >=2/3 replicates.",
        fontsize=6.6,
        color="#9E2F2F",
        fontweight="bold",
    )
    fig.subplots_adjust(left=0.12, right=0.98, top=0.88, bottom=0.25)
    save_figure(fig, "Figure2f_AUDIT_ONLY_FAD2_source_comparison")
    return comparison


def write_provenance(inputs: list[Path], stats: dict[str, int], comparison: pd.DataFrame) -> None:
    manifest = pd.DataFrame(
        [{"path": str(path), "sha256": sha256(path), "size_bytes": path.stat().st_size} for path in inputs]
    )
    manifest.to_csv(LOG_DIR / "fig2_redraw_inputs.sha256.tsv", sep="\t", index=False)
    comparison_text = comparison.to_csv(sep="\t", index=False)
    (LOG_DIR / "Figure2_source_decision.txt").write_text(
        "Figure 2 source decision\n"
        "========================\n"
        "Figure 2e: all four upper trajectories and the replacement lower heatmap use all 19 measured stages.\n"
        "P02 Storage lipids has 0 Level-2 compounds and 11 putative MS1 proxies.\n"
        "P01 Oleic balance has 1 Level-2 compound and 32 putative MS1 proxies.\n"
        "No values were interpolated.\n\n"
        "Figure 2f: 36.8% FL and 42.1% TN come from a legacy timsTOF strict-OAU rebuild, not Astral.\n"
        "The RNA heatmap uses only evm.TU.chr08B.792 for FAD2; its 170d FL/TN ratio is 4.387189.\n"
        "The Astral-114 mapping is restricted to the chr08 FAD2 orthogroup: UFTN046834/046835 (FL),\n"
        "UFTN046836/046837 (TN); UFTN046837 was not quantified. At >=1/3 replicates, detection is\n"
        "15/19 FL and 12/19 TN. The recommended >=2/3 rule gives 14/19 FL and 11/19 TN.\n"
        "The >=2/3 Astral result is used in Figure2f_Astral114_RECOMMENDED_redraw.\n\n"
        + comparison_text,
        encoding="utf-8",
    )
    (LOG_DIR / "fig2_redraw.log").write_text(
        "status\tPASS\n"
        f"python\t{sys.version.split()[0]}\n"
        f"platform\t{platform.platform()}\n"
        f"pandas\t{pd.__version__}\n"
        f"matplotlib\t{mpl.__version__}\n"
        + "".join(f"{key}\t{value}\n" for key, value in sorted(stats.items())),
        encoding="utf-8",
    )
    output_files = [
        FIG_DIR / f"{stem}.{suffix}"
        for stem in (
            "Figure2e_full19_putative_MS1_redraw",
            "Figure2f_Astral114_RECOMMENDED_redraw",
            "Figure2f_AUDIT_ONLY_legacy_timsTOF_OAU_redraw",
            "Figure2f_AUDIT_ONLY_FAD2_source_comparison",
        )
        for suffix in ("pdf", "png")
    ]
    output_files.extend(
        [
            TAB_DIR / "Figure2e_axis_coverage.tsv",
            TAB_DIR / "Figure2e_axis_summary_19stage.tsv",
            TAB_DIR / "Figure2e_representative_metabolites_19stage.tsv",
            TAB_DIR / "Figure2f_Astral114_chr08_FAD2_mapping.tsv",
            TAB_DIR / "Figure2f_Astral114_marker_summary.tsv",
            TAB_DIR / "Figure2f_Astral114_replicate_detection.tsv",
            TAB_DIR / "Figure2f_Astral114_stage_detection.tsv",
            TAB_DIR / "Figure2f_RNA_chr08_aligned_source.tsv",
            TAB_DIR / "Figure2f_FAD2_source_comparison.tsv",
        ]
    )
    pd.DataFrame(
        [
            {
                "path": str(path),
                "sha256": sha256(path),
                "size_bytes": path.stat().st_size,
            }
            for path in output_files
        ]
    ).to_csv(LOG_DIR / "fig2_redraw_outputs.sha256.tsv", sep="\t", index=False)


def main() -> None:
    for directory in (FIG_DIR, TAB_DIR, LOG_DIR):
        directory.mkdir(parents=True, exist_ok=True)
    inputs = [
        SCORES,
        COVERAGE,
        COMPOUND_MATRIX,
        SELECTED_METABOLITES,
        LEGACY_PANEL_SOURCE,
        LEGACY_OAU_MATRIX,
        FA_GENE_LIST,
        RNA_FA_VALUES,
        AFRICA_ANNOTATION,
        OAU_METRICS,
        UNIFIED_MAP,
        ASTRAL_MATRIX,
        ASTRAL_RUN / "qc/directlfq_validation_summary.json",
    ]
    missing = [path for path in inputs if not path.is_file()]
    if missing:
        raise FileNotFoundError("Missing frozen inputs: " + ", ".join(map(str, missing)))

    stats = {}
    stats.update({f"fig2e_{key}": value for key, value in draw_figure2e().items()})
    rna = build_figure2f_rna_source()
    fad2_rna = rna[rna["marker"].eq("FAD2-like")]
    stats["fig2f_rna_rows"] = len(rna)
    stats["fig2f_fad2_rna_feature_n"] = int(
        fad2_rna["matched_features"].nunique()
    )
    stats["fig2f_fad2_rna_170d_ratio"] = float(
        fad2_rna.query("genotype == 'FL' and stage == '170d'")["raw_aggregate"].iloc[0]
        / fad2_rna.query("genotype == 'TN' and stage == '170d'")["raw_aggregate"].iloc[0]
    )
    astral_matrix, group_audit = build_astral_family_assignments()
    fad2_groups = set(
        group_audit.loc[
            group_audit["assigned_family"].eq("FAD2"), "unified_protein_group"
        ]
    )
    expected_fad2_groups = {"UFTN046834;UFTN046835", "UFTN046836"}
    if fad2_groups != expected_fad2_groups:
        raise AssertionError(f"Astral chr08 FAD2 groups changed: {fad2_groups}")
    stage_table, replicate_table = aggregate_astral_markers(astral_matrix)
    stats["astral_assigned_group_n"] = int(group_audit["assigned_family"].notna().sum())
    stats["astral_stage_rows"] = len(stage_table)
    stats["astral_replicate_rows"] = len(replicate_table)
    fig2f_stats = draw_figure2f(
        stage_table,
        rna,
        stem="Figure2f_Astral114_RECOMMENDED_redraw",
        panel_label="Astral-114 directLFQ\nFAD2: chr08 orthogroup",
        summary_name="Figure2f_Astral114_marker_summary.tsv",
        footer="Protein stage detected when >=2 of 3 biological replicates have nonzero directLFQ abundance; no imputation.",
        detection_column="detected_at_least_two_replicates",
        detection_rule=">=2/3 biological replicates with nonzero directLFQ abundance",
    )
    stats.update({f"fig2f_{key}": value for key, value in fig2f_stats.items()})
    legacy_stage_table = build_legacy_oau_stage_table()
    legacy_stats = draw_figure2f(
        legacy_stage_table,
        rna,
        stem="Figure2f_AUDIT_ONLY_legacy_timsTOF_OAU_redraw",
        panel_label="Legacy timsTOF strict-OAU\nrebuild",
        summary_name="Figure2f_AUDIT_ONLY_legacy_timsTOF_OAU_marker_summary.tsv",
        footer="Legacy stage detection is based on finite strict-OAU aggregate abundance; OLE16/LOX were not mapped in this rebuild.",
        detection_column="detected_any_replicate",
        detection_rule="finite legacy strict-OAU aggregate abundance",
        audit_only=True,
    )
    stats.update({f"fig2f_legacy_{key}": value for key, value in legacy_stats.items()})
    comparison = fad2_source_comparison(stage_table)

    observed = comparison.pivot(index="provenance", columns="genotype", values="detected_stage_n")
    expected = {
        "legacy exact-gene MaxLFQ": {"FL": 10, "TN": 7},
        "legacy timsTOF directLFQ/OAU": {"FL": 7, "TN": 8},
        "Astral-114 directLFQ >=1/3": {"FL": 15, "TN": 12},
        "Astral-114 directLFQ >=2/3": {"FL": 14, "TN": 11},
    }
    for provenance, genotype_counts in expected.items():
        for genotype, expected_count in genotype_counts.items():
            count = int(observed.loc[provenance, genotype])
            if count != expected_count:
                raise AssertionError(
                    f"FAD2 source audit changed for {provenance} {genotype}: {count} != {expected_count}"
                )
    if (fig2f_stats["fad2_fl"], fig2f_stats["fad2_tn"]) != (14, 11):
        raise AssertionError(
            "Recommended Astral FAD2 panel must report 14/19 FL and 11/19 TN"
        )
    write_provenance(inputs, stats, comparison)
    print("PASS: wrote audited Figure 2e/2f redraws")
    print(f"Output root: {OUT_ROOT}")


if __name__ == "__main__":
    main()
