#!/usr/bin/env python3
"""Redraw Figure 2f with the joint-normalized 114-sample RNA data."""

from __future__ import annotations

import hashlib
import importlib.util
import math
from pathlib import Path

import numpy as np
import pandas as pd
from PIL import Image


HERE = Path(__file__).resolve().parent
AUDIT = HERE / "ST9_latest_reanalysis_audit"
ORIGINAL_SCRIPT = (
    HERE.parent
    / "Revised_Panels_20260831/scripts/fig2_redraw.py"
)
ORIGINAL_FIGURE = (
    HERE.parent
    / "Revised_Panels_20260831/figures/Figure2f_Astral114_RECOMMENDED_redraw.pdf"
)
LEGACY_RNA_LAYOUT = (
    Path("${ANALYSIS_DIR}")
    / "22_answer_reviews/00_ms/03_V3/02_figure"
    / "05_Figure2_evolution_multiomics_panels_flat_20260727"
    / "Fig2G_FL_TN_push_pull_package_protect_source.tsv"
)
RAW_COUNTS = (
    Path("${ANALYSIS_DIR}")
    / "22_answer_reviews/00_ms/03_V3/03_figure3/00_minipan"
    / "03_rnaseq_mapping/count_matrices/pangraphrna_hisat2_graph.gene_counts.tsv"
)
CROSSWALK = (
    Path("${ANALYSIS_DIR}")
    / "22_answer_reviews/00_ms/03_V3/03_figure3/05_multiomics_integration/runs"
    / "RUN-MULTIOMICS-INTEGRATION-20260721-001/outputs/stage1_identity"
    / "sample_crosswalk_114.tsv"
)
SIZE_FACTORS = AUDIT / "RNA_114_joint_size_factors.tsv"
ST9_RNA = AUDIT / "ST9_RNA_pathway_stage_means_latest.tsv"
ASTRAL_STAGE = (
    HERE.parent
    / "Revised_Panels_20260831/tables/Figure2f_Astral114_stage_detection.tsv"
)

STEM = "Figure2f_LATEST114_RECOMMENDED_redraw"
RNA_OUTPUT = AUDIT / "Figure2f_RNA_latest114_joint_source.tsv"
SUMMARY_OUTPUT = "Figure2f_LATEST114_marker_summary.tsv"
QA_OUTPUT = AUDIT / "Figure2f_LATEST114_QA.tsv"
MANIFEST_OUTPUT = AUDIT / "Figure2f_LATEST114_input_sha256.tsv"
REFERENCE_GENE_ADDITIONS = {
    "OLE16": ["evm.TU.chr03B.2346"],
    "LOX (loss risk)": ["evm.TU.chr01B.2493"],
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_plotter():
    spec = importlib.util.spec_from_file_location("original_fig2_redraw", ORIGINAL_SCRIPT)
    if spec is None or spec.loader is None:
        raise RuntimeError("Cannot load the original Figure 2 redraw script")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    module.FIG_DIR = HERE
    module.TAB_DIR = AUDIT
    module.LOG_DIR = AUDIT
    return module


def parse_features(value: str) -> list[str]:
    return [token.strip() for token in str(value).split(";") if token.strip()]


def marker_gene_map(plotter, legacy: pd.DataFrame) -> dict[str, list[str]]:
    mapping = {}
    for _process, marker, _function, _family in plotter.MARKERS:
        if marker == "FAD2-like":
            mapping[marker] = [plotter.FAD2_CHR08_RNA_GENE_ID]
            continue
        rows = legacy[
            legacy["marker"].eq(marker)
            & legacy["genotype"].eq("FL")
            & legacy["omics"].eq("RNA")
        ]
        feature_sets = {
            tuple(parse_features(value))
            for value in rows["matched_features"].dropna().astype(str)
        }
        if len(feature_sets) != 1:
            raise ValueError(f"Expected one FL feature set for {marker}: {feature_sets}")
        mapping[marker] = list(next(iter(feature_sets)))
        for gene in REFERENCE_GENE_ADDITIONS.get(marker, []):
            if gene not in mapping[marker]:
                mapping[marker].append(gene)
    return mapping


def build_latest_rna(plotter) -> tuple[pd.DataFrame, dict[str, list[str]]]:
    legacy = pd.read_csv(LEGACY_RNA_LAYOUT, sep="\t")
    plotted_markers = {marker for _process, marker, _function, _family in plotter.MARKERS}
    rna = legacy[
        legacy["omics"].eq("RNA")
        & legacy["marker"].isin(plotted_markers)
    ].copy()
    mapping = marker_gene_map(plotter, legacy)
    required_genes = sorted({gene for genes in mapping.values() for gene in genes})

    counts = pd.read_csv(RAW_COUNTS, sep="\t", index_col=0)
    missing_genes = sorted(set(required_genes) - set(counts.index))
    if missing_genes:
        raise KeyError(f"RNA count matrix is missing genes: {missing_genes}")

    crosswalk = pd.read_csv(CROSSWALK, sep="\t", dtype=str)
    size_factors = pd.read_csv(SIZE_FACTORS, sep="\t")
    if len(crosswalk) != 114 or crosswalk["integration_key"].nunique() != 114:
        raise ValueError("RNA crosswalk is not 114 unique samples")
    if len(size_factors) != 114 or size_factors["sample"].nunique() != 114:
        raise ValueError("DESeq2 size-factor table is not 114 unique samples")

    factor = size_factors.set_index("sample")["DESeq2_size_factor"].astype(float)
    sample_order = crosswalk["rna_sample"].tolist()
    missing_samples = sorted(set(sample_order) - set(counts.columns))
    missing_factors = sorted(set(sample_order) - set(factor.index))
    if missing_samples or missing_factors:
        raise KeyError(
            f"Missing RNA samples={missing_samples[:5]} size_factors={missing_factors[:5]}"
        )
    normalized = counts.loc[required_genes, sample_order].astype(float).div(
        factor.loc[sample_order], axis=1
    )

    for marker, genes in mapping.items():
        family_by_sample = normalized.loc[genes].sum(axis=0)
        for genotype in plotter.GENOTYPES:
            values = []
            for stage in plotter.STAGES:
                samples = crosswalk.loc[
                    crosswalk["genotype"].eq(genotype)
                    & crosswalk["stage"].eq(stage),
                    "rna_sample",
                ].tolist()
                if len(samples) != 3:
                    raise ValueError(f"Expected three RNA samples for {genotype} {stage}")
                values.append(float(family_by_sample.loc[samples].mean()))
            selector = rna["marker"].eq(marker) & rna["genotype"].eq(genotype)
            target_index = rna.loc[selector].sort_values("stage_index").index
            if len(target_index) != 19:
                raise ValueError(f"RNA layout is incomplete for {marker} {genotype}")
            values_array = np.asarray(values, dtype=float)
            rna.loc[target_index, "raw_aggregate"] = values_array
            rna.loc[target_index, "within_track_zscore"] = plotter._within_track_z(values_array)
            rna.loc[target_index, "matched_feature_n"] = len(genes)
            rna.loc[target_index, "matched_features"] = ";".join(genes)
            rna.loc[target_index, "normalization"] = "joint_DESeq2_114_samples"

    rna = rna.sort_values(["process", "marker", "genotype", "stage_index"])
    if len(rna) != len(plotter.MARKERS) * len(plotter.GENOTYPES) * len(plotter.STAGES):
        raise ValueError(f"Unexpected RNA plot row count: {len(rna)}")
    RNA_OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    rna.to_csv(RNA_OUTPUT, sep="\t", index=False)
    return rna, mapping


def validate_against_st9(
    rna: pd.DataFrame,
    mapping: dict[str, list[str]],
) -> tuple[int, float]:
    st9 = pd.read_csv(ST9_RNA, sep="\t")
    checked = 0
    max_abs_difference = 0.0
    st9_genes = set(st9["Gene_ID"])
    for marker, genes in mapping.items():
        if not set(genes).issubset(st9_genes):
            continue
        for genotype in ("FL", "TN"):
            for stage in (
                "0d", "15d", "35d", "50d", "65d", "80d", "95d", "110d",
                "125d", "140d", "155d", "170d", "185d", "12h", "24h", "36h",
                "48h", "60h", "72h",
            ):
                expected = st9.loc[
                    st9["Gene_ID"].isin(genes)
                    & st9["Material"].eq(genotype)
                    & st9["Stage"].eq(stage),
                    "Mean_normalized_RNA_count",
                ].sum()
                observed = rna.loc[
                    rna["marker"].eq(marker)
                    & rna["genotype"].eq(genotype)
                    & rna["stage"].eq(stage),
                    "raw_aggregate",
                ].iloc[0]
                difference = abs(float(observed) - float(expected))
                max_abs_difference = max(max_abs_difference, difference)
                if not math.isclose(float(observed), float(expected), rel_tol=1e-12, abs_tol=1e-8):
                    raise AssertionError(
                        f"Figure 2f/ST9 RNA mismatch for {marker} {genotype} {stage}: "
                        f"{observed} != {expected}"
                    )
                checked += 1
    return checked, max_abs_difference


def main() -> None:
    inputs = [
        ORIGINAL_SCRIPT, ORIGINAL_FIGURE, LEGACY_RNA_LAYOUT, RAW_COUNTS,
        CROSSWALK, SIZE_FACTORS, ST9_RNA, ASTRAL_STAGE,
    ]
    missing = [path for path in inputs if not path.is_file()]
    if missing:
        raise FileNotFoundError(f"Missing Figure 2f input(s): {missing}")

    plotter = load_plotter()
    rna, mapping = build_latest_rna(plotter)
    st9_checked, st9_max_difference = validate_against_st9(rna, mapping)
    stage_table = pd.read_csv(ASTRAL_STAGE, sep="\t")
    if len(stage_table) != 380:
        raise ValueError(f"Astral family-stage table has {len(stage_table)} rows")

    stats = plotter.draw_figure2f(
        stage_table,
        rna,
        stem=STEM,
        panel_label="Astral-114 directLFQ\nFAD2: chr08 orthogroup",
        summary_name=SUMMARY_OUTPUT,
        footer=(
            "Protein stage detected when >=2 of 3 biological replicates have nonzero "
            "directLFQ abundance; no imputation."
        ),
        detection_column="detected_at_least_two_replicates",
        detection_rule=">=2/3 biological replicates with nonzero directLFQ abundance",
    )

    fad2 = rna[rna["marker"].eq("FAD2-like")]
    fad2_fl_170 = float(
        fad2.query("genotype == 'FL' and stage == '170d'")["raw_aggregate"].iloc[0]
    )
    fad2_tn_170 = float(
        fad2.query("genotype == 'TN' and stage == '170d'")["raw_aggregate"].iloc[0]
    )
    fad2_ratio = fad2_fl_170 / fad2_tn_170
    if not math.isclose(fad2_ratio, 0.7302864126, rel_tol=1e-8):
        raise AssertionError(f"Unexpected latest FAD2 RNA ratio: {fad2_ratio}")
    if (stats["fad2_fl"], stats["fad2_tn"]) != (14, 11):
        raise AssertionError(f"Unexpected Astral FAD2 detection: {stats}")

    png = HERE / f"{STEM}.png"
    pdf = HERE / f"{STEM}.pdf"
    with Image.open(png) as image:
        rgb = np.asarray(image.convert("RGB"), dtype=np.uint8)
        nonwhite_fraction = float(np.any(rgb < 250, axis=2).mean())
        width, height = image.size
    if nonwhite_fraction < 0.02:
        raise AssertionError("Rendered Figure 2f is unexpectedly blank")

    summary = pd.read_csv(AUDIT / SUMMARY_OUTPUT, sep="\t")
    qa_rows = [
        ("rna_rows", len(rna), "PASS"),
        ("rna_markers", rna["marker"].nunique(), "PASS"),
        ("rna_material_stage_cells", len(rna), "PASS"),
        ("st9_cells_checked", st9_checked, "PASS"),
        ("st9_max_abs_difference", f"{st9_max_difference:.12g}", "PASS"),
        ("fad2_170d_FL", f"{fad2_fl_170:.6f}", "PASS"),
        ("fad2_170d_TN", f"{fad2_tn_170:.6f}", "PASS"),
        ("fad2_170d_FL_TN_ratio", f"{fad2_ratio:.9f}", "PASS"),
        ("fad2_protein_detected_FL", stats["fad2_fl"], "PASS"),
        ("fad2_protein_detected_TN", stats["fad2_tn"], "PASS"),
        ("summary_rows", len(summary), "PASS"),
        ("png_dimensions", f"{width}x{height}", "PASS"),
        ("png_nonwhite_fraction", f"{nonwhite_fraction:.6f}", "PASS"),
        ("pdf_bytes", pdf.stat().st_size, "PASS"),
        ("sample_label_mapping", "line1=FL; line4=TN; author confirmed 2026-09-02", "PASS"),
    ]
    pd.DataFrame(qa_rows, columns=["check", "observed", "status"]).to_csv(
        QA_OUTPUT, sep="\t", index=False
    )
    pd.DataFrame(
        [
            {"file": str(path), "sha256": sha256(path), "size_bytes": path.stat().st_size}
            for path in inputs
        ]
    ).to_csv(MANIFEST_OUTPUT, sep="\t", index=False)
    print(f"PASS: {pdf}")
    print(f"PASS: {png}")
    print(f"FAD2 RNA 170d FL/TN: {fad2_ratio:.9f}")
    print(f"ST9 cells checked: {st9_checked}; max abs difference: {st9_max_difference:.3g}")


if __name__ == "__main__":
    main()
