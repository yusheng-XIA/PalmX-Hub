#!/usr/bin/env python3
"""Draw the original Figure 4 RGA panel with 33 biological varieties."""

from __future__ import annotations

import importlib.util
from pathlib import Path


RUN = Path(
    "${ANALYSIS_DIR}/22_answer_reviews/00_ms/03_V3/04_figure4/"
    "RGA_BGC_meizhou4_33varieties_20260805"
)
SOURCE = Path(
    "${ANALYSIS_DIR}/22_answer_reviews/00_ms/03_V3/04_figure4/"
    "RGA_tree_dotmatrix_redesign_20260708/scripts/build_rga_tree_dotmatrix_redesign.py"
)


def main() -> None:
    spec = importlib.util.spec_from_file_location("rga_original", SOURCE)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot import {SOURCE}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)

    module.RGA_DIR = RUN / "tables"
    module.TREE_FILE = RUN / "results/SpeciesTree_33varieties_display.nwk"
    module.OUT_ROOT = RUN
    module.FIG_DIR = RUN / "figures"
    module.TAB_DIR = RUN / "tables/rga_plot"
    module.SYNTENY_COUNTS_FILE = RUN / "tables/RGA_synteny_counts_33varieties.tsv"

    def load_rga_summary():
        df = module.pd.read_csv(RUN / "tables/RGA_summary_33varieties.csv")
        df["NBS-class"] = df[["NBS", "CNL", "TNL", "CN", "TN", "NL", "TX"]].sum(axis=1)
        df["Other"] = df["Others"] if "Others" in df.columns else 0
        return df[["Genome"] + module.RGA_COLUMNS].copy()

    module.load_rga_summary = load_rga_summary
    original_display = module.display_genome_name
    original_group = module.genome_group

    def display_genome_name(name: str) -> str:
        if name == "meizhou4":
            return "Meizhou4"
        return original_display(name)

    def genome_group(name: str) -> str:
        if name == "FL":
            return "African/American haplotypes"
        if name == "BK":
            return "BK haplotypes"
        return original_group(name)

    module.display_genome_name = display_genome_name
    module.genome_group = genome_group
    module.main()


if __name__ == "__main__":
    main()
