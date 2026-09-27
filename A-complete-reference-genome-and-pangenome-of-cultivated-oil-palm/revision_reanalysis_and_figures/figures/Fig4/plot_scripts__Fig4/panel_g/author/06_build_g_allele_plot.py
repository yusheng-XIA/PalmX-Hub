#!/usr/bin/env python3
"""Combine accepted and newly recomputed allele summaries and draw Fig. 4g."""

from __future__ import annotations

import csv
import json
import os
from pathlib import Path

os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MPLBACKEND", "Agg")
os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib-fig4g-pan39")
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import numpy as np


RUN = Path("${ANALYSIS_DIR}/22_answer_reviews/00_ms/03_V3/04_figure4/Fig4_d_i_pan39_material33_redraw_20260808")
OLD = Path("${ANALYSIS_DIR}/10_genome_ann_contigs/13_allele/06_final_step/full_summary.tsv")
ORDER_TABLE = Path("${ANALYSIS_DIR}/22_answer_reviews/00_ms/03_V3/04_figure4/RGA_BGC_meizhou4_33varieties_20260805/tables/rga_plot/RGA_major_group_counts_and_zscores.tsv")
NEW_DIR = RUN / "tables/g_new_pairs"
CATEGORIES = ("Bi-allelic genes", "Haplotype-specific genes", "Allele with same CDS", "Unresolved RBH pairs")
COLORS = {"Bi-allelic genes": "#8FDAEF", "Haplotype-specific genes": "#E65D6D",
          "Allele with same CDS": "#F8AFA1", "Unresolved RBH pairs": "#728AC1"}
DISPLAY_NAMES = {
    "seedless": "FL (seedless)",
    "bk": "TN (boke)",
    "dura": "Dura",
    "pisifera": "Pisifera",
    "nrly": "NRLY",
    "houke_old": "Houke",
    "meizhou4": "Meizhou4",
}


def material_from_old(sample: str) -> str | None:
    if sample.endswith("_pa_genome") and sample.split("_", 1)[0].isdigit():
        return sample.split("_", 1)[0]
    return {"American_Africa": "seedless", "bk_hap1_hap2": "bk",
            "houke_pa_genome": "houke_old"}.get(sample)


def categories(row: dict[str, str]) -> dict[str, int]:
    # High-confidence haplotype-specific genes are the intersection supported by
    # both CDS-BLAST and GMAP.  The two marginal counts are method-level QC
    # quantities and must not be added to their intersection.
    return {"Bi-allelic genes": int(row["Heterozygous"]),
            "Haplotype-specific genes": int(row["Unique_inter"]),
            "Allele with same CDS": int(row["Homozygous"]),
            "Unresolved RBH pairs": int(row["No_hit"])}


def display_name(material: str) -> str:
    return DISPLAY_NAMES.get(material, material)


def load_rows() -> dict[str, dict[str, object]]:
    result: dict[str, dict[str, object]] = {}
    with OLD.open() as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            material = material_from_old(row["Sample"])
            if material is None:
                continue
            result[material] = {"Material": material, "Source": str(OLD), **categories(row)}
    new_map = {"dura_new": "dura", "pisifera_new": "pisifera", "nrly_new": "nrly", "meizhou4_new": "meizhou4"}
    for pair, material in new_map.items():
        # attempt7_uniqueid corrects the shared hap1/hap2 gene-ID namespace
        # that caused GeneTribe to discard same-name allele pairs.
        path = NEW_DIR / f"{pair}.allele_summary.attempt7_uniqueid.tsv"
        if not path.is_file() or path.stat().st_size == 0:
            raise FileNotFoundError(path)
        with path.open() as handle:
            row = next(csv.DictReader(handle, delimiter="\t"))
        result[material] = {"Material": material, "Source": str(path), **categories(row)}
    if len(result) != 33:
        raise AssertionError(f"Expected 33 material summaries, observed {len(result)}: {sorted(result)}")
    return result


def tree_order() -> list[str]:
    tips = []
    with ORDER_TABLE.open() as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            tip = row["Tree_tip"]
            if tip not in tips:
                tips.append(tip)
    mapping = {"FL": "seedless", "BK": "bk", "dura": "dura", "pisifera": "pisifera",
               "niriliya_pa_genome": "nrly", "houke_pa_genome": "houke_old",
               "meizhou4_pa_genome": "meizhou4"}
    order = []
    for tip in tips:
        if tip in mapping:
            material = mapping[tip]
        elif tip.endswith("_pa_genome") and tip.split("_", 1)[0].isdigit():
            material = tip.split("_", 1)[0]
        else:
            continue
        if material not in order:
            order.append(material)
    if len(order) != 33:
        raise AssertionError(f"Display tree order has {len(order)} materials: {order}")
    return order


def write_source(rows: dict[str, dict[str, object]], order: list[str]) -> None:
    path = RUN / "tables/Fig4g_allele_composition_material33.tsv"
    with path.open("w", newline="") as handle:
        fields = ["Display_order", "Material", "Display_name", *CATEGORIES, "Total", "Unresolved_RBH_rate_pct", "Source"]
        writer = csv.DictWriter(handle, fields, delimiter="\t")
        writer.writeheader()
        for i, material in enumerate(order, 1):
            row = rows[material]
            total = sum(int(row[c]) for c in CATEGORIES)
            unresolved = int(row["Unresolved RBH pairs"])
            writer.writerow({"Display_order": i, "Material": material, "Display_name": display_name(material),
                             **{c: row[c] for c in CATEGORIES}, "Total": total,
                             "Unresolved_RBH_rate_pct": f"{unresolved / total * 100:.4f}",
                             "Source": row["Source"]})


def plot(rows: dict[str, dict[str, object]], order: list[str], labelled: bool) -> None:
    plt.rcParams.update({"font.family": "Liberation Sans", "font.size": 9, "axes.labelsize": 11,
                         "axes.linewidth": 0.7, "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none"})
    fig, ax = plt.subplots(figsize=(9.8, 6.7), facecolor="white")
    y = np.arange(len(order))
    left = np.zeros(len(order), dtype=float)
    for category in CATEGORIES:
        values = []
        for material in order:
            row = rows[material]
            total = sum(int(row[c]) for c in CATEGORIES)
            values.append(int(row[category]) / total * 100.0)
        values = np.asarray(values)
        ax.barh(y, values, left=left, height=0.80, color=COLORS[category], edgecolor="white", linewidth=0.35)
        left += values
    ax.set_xlim(0, 100)
    ax.set_xticks([0, 25, 50, 75, 100])
    ax.set_xlabel("Percentage (%)")
    ax.set_yticks(y)
    ax.set_yticklabels([display_name(material) for material in order], fontsize=7.3)
    ax.invert_yaxis()
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.spines["bottom"].set_color("#808080")
    ax.tick_params(axis="y", length=0)
    handles = [Patch(facecolor=COLORS[c], edgecolor="none", label=c) for c in CATEGORIES]
    ax.legend(handles=handles, frameon=False, ncol=2, loc="lower left", bbox_to_anchor=(0, 1.005),
              borderaxespad=0, handlelength=0.9, handletextpad=0.4, columnspacing=1.2, fontsize=8)
    fig.subplots_adjust(left=0.14, right=0.985, bottom=0.11, top=0.88)
    if labelled:
        fig.text(0.012, 0.985, "g", ha="left", va="top", fontsize=17, fontweight="bold")
    stem = f"Fig4g_allele_composition_material33_{'labelled' if labelled else 'no_label'}"
    fig.savefig(RUN / "figures" / f"{stem}.pdf", facecolor="white")
    fig.savefig(RUN / "figures" / f"{stem}.svg", facecolor="white")
    fig.savefig(RUN / "figures" / f"{stem}_600dpi.png", dpi=600, facecolor="white")
    plt.close(fig)


def main() -> None:
    rows = load_rows()
    order = tree_order()
    if set(rows) != set(order):
        raise AssertionError(f"Summary/order mismatch: summaries-only={set(rows)-set(order)}, order-only={set(order)-set(rows)}")
    write_source(rows, order)
    for labelled in (False, True):
        plot(rows, order, labelled)
    unresolved = {}
    for material in order:
        total = sum(int(rows[material][c]) for c in CATEGORIES)
        count = int(rows[material]["Unresolved RBH pairs"])
        unresolved[material] = {"display_name": display_name(material), "count": count,
                                "rate_pct": round(count / total * 100, 4)}
    audit = {
        "materials": len(rows),
        "new_recomputed_pairs": 4,
        "new_pair_attempt": "attempt7_uniqueid",
        "retained_accepted_summaries": 29,
        "display_order_source": str(ORDER_TABLE),
        "category_definitions": {
            "Bi-allelic genes": "GeneTribe RBH pairs with a non-identical full-gene BLASTN match",
            "Haplotype-specific genes": "a-haplotype genes absent by both CDS-BLAST and GMAP (Unique_inter)",
            "Allele with same CDS": "GeneTribe RBH pairs with 100% identity across the complete query gene",
            "Unresolved RBH pairs": "GeneTribe RBH pairs without the exact corresponding a-gene to p-gene BLASTN pair; formerly labelled Others",
        },
        "method_counting_correction": "Unique_inter alone is plotted as the high-confidence haplotype-specific count; CDS_unique and GMAP_unique are method-level marginal counts and are not added to their intersection.",
        "unresolved_RBH_by_material": unresolved,
        "status": "PASS",
    }
    (RUN / "provenance/Fig4g_audit.json").write_text(json.dumps(audit, indent=2) + "\n")
    (RUN / "provenance/Fig4g.SUCCESS").touch()
    print(json.dumps(audit, indent=2))


if __name__ == "__main__":
    main()
