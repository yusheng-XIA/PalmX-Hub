#!/usr/bin/env python3
"""Redraw Figure 4d at its final size with clearer haplotype line styles.

Scientific values are read verbatim from the reviewed Fig4d_Nx_pan39.tsv table.
No Nx statistic is recalculated. The only intended changes are final-size
layout and stronger visual encoding of hap1 versus hap2.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd


INPUT = Path(
    "${ANALYSIS_DIR}/"
    "22_answer_reviews/00_ms/03_V3/04_figure4/"
    "Fig4_d_i_pan39_material33_singletons_20260811/tables/Fig4d_Nx_pan39.tsv"
)
OUT_DIR = Path(
    "${ANALYSIS_DIR}/"
    "22_answer_reviews/00_ms/05_MS/0918_revision/final_ms/02_figure4"
)
OUTPUT = OUT_DIR / "Figure4d_60x38.5mm_vector_solid_dashed_v3.pdf"
QC_PNG = OUT_DIR / "qc/Figure4d_60x38.5mm_solid_dashed_v3_600dpi.png"
AUDIT = OUT_DIR / "provenance/Figure4d_solid_dashed_v3_audit.json"
EXPECTED_INPUT_SHA256 = "698b97e3bcdcb7619ef4192cdeb9e83e909019733323ae562bb87730364deba1"

NX_COLS = [f"N{x}" for x in range(10, 101, 10)]
X = np.arange(10, 101, 10, dtype=float)

# Mapping is retained from the reviewed panel and its final Illustrator legend.
PAIR_META = {
    "American_hap1": ("FL", "hap1", "#E65D6D"),
    "Africa_hap2": ("FL", "hap2", "#E65D6D"),
    "bk_hap1": ("TN", "hap1", "#BBA8D8"),
    "bk_hap2": ("TN", "hap2", "#BBA8D8"),
    "dura_hap1": ("Dura", "hap1", "#4667AE"),
    "dura_hap2": ("Dura", "hap2", "#4667AE"),
    "pisifera_hap1": ("Pisifera", "hap1", "#F2B76A"),
    "pisifera_hap2": ("Pisifera", "hap2", "#F2B76A"),
    "nrly_hap1": ("Nigerian", "hap1", "#29AFD4"),
    "nrly_hap2": ("Nigerian", "hap2", "#29AFD4"),
    "meizhou4_hap1": ("Oleifera", "hap1", "#A7D9DD"),
    "meizhou4_hap2": ("Oleifera", "hap2", "#A7D9DD"),
}
MATERIAL_ORDER = ["FL", "TN", "Dura", "Pisifera", "Nigerian", "Oleifera"]
MATERIAL_COLORS = {
    label: color for label, _hap, color in PAIR_META.values()
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_reviewed_values() -> tuple[pd.DataFrame, np.ndarray]:
    observed_sha = sha256(INPUT)
    if observed_sha != EXPECTED_INPUT_SHA256:
        raise RuntimeError(
            f"Input checksum changed: {observed_sha}; expected {EXPECTED_INPUT_SHA256}"
        )
    frame = pd.read_csv(INPUT, sep="\t")
    if len(frame) != 39 or frame["Sample"].nunique() != 39:
        raise RuntimeError("Reviewed table must contain exactly 39 unique assemblies")
    if set(PAIR_META) - set(frame["Sample"]):
        raise RuntimeError("Reviewed table is missing one or more highlighted haplotypes")
    numeric = frame[NX_COLS].apply(pd.to_numeric, errors="raise").to_numpy(dtype=float)
    if not np.all(np.isfinite(numeric)) or np.any(numeric <= 0):
        raise RuntimeError("Nx values must be finite positive numbers")
    return frame, numeric / 1e6


def configure_style() -> None:
    mpl.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Liberation Sans", "DejaVu Sans"],
            "font.size": 6.0,
            "axes.labelsize": 6.5,
            "axes.linewidth": 0.55,
            "xtick.labelsize": 5.5,
            "ytick.labelsize": 5.5,
            "xtick.major.width": 0.45,
            "ytick.major.width": 0.45,
            "xtick.major.size": 2.2,
            "ytick.major.size": 2.2,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
        }
    )


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    QC_PNG.parent.mkdir(parents=True, exist_ok=True)
    AUDIT.parent.mkdir(parents=True, exist_ok=True)
    frame, values_mb = load_reviewed_values()
    configure_style()

    mm = 1.0 / 25.4
    fig, ax = plt.subplots(figsize=(60.0 * mm, 38.5 * mm), facecolor="white")
    fig.subplots_adjust(left=0.145, right=0.992, bottom=0.185, top=0.975)

    samples = frame["Sample"].astype(str).tolist()
    highlighted = set(PAIR_META)
    plotted_values: dict[str, list[float]] = {}

    # Preserve all 27 non-highlighted assembly curves as subdued context.
    for row_index, sample in enumerate(samples):
        if sample in highlighted:
            continue
        y = values_mb[row_index]
        ax.plot(X, y, color="#B3B3B3", lw=0.28, alpha=0.34, zorder=1)
        plotted_values[sample] = y.tolist()

    # Preserve the across-assembly median, but use a faint dotted style so it
    # cannot be mistaken for the long-dashed hap2 encoding.
    median = np.median(values_mb, axis=0)
    ax.plot(X, median, color="#626262", lw=0.62, ls=(0, (1.0, 1.6)), alpha=0.82, zorder=2)

    # Colour identifies material; line type alone identifies haplotype.
    # Draw hap2 last so its dashed pattern remains visible above overlaps.
    style = {
        "hap1": "-",
        "hap2": (0, (4.0, 1.6)),
    }
    for desired_haplotype in ("hap1", "hap2"):
        for row_index, sample in enumerate(samples):
            if sample not in PAIR_META:
                continue
            _material, haplotype, color = PAIR_META[sample]
            if haplotype != desired_haplotype:
                continue
            y = values_mb[row_index]
            ax.plot(
                X,
                y,
                color=color,
                lw=1.25,
                ls=style[haplotype],
                alpha=1.0,
                solid_capstyle="round",
                dash_capstyle="butt",
                zorder=4 if haplotype == "hap1" else 5,
            )
            plotted_values[sample] = y.tolist()

    # Assert that every plotted assembly curve equals the reviewed table.
    if set(plotted_values) != set(samples):
        raise RuntimeError("Not every reviewed assembly was plotted exactly once")
    for row_index, sample in enumerate(samples):
        if not np.array_equal(np.asarray(plotted_values[sample]), values_mb[row_index]):
            raise RuntimeError(f"Plotted values differ from reviewed table for {sample}")

    hap_handles = [
        Line2D([0], [0], color="#3F3F3F", lw=1.40, ls="-", label="hap1"),
        Line2D([0], [0], color="#3F3F3F", lw=1.40, ls=(0, (4.0, 1.6)), label="hap2"),
    ]
    material_handles = [
        Line2D([0], [0], color=MATERIAL_COLORS[name], lw=1.15, label=name)
        for name in MATERIAL_ORDER
    ]
    hap_legend = ax.legend(
        handles=hap_handles,
        loc="upper right",
        bbox_to_anchor=(0.69, 0.995),
        frameon=False,
        fontsize=5.5,
        handlelength=3.8,
        handletextpad=0.45,
        borderaxespad=0.0,
        labelspacing=0.35,
    )
    ax.add_artist(hap_legend)
    ax.legend(
        handles=material_handles,
        loc="upper right",
        bbox_to_anchor=(1.0, 0.995),
        frameon=False,
        fontsize=5.35,
        handlelength=2.5,
        handletextpad=0.45,
        borderaxespad=0.0,
        labelspacing=0.25,
    )

    ax.set_xlim(8, 102)
    ax.set_xticks(X)
    ax.set_ylim(bottom=0)
    ax.set_xlabel(r"$N_x$ (%)", labelpad=1.5)
    ax.set_ylabel("Contig length (Mb)", labelpad=1.5)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_color("#666666")
    ax.spines["bottom"].set_color("#666666")
    ax.tick_params(colors="#202020", pad=1.4)
    fig.text(0.008, 0.988, "d", ha="left", va="top", fontsize=10.0, fontweight="bold")

    # Do not use bbox_inches='tight': preserve the exact physical page size.
    fig.savefig(OUTPUT, facecolor="white")
    fig.savefig(QC_PNG, dpi=600, facecolor="white")
    plt.close(fig)

    audit = {
        "status": "PASS",
        "operation": "visual-style redraw from reviewed terminal table",
        "scientific_values_recalculated": False,
        "input": str(INPUT),
        "input_sha256": sha256(INPUT),
        "input_rows": len(frame),
        "input_unique_samples": int(frame["Sample"].nunique()),
        "nx_columns": NX_COLS,
        "highlighted_haplotypes": list(PAIR_META),
        "background_assemblies": len(frame) - len(PAIR_META),
        "hap1_style": "1.25 pt solid line; no markers",
        "hap2_style": "1.25 pt long-dashed line (4.0 pt on, 1.6 pt off); no markers; drawn above hap1",
        "output_page_mm": [60.0, 38.5],
        "output": str(OUTPUT),
        "output_sha256": sha256(OUTPUT),
        "all_plotted_values_exactly_equal_input": True,
    }
    AUDIT.write_text(json.dumps(audit, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(audit, indent=2))


if __name__ == "__main__":
    main()
