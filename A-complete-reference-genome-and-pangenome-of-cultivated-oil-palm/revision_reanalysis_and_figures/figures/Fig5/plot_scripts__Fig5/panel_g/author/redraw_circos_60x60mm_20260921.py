#!/usr/bin/env python3
"""Redraw the final SV circos panel at an exact 60 x 60 mm physical size.

Everything is scaled uniformly from the original 9.8 in design
(scale = (60/25.4)/9.8), so the 60 mm version is a faithful
proportional miniature of 9_SV_circos_finalSV_syriLR_with_dSV_20260630.pdf.
"""
from __future__ import annotations

import importlib.util
import math
import re
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import Patch

FIG_ROOT = Path("${ANALYSIS_DIR}/22_answer_reviews/00_ms/03_V3/05_figure")
OUT_DIR = FIG_ROOT / "06_Figure5_finalSV_syriLR_20260630"
FIG_DIR = OUT_DIR / "figures"
SCRIPT = OUT_DIR / "scripts/plot_final_syriLR_sv_figures_20260630.py"

spec = importlib.util.spec_from_file_location("plotmod", SCRIPT)
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)

MM = 25.4
SIZE_MM = 60.0
# scale relative to the original tight-bbox PDF page width (557.712 pt)
SCALE = (SIZE_MM / MM) / (557.712 / 72.0)

STEM = "9_SV_circos_finalSV_syriLR_with_dSV_20260630_60x60mm"


def main() -> None:
    catalog = mod.load_final_catalog()
    genome = mod.load_genome()
    dsv = mod.load_dsv_positions(genome)
    density = mod.make_density(catalog, genome, dsv)

    gap = 0.016
    total_size = float(genome["Size"].sum())
    usable = 2 * math.pi - gap * len(genome)
    offsets = {}
    current = math.pi / 2
    for row in genome.itertuples(index=False):
        span = usable * float(row.Size) / total_size
        offsets[row.Chrom] = (current, current + span, span, int(row.Size))
        current += span + gap

    def theta(chrom: str, pos0: float) -> float:
        start, _, span, size = offsets[chrom]
        return start + span * pos0 / size

    fig = plt.figure(figsize=(SIZE_MM / MM, SIZE_MM / MM), facecolor="white")
    # Use essentially the full 60 mm canvas. All chromosome labels remain
    # inside the polar limits, so no external margin is needed.
    ax = fig.add_axes([0.002, 0.002, 0.996, 0.996], projection="polar")
    ax.set_theta_direction(-1)
    ax.set_theta_offset(0)
    ax.set_ylim(0, 1.18)
    ax.axis("off")

    for chrom, (a0, a1, span, size) in offsets.items():
        center = (a0 + a1) / 2
        ax.bar(center, 0.025, width=span * 0.985, bottom=1.04, color="#DDE3EA", edgecolor="none", align="center")
        label = re.sub(r"chr0?([0-9]+)B", r"\1", chrom)
        ax.text(center, 1.115, label, ha="center", va="center",
                rotation=math.degrees(math.pi / 2 - center), rotation_mode="anchor",
                fontsize=6.5)

    # Compress the five tracks toward the outside to provide a clean central
    # legend area while retaining clear separation between tracks.
    rings = [
        ("All", 0.91, 0.100, mod.TYPE_COLORS["All"]),
        ("INS", 0.79, 0.095, mod.TYPE_COLORS["INS"]),
        ("DEL", 0.67, 0.095, mod.TYPE_COLORS["DEL"]),
        ("Minor", 0.55, 0.080, mod.TYPE_COLORS["Minor"]),
        ("dSV", 0.45, 0.070, "#C51B8A"),
    ]
    for column, base, height, color in rings:
        max_count = max(1, int(density[column].max()))
        for row in density.itertuples(index=False):
            width = theta(row.Chrom, row.End0) - theta(row.Chrom, row.Start0)
            center = theta(row.Chrom, (row.Start0 + row.End0) / 2)
            h = height * getattr(row, column) / max_count
            ax.bar(center, h, width=width * 0.96, bottom=base, color=color, edgecolor=color, linewidth=0, align="center")

    hotspot = density[density["Hotspot"]]
    for row in hotspot.itertuples(index=False):
        center = theta(row.Chrom, (row.Start0 + row.End0) / 2)
        ax.scatter(center, 1.015, s=5.0, color="#CB4B4B", edgecolor="none", zorder=5)

    handles = [
        Patch(facecolor=mod.TYPE_COLORS["All"], label="All SVs"),
        Patch(facecolor=mod.TYPE_COLORS["INS"], label="INS"),
        Patch(facecolor=mod.TYPE_COLORS["DEL"], label="DEL"),
        Patch(facecolor=mod.TYPE_COLORS["Minor"], label="INV+DUP+TRA"),
        Patch(facecolor="#CB4B4B", label="Hotspot"),
        Patch(facecolor="#C51B8A", label="dSV density"),
    ]
    legend = ax.legend(
        handles=handles, title="Final SV density",
        loc="center", bbox_to_anchor=(0.5, 0.50),
        ncol=1, frameon=False, fontsize=6.0, title_fontsize=7.0,
        handlelength=0.9, handleheight=0.7, handletextpad=0.4,
        labelspacing=0.12, borderaxespad=0,
    )
    legend.get_title().set_fontweight("bold")

    # Keep the requested physical canvas size; tight cropping would alter it.
    for ext, kwargs in (
        ("pdf", dict(facecolor="white")),
        ("svg", dict(facecolor="white")),
        ("png", dict(dpi=600, facecolor="white")),
    ):
        fig.savefig(FIG_DIR / f"{STEM}.{ext}", **kwargs)
    plt.close(fig)
    print(f"[DONE] {FIG_DIR / (STEM + '.pdf')}")


if __name__ == "__main__":
    main()
