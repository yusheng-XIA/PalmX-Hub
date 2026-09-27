#!/usr/bin/env python3
"""Render current Stage 09 load-only IPHs with legacy GWAS marks projected.

Visual source copied from frozen legacy script:
results/09_ideal_parent_haplotypes/scripts/16_plot_composite_step9_v5.py

The current All38/African35 paths remain the accepted GWAS_Weight=0 paths.
GWAS coordinates and gold/grey states are copied only as a clearly labelled
legacy W=4 status projection; they are not recomputed current-path genotypes.
"""

import argparse
import csv
import hashlib
import os
import re
from collections import Counter
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.font_manager as fm
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D
from matplotlib.patches import FancyBboxPatch, Patch, Rectangle
import numpy as np
import pandas as pd
from PIL import Image


WINDOW = 500_000
ZC = "chr01B"
SOURCE_SCRIPT = Path(
    "${ANALYSIS_DIR}/21_MS/06_result/dSVs/"
    "results/09_ideal_parent_haplotypes/scripts/16_plot_composite_step9_v5.py"
)

# Copied verbatim from legacy v5.
MOSAIC_COLORS = [
    "#E9F3FA", "#D8E8F8", "#C8D8E8", "#AFCFE3",
    "#8FBAD6", "#68A8D8", "#5A94BF", "#7898C8",
    "#627FB0", "#4D6F9A", "#3D638E", "#2D536F", "#213F58",
]
HEATMAP_COLORS = [
    "#FBFCFE", "#EFF6FB", "#E2EFF7", "#CDE2F1",
    "#AFCFE6", "#88BBDC", "#68A8D8", "#3F83B0", "#245E7A",
]
CHROM_BG = "#F3F7FA"
CAPTURED_COLOR = "#D9A441"
MISSED_COLOR = "#B8C5CC"
SELECT_RED = "#C84848"
SELECT_TEXT_RED = "#B83F3E"


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--reference-fai", required=True)
    parser.add_argument("--legacy-loci", required=True)
    parser.add_argument("--legacy-by-chrom", required=True)
    parser.add_argument("--panels", default="all38,african35")
    return parser.parse_args()


def chrom_sort_key(chrom):
    match = re.fullmatch(r"chr(\d+)B", chrom)
    return int(match.group(1)) if match else 10**9


def setup_style():
    for candidate in ("Arial", "Helvetica", "DejaVu Sans"):
        if any(candidate == font.name for font in fm.fontManager.ttflist):
            plt.rcParams["font.family"] = candidate
            break
    plt.rcParams.update(
        {
            "axes.edgecolor": "#9aa3ab",
            "axes.linewidth": 0.6,
            "xtick.color": "#4d555c",
            "ytick.color": "#4d555c",
            "text.color": "#2b3136",
            "axes.labelcolor": "#2b3136",
            "xtick.major.width": 0.6,
            "mathtext.fontset": "dejavusans",
            "mathtext.default": "regular",
        }
    )
    return plt.rcParams["font.family"]


def read_fai(path):
    lengths = {}
    with open(path) as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 2:
                lengths[fields[0]] = int(fields[1])
    expected = {f"chr{i:02d}B" for i in range(1, 17)}
    if set(lengths) != expected:
        raise ValueError("Reference chromosome set mismatch")
    chroms = sorted(lengths, key=chrom_sort_key)
    return lengths, chroms


def read_legacy_gwas(loci_path, by_chrom_path, chroms, chromlen):
    loci_mark = {}
    row_states = Counter()
    row_count = 0
    with open(loci_path, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"SV", "Chrom", "Pos", "Captured"}
        if not required.issubset(reader.fieldnames or []):
            raise ValueError(f"Legacy loci columns missing: {required}")
        for row in reader:
            chrom = row["Chrom"]
            pos = int(row["Pos"])
            captured = int(row["Captured"])
            if chrom not in chromlen or not 0 <= pos <= chromlen[chrom]:
                raise ValueError(f"Legacy GWAS coordinate out of range: {chrom}:{pos}")
            if captured not in (0, 1):
                raise ValueError(f"Invalid legacy Captured value: {captured}")
            key = (chrom, pos)
            loci_mark[key] = max(loci_mark.get(key, 0), captured)
            row_states[captured] += 1
            row_count += 1

    by_chrom = pd.read_csv(by_chrom_path, sep="\t")
    by_chrom = by_chrom.loc[by_chrom["Chrom"] != "TOTAL"].copy()
    legacy_counts = {
        row.Chrom: (int(row.Fav_Captured), int(row.Fav_Total))
        for row in by_chrom.itertuples(index=False)
    }
    if set(legacy_counts) != set(chroms):
        raise ValueError("Legacy GWAS chromosome-count set mismatch")
    marker_states = Counter(loci_mark.values())
    summary = {
        "Legacy_Locus_Rows": row_count,
        "GWAS_Position_Marks": len(loci_mark),
        "Legacy_Row_Captured": row_states[1],
        "Legacy_Row_Missed": row_states[0],
        "Projected_Gold_Marks": marker_states[1],
        "Projected_Grey_Marks": marker_states[0],
    }
    expected = {
        "Legacy_Locus_Rows": 284,
        "GWAS_Position_Marks": 282,
        "Legacy_Row_Captured": 184,
        "Legacy_Row_Missed": 100,
        "Projected_Gold_Marks": 183,
        "Projected_Grey_Marks": 99,
    }
    if summary != expected:
        raise ValueError(f"Legacy GWAS input fingerprint mismatch: {summary}")
    return loci_mark, legacy_counts, summary


def read_panel_data(outroot, panel_key, chroms, chromlen):
    panel_dir = outroot / panel_key
    matrix = pd.read_csv(panel_dir / "load_matrix_500kb.tsv", sep="\t")
    path_df = pd.read_csv(panel_dir / "ideal_loadonly_path.tsv", sep="\t")
    by_chrom = pd.read_csv(panel_dir / "ideal_loadonly_by_chrom.tsv", sep="\t")
    by_chrom = by_chrom.loc[by_chrom["Chrom"] != "TOTAL"].copy()
    path = {chrom: {} for chrom in chroms}
    for row in path_df.itertuples(index=False):
        path[row.Chrom][int(row.Window_Index)] = row.Donor_ID
    for chrom in chroms:
        expected_windows = chromlen[chrom] // WINDOW + 1
        observed_windows = len(path[chrom])
        if expected_windows != observed_windows:
            raise ValueError(
                f"Legacy-v5 window count mismatch for {panel_key} {chrom}: "
                f"expected={expected_windows} observed={observed_windows}"
            )
    counts = {
        row.Chrom: (int(row.Residual_DSV), int(row.Residual_DSNP))
        for row in by_chrom.itertuples(index=False)
    }
    return panel_dir, matrix, path, counts


def contiguous_segments(path_for_chrom, n_win):
    segments = []
    current = None
    start = 0
    for window in range(n_win):
        donor = path_for_chrom.get(window)
        if donor != current:
            if current is not None:
                segments.append((start, window, current))
            current = donor
            start = window
    if current is not None:
        segments.append((start, n_win, current))
    return segments


def render_panel(
    outroot, panel_key, chromlen, chroms, source_sha256,
    loci_mark, legacy_counts, gwas_summary, loci_sha256,
):
    display = "All38" if panel_key == "all38" else "African35"
    panel_dir, matrix, path, counts = read_panel_data(outroot, panel_key, chroms, chromlen)
    figure_dir = panel_dir / "figures_legacy_v5_gwas_projection"
    figure_dir.mkdir(parents=True, exist_ok=False)

    contribution = {}
    for chrom in chroms:
        for donor in path[chrom].values():
            contribution[donor] = contribution.get(donor, 0) + 1
    used = sorted(contribution, key=lambda donor: -contribution[donor])
    donor_cmap = LinearSegmentedColormap.from_list("sv3_mosaic_v5_blue", MOSAIC_COLORS)
    lo, span = 0.06, 0.88
    shade = {
        donor: donor_cmap(lo + span * (index / max(1, len(used) - 1)))
        for index, donor in enumerate(used)
    }

    # Legacy v5 chr01B donor ordering: ascending total load on the zoom chromosome.
    zoom = matrix.loc[matrix["Chrom"] == ZC].copy()
    donors = sorted(zoom["Sample_ID"].unique())
    n_win = chromlen[ZC] // WINDOW + 1
    total_pivot = zoom.pivot(index="Sample_ID", columns="Window_Index", values="Total_Load")
    dsv_pivot = zoom.pivot(index="Sample_ID", columns="Window_Index", values="DSV_Count")
    dsnp_pivot = zoom.pivot(index="Sample_ID", columns="Window_Index", values="DSNP_Count")
    columns = list(range(n_win))
    load = total_pivot.reindex(index=donors, columns=columns, fill_value=0).to_numpy(dtype=int)
    dsv_load = dsv_pivot.reindex(index=donors, columns=columns, fill_value=0).to_numpy(dtype=int)
    dsnp_load = dsnp_pivot.reindex(index=donors, columns=columns, fill_value=0).to_numpy(dtype=int)
    order = np.argsort(load.sum(1))
    donors = [donors[index] for index in order]
    load = load[order]
    dsv_load = dsv_load[order]
    dsnp_load = dsnp_load[order]
    row_of = {donor: index for index, donor in enumerate(donors)}

    # ===== Drawing block follows legacy v5 geometry and style. =====
    fig = plt.figure(figsize=(20, 8.8))
    fig.patch.set_facecolor("white")
    gs = GridSpec(
        2, 2, height_ratios=[8, 1.2], width_ratios=[1.02, 1.0],
        hspace=0.06, wspace=0.14, left=0.05, right=0.985, top=0.9, bottom=0.08,
    )
    axm = fig.add_subplot(gs[0, 0])
    axl = fig.add_subplot(gs[1, 0])
    axz = fig.add_subplot(gs[0, 1])
    ymax = len(chroms)
    bar_height = 0.6
    max_mb = max(chromlen.values()) / 1e6

    total_segment_count = 0
    for yi, chrom in enumerate(chroms):
        y = ymax - yi
        nw = chromlen[chrom] // WINDOW + 1
        axm.add_patch(
            FancyBboxPatch(
                (0, y - bar_height / 2), chromlen[chrom] / 1e6, bar_height,
                boxstyle="round,pad=0,rounding_size=0.2", linewidth=0,
                facecolor=CHROM_BG, zorder=1,
            )
        )
        segments = contiguous_segments(path[chrom], nw)
        total_segment_count += len(segments)
        for start, end, donor in segments:
            axm.barh(
                y, (end - start) * WINDOW / 1e6, left=start * WINDOW / 1e6,
                height=bar_height, color=shade[donor], edgecolor="white",
                linewidth=0.35, zorder=2,
            )
        for (marker_chrom, marker_pos), captured in loci_mark.items():
            if marker_chrom != chrom:
                continue
            axm.plot(
                marker_pos / 1e6, y + bar_height / 2 + 0.16,
                marker="v", markersize=4.3,
                color=CAPTURED_COLOR if captured else MISSED_COLOR,
                markeredgecolor="white", markeredgewidth=0.3, zorder=4,
            )
        dsv, dsnp = counts.get(chrom, (0, 0))
        fav_captured, fav_total = legacy_counts[chrom]
        axm.text(
            chromlen[chrom] / 1e6 + max_mb * 0.012, y,
            f"{dsv} / {dsnp}   $\\star$ {fav_captured}/{fav_total}",
            va="center", ha="left", fontsize=7.8, color="#2b3136",
        )
        axm.text(
            -max_mb * 0.012, y, chrom, va="center", ha="right",
            fontsize=8.5, color="#57616a",
        )
    axm.set_xlim(-max_mb * 0.05, max_mb * 1.17)
    axm.set_ylim(0.3, ymax + 1.1)
    axm.set_yticks([])
    for spine in ("top", "right", "left"):
        axm.spines[spine].set_visible(False)
    axm.spines["bottom"].set_bounds(0, max_mb)
    axm.set_xticks(np.arange(0, max_mb + 1, 25))
    axm.set_xlabel("Chromosomal position (Mb)", fontsize=9)
    axm.text(
        max_mb + max_mb * 0.012, ymax + 0.62,
        r"dSVs / dSNPs   $\star$ legacy captured/total",
        fontsize=7.8, color="#57616a", ha="left", style="italic",
    )
    axm.set_title(
        f"a   Ideal parental haplotype (IPHs): load-only path + legacy GWAS projection — {display}",
        fontsize=12.5, loc="left", pad=12, color="#1c2126", fontweight="bold",
    )

    axl.axis("off")
    handles = [
        Patch(facecolor=shade[donor], edgecolor="white", linewidth=0.3, label=donor)
        for donor in used
    ]
    legend = axl.legend(
        handles=handles, loc="center left", bbox_to_anchor=(0.0, 0.5),
        ncol=7, fontsize=6.4, frameon=False, handlelength=1.1, handleheight=1.1,
        labelspacing=0.3, columnspacing=0.9,
        title=f"Donor haplotype (n={len(used)})", title_fontsize=7.2,
    )
    legend.get_title().set_color("#57616a")
    axl.add_artist(legend)
    marker_handles = [
        Line2D(
            [0], [0], marker="v", color="white", markerfacecolor=CAPTURED_COLOR,
            markersize=6, label="captured in legacy W=4 path",
        ),
        Line2D(
            [0], [0], marker="v", color="white", markerfacecolor=MISSED_COLOR,
            markersize=6, label="missed in legacy W=4 path",
        ),
    ]
    marker_legend = axl.legend(
        handles=marker_handles, loc="center right", bbox_to_anchor=(1.0, 0.5),
        fontsize=6.8, frameon=False, title="Legacy GWAS status projection",
        title_fontsize=7.0,
    )
    marker_legend.get_title().set_color("#57616a")

    heatmap = LinearSegmentedColormap.from_list("legacy_v5_heatmap", HEATMAP_COLORS)
    image = axz.imshow(
        np.log1p(load), aspect="auto", cmap=heatmap, interpolation="nearest",
        extent=[0, n_win * WINDOW / 1e6, len(donors) - 0.5, -0.5],
    )
    axz.set_yticks(range(len(donors)))
    axz.set_yticklabels(donors, fontsize=6.2, color="#57616a")
    axz.set_xlabel(f"{ZC} position (Mb)", fontsize=9)
    axz.set_title(
        f"     {ZC}: candidate donor haplotypes (load per {WINDOW // 1000}kb) "
        "— red = selected into IPHs; label = dSV/dSNP",
        fontsize=10, loc="left", pad=12, color="#1c2126",
    )
    zoom_segments = contiguous_segments(path[ZC], n_win)
    for start, end, donor in zoom_segments:
        row = row_of[donor]
        x0 = start * WINDOW / 1e6
        width = (end - start) * WINDOW / 1e6
        axz.add_patch(
            Rectangle(
                (x0, row - 0.5), width, 1, fill=False,
                edgecolor=SELECT_RED, linewidth=1.3,
            )
        )
        donor_row = row_of[donor]
        dsv = int(dsv_load[donor_row, start:end].sum())
        dsnp = int(dsnp_load[donor_row, start:end].sum())
        axz.text(
            x0 + width / 2, row - 0.62, f"{dsv}/{dsnp}",
            color=SELECT_TEXT_RED, fontsize=4.8, ha="center", va="bottom",
        )
    for (marker_chrom, marker_pos), captured in loci_mark.items():
        if marker_chrom != ZC:
            continue
        axz.plot(
            marker_pos / 1e6, -0.9, marker="v", markersize=4.5,
            clip_on=False, color=CAPTURED_COLOR if captured else MISSED_COLOR,
            markeredgecolor="white", markeredgewidth=0.3,
        )
    for spine in ("top", "right"): 
        axz.spines[spine].set_visible(False)
    colorbar = fig.colorbar(image, ax=axz, fraction=0.018, pad=0.01)
    colorbar.set_label("log(1+load / window)", fontsize=7.5)
    colorbar.ax.tick_params(labelsize=6.5)

    basename = f"Fig5a_composite_{display}_loadonly_legacy_v5_gwas_projection"
    png_path = figure_dir / f"{basename}.png"
    pdf_path = figure_dir / f"{basename}.pdf"
    if png_path.exists() or pdf_path.exists():
        raise FileExistsError(f"Refusing to overwrite legacy GWAS projection: {basename}")
    plt.savefig(png_path, dpi=270, facecolor="white", bbox_inches="tight")
    Image.open(png_path).convert("RGB").save(pdf_path, "PDF", resolution=300.0)
    plt.close(fig)

    expected_segments = int(
        pd.read_csv(panel_dir / "ideal_loadonly_segments.tsv", sep="\t").shape[0]
    )
    if total_segment_count != expected_segments:
        raise RuntimeError(
            f"Segment recount mismatch: {display}: {total_segment_count}!={expected_segments}"
        )
    return {
        "Panel": display,
        "PNG": str(png_path),
        "PDF": str(pdf_path),
        "PNG_DPI": 270,
        "Canvas_Inches": "20x8.8",
        "Zoom_Chrom": ZC,
        "Zoom_Red_Box_Count": len(zoom_segments),
        "Genome_Segment_Count": total_segment_count,
        "Connector_Line_Count": 0,
        "Current_Path_GWAS_Weight": 0,
        "GWAS_Position_Marks": gwas_summary["GWAS_Position_Marks"],
        "Projected_Gold_Marks": gwas_summary["Projected_Gold_Marks"],
        "Projected_Grey_Marks": gwas_summary["Projected_Grey_Marks"],
        "Marker_Status_Source": "Legacy_W4_Path_Not_Current_Path",
        "Legacy_Loci_SHA256": loci_sha256,
        "Legacy_Source_Script": str(SOURCE_SCRIPT),
        "Legacy_Source_SHA256": source_sha256,
        "Status": "LEGACY_PROJECTION_PASS",
    }


def main():
    args = parse_args()
    outroot = Path(args.outdir)
    if not SOURCE_SCRIPT.is_file():
        raise FileNotFoundError(SOURCE_SCRIPT)
    source_sha256 = hashlib.sha256(SOURCE_SCRIPT.read_bytes()).hexdigest()
    chromlen, chroms = read_fai(args.reference_fai)
    loci_path = Path(args.legacy_loci)
    by_chrom_path = Path(args.legacy_by_chrom)
    for path in (loci_path, by_chrom_path):
        if not path.is_file():
            raise FileNotFoundError(path)
    loci_sha256 = hashlib.sha256(loci_path.read_bytes()).hexdigest()
    loci_mark, legacy_counts, gwas_summary = read_legacy_gwas(
        loci_path, by_chrom_path, chroms, chromlen,
    )
    setup_style()
    panel_keys = [value.strip().lower() for value in args.panels.split(",") if value.strip()]
    if panel_keys != ["all38", "african35"]:
        raise ValueError("Formal legacy-v5 render requires panels=all38,african35")
    manifest_path = outroot / "legacy_v5_gwas_projection_manifest.tsv"
    provenance_path = outroot / "legacy_v5_gwas_projection_provenance.tsv"
    for path in (manifest_path, provenance_path):
        if path.exists():
            raise FileExistsError(f"Refusing to overwrite: {path}")
    rows = [
        render_panel(
            outroot, panel, chromlen, chroms, source_sha256,
            loci_mark, legacy_counts, gwas_summary, loci_sha256,
        )
        for panel in panel_keys
    ]
    fields = list(rows[0])
    with open(manifest_path, "x", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    provenance = [
        ("Visual_Source", str(SOURCE_SCRIPT), "Unchanged_reference"),
        ("Visual_Source_SHA256", source_sha256, "Recorded"),
        ("Canvas_Inches", "20x8.8", "Copied_verbatim"),
        ("GridSpec_Height_Ratios", "8,1.2", "Copied_verbatim"),
        ("GridSpec_Width_Ratios", "1.02,1.0", "Copied_verbatim"),
        ("GridSpec_Margins", "left=0.05;right=0.985;top=0.9;bottom=0.08;hspace=0.06;wspace=0.14", "Copied_verbatim"),
        ("Mosaic_Palette", ",".join(MOSAIC_COLORS), "Copied_verbatim"),
        ("Heatmap_Palette", ",".join(HEATMAP_COLORS), "Copied_verbatim"),
        ("Selected_Box_Color", SELECT_RED, "Copied_verbatim"),
        ("Selected_Box_Linewidth", "1.3", "Copied_verbatim"),
        ("Selected_Label_Color", SELECT_TEXT_RED, "Copied_verbatim"),
        ("Selected_Label_Fontsize", "4.8", "Copied_verbatim"),
        ("PNG_DPI", "270", "Copied_verbatim"),
        ("Captured_Marker_Color", CAPTURED_COLOR, "Copied_verbatim"),
        ("Missed_Marker_Color", MISSED_COLOR, "Copied_verbatim"),
        ("Data_Adapter", "Current_load_matrix,path,by_chrom", "Changed_data_only"),
        ("Current_Path_GWAS_Weight", "0", "Accepted_baseline_unchanged"),
        ("Legacy_Loci", str(loci_path), "Direct_coordinate_projection"),
        ("Legacy_Loci_SHA256", loci_sha256, "Recorded"),
        ("Legacy_Locus_Rows", str(gwas_summary["Legacy_Locus_Rows"]), "Validated"),
        ("GWAS_Position_Marks", str(gwas_summary["GWAS_Position_Marks"]), "Validated"),
        ("Projected_Gold_Marks", str(gwas_summary["Projected_Gold_Marks"]), "Legacy_status_only"),
        ("Projected_Grey_Marks", str(gwas_summary["Projected_Grey_Marks"]), "Legacy_status_only"),
        ("Marker_Status_Source", "Legacy_W4_Path_Not_Current_Path", "Explicit_caveat"),
    ]
    with open(provenance_path, "x", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["Item", "Value", "Status"])
        writer.writerows(provenance)
    print(
        "[OK] Legacy v5 GWAS-projection figures complete | "
        f"Panels={len(rows)} | Marks={gwas_summary['GWAS_Position_Marks']} | "
        "CurrentPathGWAS=0 | MarkerStatus=Legacy_W4"
    )


if __name__ == "__main__":
    main()
