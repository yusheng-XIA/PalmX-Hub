#!/usr/bin/env python3
"""Summarize vcftools pi/FST windows and draw simple SVG figures."""

import glob
import math
import os
from pathlib import Path


OUTDIR = Path("${ANALYSIS_DIR}/22_answer_reviews/00_ms/03_V3/04_figure4/K3_K4_pi_fst_pca_groups")
NA_SET = {"na", "nan", "-nan", "none", "."}


def safe_float(x):
    s = x.strip()
    if not s or s.lower() in NA_SET:
        return None
    try:
        value = float(s)
    except ValueError:
        return None
    if math.isnan(value) or math.isinf(value):
        return None
    return value


def read_rows(path):
    header = None
    rows = []
    with path.open(errors="replace") as handle:
        for line in handle:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            if header is None:
                header = line.split()
            else:
                rows.append(line.split())
    return header or [], rows


def col_index(header, names):
    for name in names:
        if name in header:
            return header.index(name)
    return None


def weighted_mean(values, weights):
    numerator = 0.0
    denominator = 0.0
    used = 0
    for value, weight in zip(values, weights):
        if value is None or weight is None or weight <= 0:
            continue
        numerator += value * weight
        denominator += weight
        used += 1
    if denominator <= 0:
        return None, used
    return numerator / denominator, used


def simple_mean(values):
    clean = [x for x in values if x is not None]
    if not clean:
        return None, 0
    return sum(clean) / len(clean), len(clean)


def fmt(value):
    return "NA" if value is None else f"{value:.12g}"


def normalize_name(path):
    name = path.name
    for suffix in ["_100kb_fst.windowed.weir.fst", "_100kb.pi.windowed.pi"]:
        if name.endswith(suffix):
            return name[: -len(suffix)]
    return name


def summarize_pi(path, k_label):
    header, rows = read_rows(path)
    pi_idx = col_index(header, ["PI"])
    nvar_idx = col_index(header, ["N_VARIANTS"])
    nmono_idx = col_index(header, ["N_MONOMORPHIC", "N_MONO", "N_INVARIANT"])
    if pi_idx is None:
        raise ValueError(f"No PI column in {path}")

    values = []
    weights = []
    for row in rows:
        values.append(safe_float(row[pi_idx]) if len(row) > pi_idx else None)
        nvar = safe_float(row[nvar_idx]) if nvar_idx is not None and len(row) > nvar_idx else None
        nmono = safe_float(row[nmono_idx]) if nmono_idx is not None and len(row) > nmono_idx else None
        weights.append((nvar + nmono) if nvar is not None and nmono is not None else None)
    sm, ns = simple_mean(values)
    wm, nw = weighted_mean(values, weights)
    return [k_label, normalize_name(path), "PI", "PI", fmt(sm), fmt(wm), "N_VARIANTS+N_MONOMORPHIC", ns, nw, len(rows)]


def summarize_fst(path, k_label):
    header, rows = read_rows(path)
    fst_idx = col_index(header, ["MEAN_FST", "WEIGHTED_FST", "WEIR_AND_COCKERHAM_FST", "FST"])
    nvar_idx = col_index(header, ["N_VARIANTS", "N_SNPS", "N_VARIANT", "N_SITES"])
    if fst_idx is None:
        raise ValueError(f"No FST column in {path}")

    values = []
    weights = []
    for row in rows:
        values.append(safe_float(row[fst_idx]) if len(row) > fst_idx else None)
        weights.append(safe_float(row[nvar_idx]) if nvar_idx is not None and len(row) > nvar_idx else None)
    sm, ns = simple_mean(values)
    wm, nw = weighted_mean(values, weights)
    return [k_label, normalize_name(path), "FST", header[fst_idx], fmt(sm), fmt(wm), "N_VARIANTS", ns, nw, len(rows)]


def read_summary_table(path):
    with path.open() as handle:
        header = handle.readline().rstrip("\n").split("\t")
        rows = []
        for line in handle:
            if line.strip():
                rows.append(dict(zip(header, line.rstrip("\n").split("\t"))))
    return rows


def write_table(path, header, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as handle:
        handle.write("\t".join(header) + "\n")
        for row in rows:
            handle.write("\t".join(str(x) for x in row) + "\n")


def build_matrices(summary_rows):
    for k_label in ["K3", "K4"]:
        pi_rows = [r for r in summary_rows if r[0] == k_label and r[2] == "PI"]
        fst_rows = [r for r in summary_rows if r[0] == k_label and r[2] == "FST"]

        pi_out = [[r[1], r[5], r[8], r[9]] for r in sorted(pi_rows, key=lambda x: x[1])]
        write_table(OUTDIR / "summary" / f"{k_label}_PI_weighted.tsv", ["group", "PI_weighted", "n_used_weighted", "n_windows"], pi_out)

        groups = sorted([r[0] for r in pi_out])
        matrix = {g: {h: "NA" for h in groups} for g in groups}
        for g in groups:
            matrix[g][g] = "0"
        long_rows = []
        for row in fst_rows:
            name = str(row[1])
            left = name.rsplit("_100kb_fst", 1)[0]
            parts = left.split("_")
            if len(parts) < 2:
                continue
            pop1 = "_".join(parts[:2]) if parts[0].startswith("K") and parts[1].startswith("Pop") else parts[0]
            pop2 = "_".join(parts[2:4]) if len(parts) >= 4 and parts[2].startswith("K") and parts[3].startswith("Pop") else parts[-1]
            value = str(row[5])
            if pop1 in matrix and pop2 in matrix:
                matrix[pop1][pop2] = value
                matrix[pop2][pop1] = value
                long_rows.append([pop1, pop2, value, row[8], row[9]])

        matrix_rows = [[g] + [matrix[g][h] for h in groups] for g in groups]
        write_table(OUTDIR / "summary" / f"{k_label}_FST_weighted_matrix.tsv", [""] + groups, matrix_rows)
        write_table(OUTDIR / "summary" / f"{k_label}_FST_weighted_long.tsv", ["pop1", "pop2", "FST_weighted", "n_used_weighted", "n_windows"], long_rows)


def svg_escape(text):
    return text.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")


def draw_bar_svg(k_label, pi_table):
    rows = read_summary_table(pi_table)
    labels = [r["group"] for r in rows]
    values = [float(r["PI_weighted"]) for r in rows]
    width, height = 760, 480
    left, right, top, bottom = 90, 30, 60, 90
    plot_w = width - left - right
    plot_h = height - top - bottom
    vmax = max(values) * 1.18 if values else 1.0
    colors = ["#4E79A7", "#F28E2B", "#59A14F", "#E15759"]
    bar_w = plot_w / max(len(values), 1) * 0.62
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="white"/>',
        f'<text x="{width/2}" y="30" text-anchor="middle" font-family="Arial" font-size="22" font-weight="700">{k_label} nucleotide diversity</text>',
        f'<line x1="{left}" y1="{top + plot_h}" x2="{left + plot_w}" y2="{top + plot_h}" stroke="#333"/>',
        f'<line x1="{left}" y1="{top}" x2="{left}" y2="{top + plot_h}" stroke="#333"/>',
    ]
    for i in range(6):
        value = vmax * i / 5
        y = top + plot_h - value / vmax * plot_h
        parts.append(f'<line x1="{left-5}" y1="{y:.1f}" x2="{left + plot_w}" y2="{y:.1f}" stroke="#e5e5e5"/>')
        parts.append(f'<text x="{left-10}" y="{y+4:.1f}" text-anchor="end" font-family="Arial" font-size="12">{value:.4f}</text>')
    for idx, (label, value) in enumerate(zip(labels, values)):
        cx = left + plot_w * (idx + 0.5) / len(values)
        bar_h = value / vmax * plot_h
        x = cx - bar_w / 2
        y = top + plot_h - bar_h
        color = colors[idx % len(colors)]
        parts.append(f'<rect x="{x:.1f}" y="{y:.1f}" width="{bar_w:.1f}" height="{bar_h:.1f}" fill="{color}"/>')
        parts.append(f'<text x="{cx:.1f}" y="{y-8:.1f}" text-anchor="middle" font-family="Arial" font-size="12">{value:.5f}</text>')
        parts.append(f'<text x="{cx:.1f}" y="{top + plot_h + 25}" text-anchor="middle" font-family="Arial" font-size="13">{svg_escape(label)}</text>')
    parts.append(f'<text x="24" y="{top + plot_h/2}" transform="rotate(-90 24 {top + plot_h/2})" text-anchor="middle" font-family="Arial" font-size="15">Weighted pi</text>')
    parts.append("</svg>\n")
    (OUTDIR / "figures" / f"{k_label}_PI_weighted_bar.svg").write_text("\n".join(parts))


def color_ramp(value, vmin, vmax):
    if vmax <= vmin:
        t = 0.0
    else:
        t = (value - vmin) / (vmax - vmin)
    start = (247, 251, 255)
    end = (8, 81, 156)
    rgb = tuple(round(start[i] + t * (end[i] - start[i])) for i in range(3))
    return f"#{rgb[0]:02x}{rgb[1]:02x}{rgb[2]:02x}"


def draw_heatmap_svg(k_label, matrix_table):
    with matrix_table.open() as handle:
        header = handle.readline().rstrip("\n").split("\t")[1:]
        rows = []
        values = []
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            label = fields[0]
            row_vals = fields[1:]
            rows.append((label, row_vals))
            for value in row_vals:
                if value != "NA" and value != "0":
                    values.append(float(value))
    n = len(header)
    cell = 94
    left, top = 130, 80
    width = left + n * cell + 40
    height = top + n * cell + 70
    vmin, vmax = (min(values), max(values)) if values else (0.0, 1.0)
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="white"/>',
        f'<text x="{width/2}" y="34" text-anchor="middle" font-family="Arial" font-size="22" font-weight="700">{k_label} pairwise FST</text>',
    ]
    for j, label in enumerate(header):
        x = left + j * cell + cell / 2
        parts.append(f'<text x="{x:.1f}" y="{top-18}" text-anchor="middle" font-family="Arial" font-size="13">{svg_escape(label)}</text>')
    for i, (label, row_vals) in enumerate(rows):
        y = top + i * cell + cell / 2
        parts.append(f'<text x="{left-12}" y="{y+4:.1f}" text-anchor="end" font-family="Arial" font-size="13">{svg_escape(label)}</text>')
        for j, value in enumerate(row_vals):
            x0 = left + j * cell
            y0 = top + i * cell
            if value == "NA":
                fill = "#f0f0f0"
                text = "NA"
            else:
                numeric = float(value)
                fill = "#ffffff" if numeric == 0 else color_ramp(numeric, vmin, vmax)
                text = "0" if numeric == 0 else f"{numeric:.4f}"
            parts.append(f'<rect x="{x0}" y="{y0}" width="{cell}" height="{cell}" fill="{fill}" stroke="#ffffff"/>')
            parts.append(f'<text x="{x0 + cell/2:.1f}" y="{y0 + cell/2 + 4:.1f}" text-anchor="middle" font-family="Arial" font-size="12" fill="#111">{text}</text>')
    parts.append(f'<text x="{left}" y="{height-24}" font-family="Arial" font-size="12">Color scale excludes diagonal zeros; darker blue indicates larger weighted FST.</text>')
    parts.append("</svg>\n")
    (OUTDIR / "figures" / f"{k_label}_FST_weighted_heatmap.svg").write_text("\n".join(parts))


def main():
    summary_rows = []
    for k_label in ["K3", "K4"]:
        for path in sorted(Path(OUTDIR / "vcftools" / k_label).glob("*_100kb.pi.windowed.pi")):
            summary_rows.append(summarize_pi(path, k_label))
        for path in sorted(Path(OUTDIR / "vcftools" / k_label).glob("*_100kb_fst.windowed.weir.fst")):
            summary_rows.append(summarize_fst(path, k_label))

    if not summary_rows:
        raise SystemExit("No vcftools pi/FST output files found yet.")

    write_table(
        OUTDIR / "summary" / "K3_K4_pi_fst_mean.tsv",
        ["K", "name", "metric", "value_column", "simple_mean", "weighted_mean", "weight_by", "n_used_simple", "n_used_weighted", "n_windows"],
        summary_rows,
    )
    build_matrices(summary_rows)
    for k_label in ["K3", "K4"]:
        draw_bar_svg(k_label, OUTDIR / "summary" / f"{k_label}_PI_weighted.tsv")
        draw_heatmap_svg(k_label, OUTDIR / "summary" / f"{k_label}_FST_weighted_matrix.tsv")
    print(f"Wrote summaries and SVG figures under {OUTDIR}")


if __name__ == "__main__":
    main()
