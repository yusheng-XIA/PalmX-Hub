#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""18_compose_a4_single_row.py — A4 横版单行整基因组拼图（39 基因组 SyRI 链）

- 16 条染色体横向单行排列（chr01-16），同一比例尺（pt/Mb 完全一致），宽度严格等比。
- 面板来自 16_plot39_a4_landscape_panels.slurm（-s 1000000，仅保留 >=1Mb 彩色变异）。
- 横向比例 fx 由 A4 行宽决定；纵向 sy 默认在 A4 高度内取最大（≈13 pt/轨道），
  可用 --track-spacing-pt 或 --sy 覆盖。
- 面板以 raster（RASTER_PPI）非等比缩放插入，保证 16 条染色体一行放得下。
- 输出 A4 横版 PDF + 300dpi PNG + 布局/标签/审计/命令记录。

依赖：PyMuPDF（base 环境 ${DATA_DIR}/miniconda3/bin/python3）
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import shlex
import sys
from datetime import datetime, timezone
from pathlib import Path

import fitz

A4_W, A4_H = 518.740, 198.425  # 183mm x 70mm
MARGIN_L = 6.0
MARGIN_R = 6.0
MIN_MARGIN_Y = 3.0
LABEL_W = 38.0
GAP = 2.0
HEADER_H = 16.0
AXIS_H = 9.0
AXIS_TOTAL = 17.0
TRACK_LW = 0.85
TITLE_FS = 6.0
LABEL_FS = 4.0
TICK_FS = 6.0
CAPTION_FS = 6.0
LEGEND = (
    ("Syntenic", "#D8DDE3"),
    ("Inversion", "#E69F00"),
    ("Translocation", "#56B4E9"),
    ("Duplication", "#D55E00"),
)
BOLD_LABELS = {"nrly hap1", "nrly hap2", "Nigerian-Hap1", "Nigerian-Hap2"}
CHROMS = [f"chr{i:02d}" for i in range(1, 17)]
RASTER_PPI = 480.0


def rgb(value: str):
    value = value.lstrip("#")
    return tuple(int(value[i:i + 2], 16) / 255 for i in (0, 2, 4))


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_labels(path: Path, label_map: Path | None = None) -> list[str]:
    with path.open(encoding="utf-8") as handle:
        labels = [row[1] for row in csv.reader(handle, delimiter="\t") if len(row) >= 2]
    if len(labels) != 39:
        raise RuntimeError(f"Expected 39 labels, found {len(labels)}")
    if label_map is not None:
        mapping = {}
        with label_map.open(encoding="utf-8") as handle:
            for row in csv.reader(handle, delimiter="\t"):
                if len(row) >= 2 and row[0] != "old_label":
                    mapping[row[0]] = row[1]
        missing = [lab for lab in labels if lab not in mapping]
        if missing:
            raise RuntimeError(f"labels without new name: {missing}")
        labels = [mapping[lab] for lab in labels]
    return labels


def panel_ticks(page: fitz.Page, last_track_y: float) -> list[float]:
    values = []
    for block in page.get_text("dict")["blocks"]:
        for line in block.get("lines", []):
            for span in line["spans"]:
                if span["bbox"][1] <= last_track_y + 1:
                    continue
                try:
                    values.append(float(span["text"].strip()))
                except ValueError:
                    continue
    return sorted(set(values))


def measure_panel(page: fitz.Page):
    segs = []
    for drawing in page.get_drawings():
        color = drawing.get("color")
        width = drawing.get("width")
        if color is None or width is None or abs(width - TRACK_LW) > 0.06:
            continue
        for item in drawing["items"]:
            if item[0] != "l":
                continue
            a, b = item[1], item[2]
            if abs(a.y - b.y) < 0.05 and abs(a.x - b.x) >= 8:
                segs.append((round(a.y, 2), min(a.x, b.x), max(a.x, b.x)))
    if not segs:
        raise RuntimeError("no track-like segments found")
    groups = []
    for y, x0, x1 in sorted(segs):
        if groups and abs(groups[-1]["y"] - y) <= 0.6:
            g = groups[-1]
            g["segs"].append((x0, x1))
        else:
            groups.append({"y": y, "segs": [(x0, x1)]})
    for g in groups:
        g["x0"] = min(s[0] for s in g["segs"])
        g["x1"] = max(s[1] for s in g["segs"])
        g["maxseg"] = max(s[1] - s[0] for s in g["segs"])
    if len(groups) < 39:
        raise RuntimeError(f"only {len(groups)} line groups found; expected >= 39")
    ordered = sorted(groups, key=lambda g: -g["maxseg"])
    tracks = ordered[:39]
    if len(ordered) > 39 and ordered[39]["maxseg"] > 0.5 * ordered[38]["maxseg"]:
        raise RuntimeError("ambiguous track/decoration separation")
    tracks.sort(key=lambda g: g["y"])
    ys = [g["y"] for g in tracks]
    if max(g["x0"] for g in tracks) - min(g["x0"] for g in tracks) > 0.2:
        raise RuntimeError("track left edges are not aligned")
    x0 = min(g["x0"] for g in tracks)
    if len(tracks[0]["segs"]) != 1:
        raise RuntimeError("reference track line is not a single segment")
    ref_width = tracks[0]["segs"][0][1] - tracks[0]["segs"][0][0]
    full_width = max(g["x1"] - g["x0"] for g in tracks)
    diffs = [ys[i + 1] - ys[i] for i in range(38)]
    t = sum(diffs) / 38
    if (max(diffs) - min(diffs)) / t > 5e-3:
        raise RuntimeError(f"uneven track spacing: {min(diffs):.3f}-{max(diffs):.3f}")
    return {"ys": ys, "t": t, "x0": x0, "ref_width": ref_width,
            "full_width": full_width, "ticks": panel_ticks(page, ys[-1])}


def nice_step(len_mb: float, pt_per_mb: float) -> int:
    for step in (5, 10, 20, 25, 50, 100):
        if step * pt_per_mb >= 12.0 and 3 <= len_mb // step + 1 <= 6:
            return step
    candidates = []
    for step in (5, 10, 20, 25, 50, 100):
        n = len_mb // step + 1
        if n >= 2:
            candidates.append((step * pt_per_mb, step))
    return max(candidates)[1] if candidates else 10


def fit_ticks(measured: list[float], len_mb: float, scale_panel: float) -> list[float]:
    ticks = [v for v in measured if 0.0 <= v <= len_mb]
    if len(ticks) >= 2:
        widths = [fitz.get_text_length(str(int(round(v))), fontname="helv",
                                       fontsize=TICK_FS) for v in ticks]
        crowded = any((ticks[i + 1] - ticks[i]) * scale_panel
                      < max(widths[i], widths[i + 1]) + 2.0
                      for i in range(len(ticks) - 1))
        if not crowded:
            return ticks
    step = nice_step(len_mb, scale_panel)
    return [k * step for k in range(int(len_mb // step) + 1)]


def draw_legend(page: fitz.Page, x: float, y: float) -> None:
    page.insert_text((x, y), "Structural relationship", fontsize=6.0,
                     fontname="hebo", color=(0.12, 0.12, 0.12))
    x2 = x + 78
    for label, color in LEGEND:
        page.draw_rect(fitz.Rect(x2, y - 4.5, x2 + 8, y + 1.5),
                       color=rgb(color), fill=rgb(color))
        page.insert_text((x2 + 11, y), label, fontsize=6.0, fontname="helv",
                         color=(0.12, 0.12, 0.12))
        x2 += 50


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", required=True, type=Path)
    ap.add_argument("--panels", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--dpi", type=int, default=300)
    ap.add_argument("--scale-tol", type=float, default=5e-3,
                    help="relative tolerance for cross-panel pt/Mb consistency")
    ap.add_argument("--labels-tsv", type=Path, default=None,
                    help="old_label->new_label TSV (latest display names)")
    ap.add_argument("--sy", type=float, default=None,
                    help="vertical scale (default: maximum that fits A4 height)")
    ap.add_argument("--track-spacing-pt", type=float, default=None,
                    help="target track spacing in pt (overrides --sy)")
    args = ap.parse_args()

    args.out.mkdir(parents=True, exist_ok=True)
    (args.out / "COMMAND.txt").write_text(
        shlex.join([sys.executable, *sys.argv]) + "\n", encoding="utf-8")

    labels = read_labels(args.run / "results/07_plot39_plotsr_inputs_linked/genomes.txt",
                         args.labels_tsv)
    lengths = {}
    max_lengths = {}
    panels = {}
    for chrom in CHROMS:
        src_dir = args.panels / chrom
        candidates = sorted(src_dir.glob("*.pdf"))
        if not candidates:
            raise RuntimeError(f"missing panel pdf: {src_dir}")
        doc = fitz.open(candidates[0])
        measured = measure_panel(doc[0])
        measured["doc"] = doc
        measured["path"] = candidates[0]
        panels[chrom] = measured
        lengths[chrom] = int((src_dir / "length_ref_bp.txt").read_text().strip())
        max_lengths[chrom] = int((src_dir / "length_max_bp.txt").read_text().strip())

    ts = [panels[c]["t"] for c in panels]
    t_mean = sum(ts) / len(ts)
    if (max(ts) - min(ts)) / t_mean > 5e-3:
        raise RuntimeError(f"track spacing mismatch across panels: {min(ts):.3f}-{max(ts):.3f}")

    scales = {c: panels[c]["ref_width"] / (lengths[c] / 1e6) for c in panels}
    s_mean = sum(scales.values()) / len(scales)
    if (max(scales.values()) - min(scales.values())) / s_mean > args.scale_tol:
        raise RuntimeError(f"per-panel pt/Mb mismatch: {min(scales.values()):.4f}-{max(scales.values()):.4f}")

    row_avail = A4_W - MARGIN_L - MARGIN_R - LABEL_W
    ref_sum = sum(panels[c]["ref_width"] for c in CHROMS)
    fx = (row_avail - GAP * (len(CHROMS) - 1)) / ref_sum
    sy_max = (A4_H - 2 * MIN_MARGIN_Y - HEADER_H - AXIS_TOTAL) / (39 * t_mean)
    if args.track_spacing_pt is not None:
        sy = args.track_spacing_pt / t_mean
    elif args.sy is not None:
        sy = args.sy
    else:
        sy = sy_max
    body_h = 39 * t_mean * sy
    content_h = HEADER_H + body_h + AXIS_TOTAL
    if content_h > A4_H:
        raise RuntimeError(f"body does not fit A4: content {content_h:.2f} > {A4_H}")
    y0 = max(MIN_MARGIN_Y, (A4_H - content_h) / 2)

    out = fitz.open()
    page = out.new_page(width=A4_W, height=A4_H)
    subtitle = "39 genomes; two references excluded; minimum SR 1 Mb"
    sub_w = fitz.get_text_length(subtitle, fontname="helv", fontsize=6.0)
    page.insert_text((A4_W - MARGIN_R - sub_w, y0 + 6),
                     subtitle,
                     fontsize=6.0, fontname="helv", color=(0.35, 0.35, 0.35))
    draw_legend(page, MARGIN_L + LABEL_W, y0 + 6)

    y_band = y0 + HEADER_H
    body_bottom = y_band + body_h
    axis_y = body_bottom + 2.5
    x = MARGIN_L + LABEL_W
    layout = []
    track_map = []
    for j, chrom in enumerate(CHROMS):
        m = panels[chrom]
        clip = fitz.Rect(m["x0"], m["ys"][0] - m["t"] / 2,
                         m["x0"] + m["ref_width"], m["ys"][-1] + m["t"] / 2)
        w = m["ref_width"] * fx
        h = clip.height * sy
        dest = fitz.Rect(x, y_band, x + w, y_band + h)
        pix = m["doc"][0].get_pixmap(
            matrix=fitz.Matrix(fx * RASTER_PPI / 72, sy * RASTER_PPI / 72),
            clip=clip, alpha=False)
        page.insert_image(dest, pixmap=pix)
        page.insert_text((x + w / 2 - 7, y_band - 1.5), chrom, fontsize=6.0,
                         fontname="hebo", color=(0.08, 0.08, 0.08))
        if j > 0:
            page.draw_line((x - GAP / 2, y_band), (x - GAP / 2, body_bottom),
                           color=(0.85, 0.85, 0.85), width=0.35)
        len_mb = lengths[chrom] / 1e6
        scale_panel = w / len_mb
        ticks = fit_ticks(m["ticks"], len_mb, scale_panel)
        page.draw_line((x, axis_y), (x + w, axis_y), color=(0.55, 0.55, 0.55), width=0.4)
        for value in ticks:
            tx = x + value * scale_panel
            page.draw_line((tx, axis_y - 1.8), (tx, axis_y), color=(0.55, 0.55, 0.55), width=0.4)
            label = str(int(round(value)))
            tw = fitz.get_text_length(label, fontname="helv", fontsize=TICK_FS)
            page.insert_text((tx - tw / 2, axis_y + 4.8), label, fontsize=TICK_FS,
                             fontname="helv", color=(0.3, 0.3, 0.3))
        layout.append((j + 1, chrom, str(m["path"]),
                       round(x, 2), round(w, 3), round(len_mb, 3),
                       round(max_lengths[chrom] / 1e6, 3),
                       round(m["full_width"] / m["ref_width"], 4),
                       round(scale_panel, 6)))
        x += w + GAP

    for label, ty in zip(labels, panels[CHROMS[0]]["ys"]):
        dy = y_band + (ty - (panels[CHROMS[0]]["ys"][0] - panels[CHROMS[0]]["t"] / 2)) * sy
        font = "hebo" if label in BOLD_LABELS else "helv"
        text_w = fitz.get_text_length(label, fontname=font, fontsize=LABEL_FS)
        page.insert_text((MARGIN_L + LABEL_W - 3 - text_w, dy + LABEL_FS / 3), label,
                         fontsize=LABEL_FS, fontname=font, color=(0.08, 0.08, 0.08))
        track_map.append((label, round(ty, 3), round(dy, 3)))
    page.draw_line((MARGIN_L + LABEL_W - 1.5, y_band),
                   (MARGIN_L + LABEL_W - 1.5, body_bottom),
                   color=(0.76, 0.76, 0.76), width=0.4)

    row_w = x - GAP - (MARGIN_L + LABEL_W)
    caption = "Chromosome position (Mb)"
    cw = fitz.get_text_length(caption, fontname="helv", fontsize=CAPTION_FS)
    page.insert_text((MARGIN_L + LABEL_W + row_w / 2 - cw / 2, axis_y + 13.2),
                     caption, fontsize=CAPTION_FS, fontname="helv",
                     color=(0.35, 0.35, 0.35))

    out.set_metadata({
        "title": "Oil palm 39-genome SyRI chain, whole genome, 183x70mm, single row 16x1",
        "subject": ("chr01-16 in one row; rowfill; nrly hap2 V3; minimum SR 1 Mb; "
                    "latest display names"),
        "creator": "plotsr (scaled panels) + PyMuPDF 183x70mm single-row compositor",
    })
    pdf = args.out / ("OilPalm_39genomes_wholegenome_183x70mm_16x1_rowfill_"
                      "MAJOR1Mb_NRLY_HAP2_V3_COORD.pdf")
    png = args.out / ("OilPalm_39genomes_wholegenome_183x70mm_16x1_rowfill_"
                      "MAJOR1Mb_NRLY_HAP2_V3_COORD.png")
    out.save(pdf, garbage=4, deflate=True)
    out.close()
    for m in panels.values():
        m["doc"].close()

    check = fitz.open(pdf)
    text = check[0].get_text()
    import re
    missing_chr = [f"chr{i:02d}" for i in range(1, 17) if f"chr{i:02d}" not in text]
    label_counts = {label: len(re.findall(
        r"(?<![\w.\-])" + re.escape(label) + r"(?![\w.\-])", text)) for label in labels}
    wrong = {label: count for label, count in label_counts.items() if count != 1}
    if missing_chr or wrong:
        raise RuntimeError(f"audit failure missing_chr={missing_chr} wrong_label_counts={wrong}")
    pix = check[0].get_pixmap(dpi=args.dpi, alpha=False)
    pix.save(png)
    check.close()

    with (args.out / "PANEL_LAYOUT.tsv").open("w", encoding="utf-8") as handle:
        handle.write("order\tchromosome\tsource_pdf\tx_pt\twidth_pt\tref_length_mb\t"
                     "max_homolog_mb\toverhang_ratio\tpt_per_mb\n")
        for item in layout:
            handle.write("\t".join(map(str, item)) + "\n")
    with (args.out / "TRACK_LABEL_MAPPING.tsv").open("w", encoding="utf-8") as handle:
        handle.write("genome\tsource_track_y_pt\tcomposite_label_y_pt\n")
        for item in track_map:
            handle.write("\t".join(map(str, item)) + "\n")
    with (args.out / "DISPLAY_LABELS_39.tsv").open("w", encoding="utf-8") as handle:
        handle.write("order\tlabel\n")
        for i, lab in enumerate(labels, 1):
            handle.write(f"{i}\t{lab}\n")
    pt_per_mb = [item[8] for item in layout]
    rel = (max(pt_per_mb) - min(pt_per_mb)) / (sum(pt_per_mb) / len(pt_per_mb))
    with (args.out / "FINAL_AUDIT.tsv").open("w", encoding="utf-8") as handle:
        handle.write("metric\tvalue\n")
        for key, value in (
            ("status", "passed"), ("page", "183x70mm 518.740x198.425pt"),
            ("layout", "16x1 single row chr01-16"), ("displayed_genomes", 39),
            ("bands_mode", "16x1_rowfill"),
            ("min_sr_bp", 1000000), ("track_spacing_pt", round(t_mean * sy, 3)),
            ("composite_scale_fx", round(fx, 6)),
            ("composite_scale_sy", round(sy, 6)),
            ("composite_scale_sy_max", round(sy_max, 6)),
            ("pt_per_mb", round(sum(pt_per_mb) / len(pt_per_mb), 6)),
            ("pt_per_mb_rel_diff", round(rel, 6)),
            ("row_width_pt", round(row_w, 3)),
            ("reference", "EG_146; panels cropped at reference chromosome end"),
            ("body_height_pt", round(body_h, 3)), ("label_font_pt", LABEL_FS),
            ("labels_tsv_source", str(args.labels_tsv) if args.labels_tsv else "none"),
            ("labels_tsv_sha256", sha256(args.labels_tsv) if args.labels_tsv else "none"),
            ("pymupdf_version", fitz.VersionBind),
            ("composed_utc", datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")),
            ("png_dpi", args.dpi),
            ("pdf_sha256", sha256(pdf)), ("png_sha256", sha256(png)),
        ):
            handle.write(f"{key}\t{value}\n")
    (args.out / "SUCCESS").touch()
    print(f"PASS {pdf}")
    print(f"PASS {png}")
    print(f"fx={fx:.6f} sy={sy:.6f} (sy_max={sy_max:.6f}) track={t_mean * sy:.3f}pt "
          f"body={body_h:.2f}pt row_w={row_w:.2f}pt")


if __name__ == "__main__":
    main()
