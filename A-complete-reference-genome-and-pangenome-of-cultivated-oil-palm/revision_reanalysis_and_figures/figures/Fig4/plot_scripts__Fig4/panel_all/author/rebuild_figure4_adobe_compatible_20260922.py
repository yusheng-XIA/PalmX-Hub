#!/usr/bin/env python3
"""Build clean standalone Figure 4 panels for Adobe interoperability.

The reviewed independent panel artwork is placed on a fresh SVG canvas and
converted to PDF with CairoSVG. This avoids the previous PyPDF2 operation that
retained the complete composite page inside every clipped panel. Scientific
values and encodings are not recalculated or edited.
"""

from __future__ import annotations

import copy
import hashlib
import json
import os
import platform
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

os.environ.setdefault("MPLBACKEND", "Agg")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")

import cairosvg
from lxml import etree
from PIL import Image, ImageDraw


RUN = Path(
    "${ANALYSIS_DIR}/"
    "22_answer_reviews/00_ms/05_MS/0918_revision/final_ms/02_figure4/"
    "redrawn_adobe_compatible_20260922"
)
PANELS = RUN / "panels"
SOURCE_NORMALIZED = RUN / "source_panels"
QC = RUN / "qc"
PROV = RUN / "provenance"

BASE = Path(
    "${ANALYSIS_DIR}/"
    "22_answer_reviews/00_ms/03_V3/04_figure4"
)
FINAL = RUN.parent

SOURCES = {
    "a": BASE / "K3_K4_pi_fst_pca_groups/figures/requested_01_tree_structure_K3_K4_K6_K4treecolor.pdf",
    "b": BASE / "K3_K4_pi_fst_pca_groups/figures/requested_02_K4_population_PCA_PC1_PC2.pdf",
    "c": BASE / "K3_K4_pi_fst_pca_groups/figures/requested_03_K4_pi_fst_network.pdf",
    "d": FINAL / "Figure4d_60x38.5mm_vector_solid_dashed_v3.pdf",
    "e": BASE / "Fig4_d_i_pan39_material33_singletons_20260811/figures/Fig4e_pan_core_genome39_labelled.svg",
    "f": BASE / "Fig4_d_i_pan39_material33_singletons_20260811/figures/Fig4f_frequency_donut_genome39_labelled.svg",
    "g": BASE / "Fig4_d_i_pan39_material33_singletons_20260811/figures/Fig4g_allele_composition_material33_labelled.svg",
    "h": BASE / "Fig4_d_i_pan39_material33_singletons_20260811/figures/Fig4h_WGD_Ks_KaKs_material33_labelled.svg",
    "i": BASE / "Fig4_d_i_pan39_material33_singletons_20260811/figures/Fig4i_TE_density_pan39_material33_labelled.svg",
    "j": BASE / "RGA_tree_dotmatrix_redesign_20260708/figures/RGA_tree_dotmatrix_major_groups.svg",
}

TARGET_MM = {panel: (60.0, 38.5) for panel in "abcdefghi"}
TARGET_MM["j"] = (91.0, 72.2)
ADD_LETTER = {"a", "b", "c", "j"}
SVG_NS = "http://www.w3.org/2000/svg"
NSMAP = {None: SVG_NS, "xlink": "http://www.w3.org/1999/xlink"}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def run_command(args: list[str]) -> None:
    subprocess.run(args, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)


def pdf_to_svg(panel: str, source: Path) -> Path:
    output = SOURCE_NORMALIZED / f"Figure4{panel}_independent_source.svg"
    run_command(["pdftocairo", "-svg", str(source), str(output)])
    if not output.is_file() or output.stat().st_size == 0:
        raise RuntimeError(f"pdftocairo did not create {output}")
    return output


def parse_viewbox(root: etree._Element) -> tuple[float, float, float, float]:
    raw = root.get("viewBox")
    if not raw:
        raise ValueError("Source SVG lacks a viewBox")
    values = tuple(float(value) for value in raw.replace(",", " ").split())
    if len(values) != 4 or values[2] <= 0 or values[3] <= 0:
        raise ValueError(f"Invalid source viewBox: {raw}")
    return values


def remove_j_heading(old_root: etree._Element) -> None:
    # IDs are stable in the reviewed Matplotlib export. Removing only these
    # three groups drops the old panel-a letter, title and subtitle. Plot data,
    # legends and annotations remain untouched.
    for element_id in ("text_54", "text_55", "text_56"):
        matches = old_root.xpath(f".//*[@id='{element_id}']")
        if len(matches) != 1:
            raise RuntimeError(f"Expected one {element_id} in panel-j source; found {len(matches)}")
        node = matches[0]
        node.getparent().remove(node)


def build_outer_svg(panel: str, source_svg: Path, output_svg: Path) -> dict[str, object]:
    parser = etree.XMLParser(resolve_entities=False, no_network=True, remove_blank_text=False)
    old_root = etree.parse(str(source_svg), parser).getroot()
    x0, y0, source_w, source_h = parse_viewbox(old_root)
    if panel == "j":
        remove_j_heading(old_root)
        # Remove the now-empty title band while preserving the plot and legends.
        y0 += 32.0
        source_h -= 32.0

    width_mm, height_mm = TARGET_MM[panel]
    canvas_w = width_mm * 10.0
    canvas_h = height_mm * 10.0
    root = etree.Element(
        f"{{{SVG_NS}}}svg",
        nsmap=NSMAP,
        width=f"{width_mm:g}mm",
        height=f"{height_mm:g}mm",
        viewBox=f"0 0 {canvas_w:g} {canvas_h:g}",
        version="1.1",
    )
    title = etree.SubElement(root, f"{{{SVG_NS}}}title")
    title.text = f"Figure 4{panel}: clean Adobe-compatible standalone vector panel"
    desc = etree.SubElement(root, f"{{{SVG_NS}}}desc")
    desc.text = "Rebuilt from reviewed independent vector artwork; scientific values unchanged."
    etree.SubElement(
        root,
        f"{{{SVG_NS}}}rect",
        x="0",
        y="0",
        width=f"{canvas_w:g}",
        height=f"{canvas_h:g}",
        fill="#ffffff",
    )

    nested = copy.deepcopy(old_root)
    nested.set("x", "0")
    nested.set("y", "0")
    nested.set("width", f"{canvas_w:g}")
    nested.set("height", f"{canvas_h:g}")
    nested.set("viewBox", f"{x0:g} {y0:g} {source_w:g} {source_h:g}")
    # Match the accepted Illustrator assembly geometry: fill the exact panel
    # canvas. This does not alter any plotted coordinate or value.
    nested.set("preserveAspectRatio", "none")
    root.append(nested)

    if panel in ADD_LETTER:
        letter = etree.SubElement(
            root,
            f"{{{SVG_NS}}}text",
            x=f"{canvas_w * 0.012:g}",
            y=f"{canvas_h * 0.085:g}",
            fill="#111111",
            style="font-family:Arial,'Liberation Sans',sans-serif;font-size:30px;font-weight:700",
        )
        letter.text = panel

    tree = etree.ElementTree(root)
    tree.write(str(output_svg), encoding="utf-8", xml_declaration=True, pretty_print=True)
    return {
        "source_svg": str(source_svg),
        "source_viewbox": [x0, y0, source_w, source_h],
        "target_mm": [width_mm, height_mm],
        "geometry": "independent x/y scaling to exact final canvas, matching accepted assembly",
    }


def make_contact_sheet(pngs: list[Path]) -> Path:
    cell_w, cell_h = 900, 650
    sheet = Image.new("RGB", (cell_w * 2, cell_h * 5), "#dddddd")
    font = None
    for index, path in enumerate(pngs):
        image = Image.open(path).convert("RGB")
        image.thumbnail((cell_w - 50, cell_h - 90))
        cell = Image.new("RGB", (cell_w, cell_h), "white")
        x = (cell_w - image.width) // 2
        y = 55 + (cell_h - 75 - image.height) // 2
        cell.paste(image, (x, y))
        ImageDraw.Draw(cell).text((15, 15), path.stem, fill="black", font=font)
        sheet.paste(cell, ((index % 2) * cell_w, (index // 2) * cell_h))
    output = QC / "Figure4_redrawn_adobe_compatible_contact_sheet.png"
    sheet.save(output)
    return output


def main() -> None:
    started = datetime.now(timezone.utc)
    for directory in (PANELS, SOURCE_NORMALIZED, QC, PROV):
        directory.mkdir(parents=True, exist_ok=True)
    for panel, source in SOURCES.items():
        if not source.is_file() or source.stat().st_size == 0:
            raise FileNotFoundError(source)

    input_manifest = []
    output_manifest = []
    panel_records: dict[str, object] = {}
    pngs: list[Path] = []

    for panel in "abcdefghij":
        source = SOURCES[panel]
        input_manifest.append((sha256(source), str(source)))
        if source.suffix.lower() == ".pdf":
            source_svg = pdf_to_svg(panel, source)
        else:
            source_svg = source

        stem = (
            f"Figure4{panel}_60x38.5mm_redrawn_adobe"
            if panel != "j"
            else "Figure4j_91x72.2mm_redrawn_adobe"
        )
        svg_path = PANELS / f"{stem}.svg"
        pdf_path = PANELS / f"{stem}.pdf"
        png_path = PANELS / f"{stem}_600dpi.png"
        record = build_outer_svg(panel, source_svg, svg_path)
        width_mm, height_mm = TARGET_MM[panel]
        cairosvg.svg2pdf(url=str(svg_path), write_to=str(pdf_path))
        cairosvg.svg2png(
            url=str(svg_path),
            write_to=str(png_path),
            output_width=round(width_mm / 25.4 * 600),
            output_height=round(height_mm / 25.4 * 600),
        )
        for output in (svg_path, pdf_path, png_path):
            if not output.is_file() or output.stat().st_size == 0:
                raise RuntimeError(f"Missing output: {output}")
            output_manifest.append((sha256(output), str(output)))
        record.update(
            {
                "source": str(source),
                "source_sha256": sha256(source),
                "svg": str(svg_path),
                "pdf": str(pdf_path),
                "png_600dpi": str(png_path),
                "svg_sha256": sha256(svg_path),
                "pdf_sha256": sha256(pdf_path),
                "png_sha256": sha256(png_path),
            }
        )
        panel_records[panel] = record
        pngs.append(png_path)

    contact_sheet = make_contact_sheet(pngs)
    output_manifest.append((sha256(contact_sheet), str(contact_sheet)))
    completed = datetime.now(timezone.utc)
    record = {
        "status": "EXECUTION_COMPLETE_PENDING_VISUAL_REVIEW",
        "started_utc": started.isoformat(),
        "completed_utc": completed.isoformat(),
        "elapsed_seconds": (completed - started).total_seconds(),
        "command": f"{sys.executable} {Path(__file__).resolve()}",
        "host": platform.node(),
        "platform": platform.platform(),
        "python": sys.version,
        "cairosvg": cairosvg.__version__,
        "operation": "clean standalone vector rebuild from reviewed independent panel artwork",
        "scientific_values_recalculated": False,
        "existing_outputs_overwritten": False,
        "panels": panel_records,
        "contact_sheet": str(contact_sheet),
    }
    (PROV / "execution_record.json").write_text(json.dumps(record, indent=2) + "\n", encoding="utf-8")
    with (PROV / "input_manifest.sha256.tsv").open("w", encoding="utf-8") as handle:
        handle.write("sha256\tpath\n")
        for digest, path in input_manifest:
            handle.write(f"{digest}\t{path}\n")
    with (PROV / "output_manifest.sha256.tsv").open("w", encoding="utf-8") as handle:
        handle.write("sha256\tpath\n")
        for digest, path in output_manifest:
            handle.write(f"{digest}\t{path}\n")
    print(json.dumps(record, indent=2))


if __name__ == "__main__":
    main()
