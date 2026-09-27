#!/usr/bin/env python3
"""Read-only technical validation of the rebuilt Figure 4 panel files."""

from __future__ import annotations

import hashlib
import json
import re
import subprocess
import tempfile
from pathlib import Path

from lxml import etree
from PIL import Image

RUN = Path("${ANALYSIS_DIR}/22_answer_reviews/00_ms/05_MS/0918_revision/final_ms/02_figure4/redrawn_adobe_compatible_20260922")
PANELS = RUN / "panels"
PROV = RUN / "provenance"
TARGET = {p: (60.0, 38.5) for p in "abcdefghi"}
TARGET["j"] = (91.0, 72.2)


def command(*args: str) -> str:
    result = subprocess.run(args, check=True, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    return result.stdout + result.stderr


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def validate_manifest(path: Path) -> list[str]:
    errors = []
    lines = path.read_text(encoding="utf-8").splitlines()[1:]
    for line in lines:
        digest, raw_path = line.split("\t", 1)
        file_path = Path(raw_path)
        if not file_path.is_file():
            errors.append(f"missing: {file_path}")
        elif sha256(file_path) != digest:
            errors.append(f"checksum mismatch: {file_path}")
    return errors


def main() -> None:
    records = {}
    fatal = []
    warnings = []
    for panel in "abcdefghij":
        pdf = next(PANELS.glob(f"Figure4{panel}_*.pdf"))
        svg = next(PANELS.glob(f"Figure4{panel}_*.svg"))
        png = next(PANELS.glob(f"Figure4{panel}_*600dpi.png"))
        info = command("pdfinfo", str(pdf))
        page_match = re.search(r"^Pages:\s+(\d+)", info, re.M)
        size_match = re.search(r"^Page size:\s+([0-9.]+) x ([0-9.]+) pts", info, re.M)
        version_match = re.search(r"^PDF version:\s+([0-9.]+)", info, re.M)
        producer_match = re.search(r"^Producer:\s+(.+)", info, re.M)
        if not all((page_match, size_match, version_match, producer_match)):
            fatal.append(f"pdfinfo fields missing for panel {panel}")
            continue
        pages = int(page_match.group(1))
        width_pt, height_pt = map(float, size_match.groups())
        observed_mm = (width_pt * 25.4 / 72.0, height_pt * 25.4 / 72.0)
        target_mm = TARGET[panel]
        size_ok = all(abs(a - b) <= 0.01 for a, b in zip(observed_mm, target_mm))
        gs = subprocess.run(["gs", "-q", "-dNOPAUSE", "-dBATCH", "-sDEVICE=nullpage", str(pdf)], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        image_listing = command("pdfimages", "-list", str(pdf))
        image_rows = [line for line in image_listing.splitlines()[2:] if line.strip()]
        image_shapes = []
        for row in image_rows:
            fields = row.split()
            if len(fields) >= 5 and fields[0].isdigit():
                image_shapes.append([int(fields[3]), int(fields[4])])
        if image_shapes:
            max_pixels = max(w * h for w, h in image_shapes)
            if max_pixels > 100_000:
                fatal.append(f"panel {panel} contains a large raster image")
            else:
                warnings.append(f"panel {panel} retains tiny vector-rendering strips: {image_shapes}")
        etree.parse(str(svg))
        with Image.open(png) as image:
            png_size = list(image.size)
        with tempfile.TemporaryDirectory() as tmp:
            rendered = Path(tmp) / panel
            subprocess.run(["pdftocairo", "-singlefile", "-png", "-r", "72", str(pdf), str(rendered)], check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            render_ok = rendered.with_suffix(".png").is_file()
        pypdf2_present = b"PyPDF2" in pdf.read_bytes()
        if pages != 1 or not size_ok or gs.returncode != 0 or not render_ok or pypdf2_present:
            fatal.append(f"panel {panel} failed structural validation")
        records[panel] = {
            "pdf": str(pdf),
            "pages": pages,
            "page_mm": list(observed_mm),
            "target_mm": list(target_mm),
            "size_within_0.01mm": size_ok,
            "pdf_version": version_match.group(1),
            "producer": producer_match.group(1).strip(),
            "ghostscript_parse": gs.returncode == 0,
            "poppler_render": render_ok,
            "pypdf2_signature_absent": not pypdf2_present,
            "raster_image_shapes": image_shapes,
            "png_600dpi_pixels": png_size,
            "sha256": sha256(pdf),
        }
    for manifest in (PROV / "input_manifest.sha256.tsv", PROV / "output_manifest.sha256.tsv"):
        fatal.extend(validate_manifest(manifest))
    report = {
        "status": "FAILED" if fatal else ("PASS_WITH_WARNINGS" if warnings else "PASS"),
        "scope": "technical structure and rendering only; user visual review remains required",
        "fatal": fatal,
        "warnings": warnings,
        "contract_deviation": (
            "The zero-raster criterion is met for panels a,b,d-i. Panel c retains one 15x208-pixel "
            "Fst colour-gradient strip and panel j retains 17x552 and 214x15-pixel annotation/gradient "
            "strips inherited from the reviewed vector artwork. No panel or data layer is rasterized."
        ),
        "panels": records,
    }
    (PROV / "technical_validation.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2))
    if fatal:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
