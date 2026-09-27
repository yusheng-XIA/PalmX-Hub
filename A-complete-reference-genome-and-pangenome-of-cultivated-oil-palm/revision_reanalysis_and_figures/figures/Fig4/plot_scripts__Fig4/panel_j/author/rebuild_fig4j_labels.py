#!/usr/bin/env python3
"""Reframe the 33-material RGA panel with current display names (no data recalculation)."""
import copy
import csv
import hashlib
import importlib.util
import json
import platform
import re
import sys
from datetime import datetime, timezone
from pathlib import Path

import cairosvg
import fitz
from lxml import etree

RUN = Path(__file__).resolve().parents[1]
FINAL = RUN.parents[1]
BASE = FINAL.parents[3] / "03_V3/04_figure4"
SOURCE = BASE / "RGA_BGC_meizhou4_33varieties_20260805"
SVG = SOURCE / "figures/RGA_tree_dotmatrix_major_groups.svg"
TABLE = SOURCE / "tables/RGA_summary_33varieties.csv"
TIPS = SOURCE / "tables/rga_plot/RGA_tree_dotmatrix_tip_order.tsv"
REBUILD = FINAL / "scripts/rebuild_figure4_adobe_compatible_20260922.py"
OUT = RUN / "outputs"


def digest(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def main():
    with TABLE.open(newline="") as f:
        rows = list(csv.DictReader(f))
    with TIPS.open(newline="") as f:
        tips = list(csv.DictReader(f, delimiter="\t"))
    assert len(rows) == len(tips) == 33
    materials = {str(r["Genome"]) for r in rows}
    assert len(materials) == 33
    assert {str(r["Genome"]) for r in tips} == materials
    assert {"FL", "BK", "meizhou4"} <= materials
    assert not ({"bk_hap1", "bk_hap2", "Africa_hap2", "American_hap1"} & materials)
    assert [int(t["Y_order_top_to_bottom"]) for t in tips] == list(range(1, 34))
    svg_tree = etree.parse(str(SVG), etree.XMLParser(resolve_entities=False, no_network=True))
    comments = {c.text.strip() for c in svg_tree.xpath("//comment()") if c.text}
    assert {t["Display"] for t in tips} <= comments
    assert not {"bk_hap1", "bk_hap2", "Africa_hap2", "American_hap1"} & comments
    spec = importlib.util.spec_from_file_location("standalone_figure4", REBUILD)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    output_svg = OUT / "Figure4j_33materials_current_names_91x72.2mm_vector.svg"
    output_pdf = OUT / "Figure4j_91mm_width_proportional_vector.pdf"
    output_png = RUN / "qc/Figure4j_33materials_600dpi.png"
    # The old assembler mistakenly removed SVG IDs text_54/55/56, which in
    # this 33-material SVG are the abundance scale, group title and FL legend.
    # The obsolete header is removed solely by the 32-unit title-band crop.
    module.remove_j_heading = lambda root: None
    geometry = module.build_outer_svg("j", SVG, output_svg)
    panel_tree = etree.parse(str(output_svg))
    svg_ns = "http://www.w3.org/2000/svg"
    for element_id, label, y in (("text_56", "FL (merged material)", 536.69175),
                                  ("text_57", "TN (merged material)", 548.43425)):
        matches = panel_tree.xpath(f".//*[@id='{element_id}']")
        assert len(matches) == 1
        old = matches[0]
        parent = old.getparent()
        at = parent.index(old)
        parent.remove(old)
        node = etree.Element(f"{{{svg_ns}}}text", x="62.648", y=str(y),
                             style="font-family:'DejaVu Sans',sans-serif;font-size:8px;fill:#262626")
        node.text = label
        parent.insert(at, node)
    # Only rename terminal display glyphs; preserve the source tree, order,
    # counts, dot encoding and immutable original tables (which use old IDs).
    # Matplotlib places shared font glyph definitions inside some label groups.
    # Hoist them before removing the groups or unrelated labels lose letters.
    # Keep glyphs in the source nested SVG's own definitions scope.
    nested_source = panel_tree.xpath(".//*[@id='figure_1']")[0].getparent()
    root_defs = nested_source.find(f"{{{svg_ns}}}defs")
    assert root_defs is not None
    names = {"text_13": ("EG_dura", "TK"),
             "text_14": ("EG_pisifera", "NS"),
             "text_15": ("BK", "TN"),
             "text_27": ("Meizhou4", "E.Oleifera")}
    for element_id, (old_name, new_name) in names.items():
        matches = panel_tree.xpath(f".//*[@id='{element_id}']")
        assert len(matches) == 1
        old = matches[0]
        assert len(old.xpath("./comment()")) == 1
        assert old.xpath("./comment()")[0].text.strip() == old_name
        transform = old.xpath("./*[local-name()='g']")[0].get("transform")
        match = re.fullmatch(r"translate\(([-\d.]+) ([-\d.]+)\) scale\(0\.072 -0\.072\)", transform)
        assert match, (element_id, transform)
        for glyph in old.xpath(".//*[local-name()='defs']/*"):
            root_defs.append(copy.deepcopy(glyph))
        parent = old.getparent()
        at = parent.index(old)
        parent.remove(old)
        node = etree.Element(f"{{{svg_ns}}}text", x="370", y=match.group(2),
                             style="font-family:'DejaVu Sans',sans-serif;font-size:7.2px;fill:#222222;text-anchor:end")
        node.text = new_name
        parent.insert(at, node)
    panel_tree.write(str(output_svg), encoding="utf-8", xml_declaration=True)
    cairosvg.svg2pdf(url=str(output_svg), write_to=str(output_pdf))
    cairosvg.svg2png(url=str(output_svg), write_to=str(output_png),
                     output_width=round(91 / 25.4 * 600),
                     output_height=round(72.2 / 25.4 * 600))
    doc = fitz.open(output_pdf)
    assert len(doc) == 1
    assert abs(doc[0].rect.width - 91 / 25.4 * 72) < .1
    assert abs(doc[0].rect.height - 72.2 / 25.4 * 72) < .1
    assert len(doc[0].get_drawings()) > 100
    record = {
        "status": "RENDERED_PENDING_REVIEW_AND_PROMOTION",
        "utc": datetime.now(timezone.utc).isoformat(),
        "command": f"{sys.executable} {Path(__file__).resolve()}",
        "host": platform.node(), "python": sys.version,
        "cairosvg": cairosvg.__version__, "pymupdf": fitz.VersionBind,
        "source_files_sha256": {str(p): digest(p) for p in (SVG, TABLE, TIPS, REBUILD)},
        "outputs_sha256": {str(p): digest(p) for p in (output_svg, output_pdf, output_png)},
        "materials": len(materials), "plotted_tips": len(tips),
        "source_material_ids_unchanged": True,
        "display_renames": {old: new for old, new in names.values()},
        "merged_display_tips": ["FL", "TN", "E.Oleifera"],
        "geometry": geometry,
        "caveat": "E.Oleifera (source ID meizhou4) is a display-only tree graft, not a new phylogenetic inference.",
    }
    (RUN / "provenance/render_record.json").write_text(json.dumps(record, indent=2) + "\n")
    print(json.dumps(record, indent=2))


if __name__ == "__main__":
    main()
