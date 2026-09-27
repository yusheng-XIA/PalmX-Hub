#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Rebuild Extended_Data_Fig_01 with the latest high-resolution originals for panels c and d.

- Panels a (study design), b (hap circos) and e (Hi-C map) stay byte-identical to
  the original page raster (unchanged pixels).
- Panel c (HiFi/ONT coverage) and panel d (IGV gap filling) are replaced by the
  full-resolution JPEG originals extracted from Nature-supplementary information.doc.
- Page size and the a-e label positions are kept exactly as in the original.
"""

from __future__ import annotations

import argparse
import hashlib
from pathlib import Path

import fitz
import numpy as np
from PIL import Image, ImageDraw, ImageFilter

Image.MAX_IMAGE_PIXELS = None

LABELS = (
    ("a", (2.6, 10.1, 8.7, 20.6), (2.2, 13.9, 7.4, 19.7)),
    ("b", (222.0, 10.1, 228.4, 20.6), (221.5, 12.2, 226.3, 19.7)),
    ("c", (2.9, 203.6, 8.2, 214.0), (2.4, 207.4, 7.4, 213.1)),
    ("d", (222.0, 203.6, 228.4, 214.0), (221.5, 205.9, 226.3, 213.1)),
    ("e", (2.6, 346.4, 8.7, 356.8), (2.2, 350.2, 7.4, 355.9)),
)
LABEL_FS = 9.0
LABEL_DESCENDER = 0.236
RASTER_RECT = (-28.3, -27.1, 508.8, 511.0)
WHITEOUT = ((1.0, 206.0, 226.0, 349.0), (227.0, 206.0, 482.88, 497.5))
WHITEOUT_MOVE_E = ((1.0, 206.0, 226.0, 497.5), (227.0, 206.0, 482.88, 497.5))
NEW_PANELS = {
    "c": (10.0, 215.0, 225.0, 349.0),
    "d": (227.0, 208.0, 482.88, 497.5),
}
E_SRC_RECT = (2.5, 349.5, 186.6, 496.3)
E_DEST_SLOT = (262.89, 215.0, 446.99, 400.0)
E_LABEL_ERASE = (1.0, 348.5, 9.6, 359.5)
A_SLOT = (2.88, 14.64, 210.48, 195.36)
B_SLOT = (222.2, 3.12, 452.36, 205.92)
B_CROPBOX = (25.9, 36.2, 5952.5, 5967.0)
WHITEOUT_FULL = ((0.0, 0.0, 482.88, 346.5), (227.0, 346.5, 482.88, 497.5))


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def fit_rect(slot, aspect):
    x0, y0, x1, y1 = slot
    sw, sh = x1 - x0, y1 - y0
    if aspect >= sw / sh:
        w, h = sw, sw / aspect
    else:
        h, w = sh, sh * aspect
    cx = x0 + (sw - w) / 2
    return (cx, y0, cx + w, y0 + h)


def fix_fillled_typo(doc, font_out: Path) -> bool:
    """Fix the 'Fillled gaps' typo in the vector coverage PDF (in-memory copy only)."""
    page = doc[0]
    hits = page.search_for("Fillled gaps")
    if not hits:
        return False
    for rect in hits:
        page.add_redact_annot(rect)
    page.apply_redactions(images=fitz.PDF_REDACT_IMAGE_NONE)
    fontfile = None
    for entry in page.get_fonts(full=True):
        xref, ext, ftype, basefont, name, enc = entry[:6]
        if "Arial" in basefont:
            _, fext, _, buf = doc.extract_font(xref)
            if buf:
                fontfile = font_out.with_suffix(f".{fext}")
                fontfile.write_bytes(buf)
            break
    origin, kwargs = (976.32, 657.15), dict(fontsize=11, color=(0, 0, 0))
    if fontfile is not None:
        page.insert_text(origin, "Filled gaps", fontname="arial",
                         fontfile=str(fontfile), **kwargs)
    else:
        page.insert_text(origin, "Filled gaps", fontname="helv", **kwargs)
    return True


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--original", required=True, type=Path)
    ap.add_argument("--new-c", type=Path,
                    help="raster replacement for panel c (not used with --c-pdf/--d-pdf)")
    ap.add_argument("--new-d", type=Path,
                    help="replacement for panel d (not used with --drop-d/--move-e-as-d/--c-pdf)")
    ap.add_argument("--c-pdf", type=Path,
                    help="vector PDF whose page 0 replaces panel c (use together with --d-pdf)")
    ap.add_argument("--d-pdf", type=Path,
                    help="PDF whose page 0 replaces panel d (use together with --c-pdf)")
    ap.add_argument("--a-png", type=Path,
                    help="high-resolution PNG that replaces panel a (use with --b-pdf/--c-pdf/--d-pdf)")
    ap.add_argument("--b-pdf", type=Path,
                    help="vector PDF whose page 0 replaces panel b (use with --a-png/--c-pdf/--d-pdf)")
    ap.add_argument("--fix-typo", action="store_true",
                    help="fix the 'Fillled gaps' typo in the embedded c PDF copy")
    ap.add_argument("--out-dir", required=True, type=Path)
    ap.add_argument("--font", type=Path,
                    default=Path("/usr/share/fonts/dejavu/DejaVuSans-Bold.ttf"))
    ap.add_argument("--dpi", type=int, default=300)
    ap.add_argument("--bg-upscale", type=int, default=1,
                    help="integer upscale factor for the original page raster (a/b/e) so the "
                         "embedded resolution matches the target render dpi (2 => ~600 dpi)")
    ap.add_argument("--bg-sharpen", action="store_true",
                    help="apply a mild unsharp mask after --bg-upscale")
    ap.add_argument("--compact", type=float, default=0.0, metavar="PAD_PT",
                    help="crop the page to the content bounding box plus PAD_PT margin (pt)")
    ap.add_argument("--drop-d", action="store_true",
                    help="omit panel d (right column below b); label d is removed as well")
    ap.add_argument("--move-e-as-d", action="store_true",
                    help="crop panel e from the original raster, move it under b and relabel it d "
                         "(a|b over c|e two-row layout; original d dropped)")
    ap.add_argument("--out-stem", default="Extended_Data_Fig_01_rebuilt_20260921")
    ap.add_argument("--record", default="REBUILD_RECORD_20260921.tsv")
    args = ap.parse_args()
    replace_cd = args.c_pdf is not None and args.d_pdf is not None
    replace_abcd = replace_cd and args.a_png is not None and args.b_pdf is not None
    if (args.c_pdf is None) != (args.d_pdf is None):
        ap.error("--c-pdf and --d-pdf must be given together")
    if (args.a_png is None) != (args.b_pdf is None):
        ap.error("--a-png and --b-pdf must be given together")
    if (args.a_png is not None or args.b_pdf is not None) and not replace_cd:
        ap.error("--a-png/--b-pdf need --c-pdf and --d-pdf as well (full a/b/c/d replacement)")
    if replace_cd and (args.drop_d or args.move_e_as_d):
        ap.error("--c-pdf/--d-pdf cannot be combined with --drop-d/--move-e-as-d")
    if args.drop_d and args.move_e_as_d:
        ap.error("--drop-d and --move-e-as-d are mutually exclusive")
    if not replace_cd and args.new_c is None:
        ap.error("--new-c is required unless --c-pdf/--d-pdf is used")
    if not replace_cd and not args.drop_d and not args.move_e_as_d and args.new_d is None:
        ap.error("--new-d is required unless --drop-d, --move-e-as-d or --c-pdf/--d-pdf is used")

    src = fitz.open(args.original)
    page = src[0]
    page_size = (page.rect.width, page.rect.height)
    infos = page.get_image_info(xrefs=True)
    if len(infos) != 1:
        raise RuntimeError(f"expected exactly 1 raster on original page, got {len(infos)}")
    info = infos[0]
    bbox = info["bbox"]
    raw = src.extract_image(info["xref"])
    raster = Image.open(__import__("io").BytesIO(raw["image"])).convert("RGB")
    if raster.size != (info["width"], info["height"]):
        raise RuntimeError("extracted raster size mismatch")
    sx = info["width"] / (bbox[2] - bbox[0])
    sy = info["height"] / (bbox[3] - bbox[1])
    if args.bg_upscale < 1:
        raise SystemExit("--bg-upscale must be >= 1")
    if args.bg_upscale > 1:
        raster = raster.resize(
            (raster.width * args.bg_upscale, raster.height * args.bg_upscale), Image.LANCZOS)
        if args.bg_sharpen:
            raster = raster.filter(ImageFilter.UnsharpMask(radius=1.2, percent=60, threshold=3))
    f = float(args.bg_upscale)

    def to_px(x, y):
        return (x - bbox[0]) * sx * f, (y - bbox[1]) * sy * f

    def px_to_pt(px, py):
        return bbox[0] + px / (sx * f), bbox[1] + py / (sy * f)

    args.out_dir.mkdir(parents=True, exist_ok=True)
    eff_dpi = 300 * args.bg_upscale
    e_crop = None
    if args.move_e_as_d:
        ex0, ey0, ex1, ey1 = E_SRC_RECT
        e_crop = args.out_dir / f"_ED_Fig_01_panel_e_crop_{eff_dpi}dpi.png"
        crop = raster.crop((*to_px(ex0, ey0), *to_px(ex1, ey1)))
        crop_draw = ImageDraw.Draw(crop)
        lx0, ly0 = to_px(E_LABEL_ERASE[0], E_LABEL_ERASE[1])
        lx1, ly1 = to_px(E_LABEL_ERASE[2], E_LABEL_ERASE[3])
        crop_draw.rectangle((lx0 - to_px(ex0, ey0)[0], ly0 - to_px(ex0, ey0)[1],
                             lx1 - to_px(ex0, ey0)[0], ly1 - to_px(ex0, ey0)[1]), fill="white")
        crop.save(e_crop)

    draw = ImageDraw.Draw(raster)
    whiteout = (WHITEOUT_FULL if replace_abcd else
                WHITEOUT_MOVE_E if args.move_e_as_d else WHITEOUT)
    for x0, y0, x1, y1 in whiteout:
        px0, py0 = to_px(x0, y0)
        px1, py1 = to_px(x1, y1)
        draw.rectangle((px0, py0, px1, py1), fill="white")
    bg = args.out_dir / f"_ED_Fig_01_background_whiteout_{eff_dpi}dpi.png"
    raster.save(bg)

    placements = {}
    pdf_pages = {}
    a_image = None
    if replace_cd:
        cdoc = fitz.open(args.c_pdf)
        ddoc = fitz.open(args.d_pdf)
        if args.fix_typo and fix_fillled_typo(cdoc, args.out_dir / "_gaps_readscov_ArialMT"):
            print("typo fixed: 'Fillled gaps' -> 'Filled gaps' in embedded c copy")
        for key, doc in (("c", cdoc), ("d", ddoc)):
            r = doc[0].rect
            placements[key] = fit_rect(NEW_PANELS[key], r.width / r.height)
        pdf_pages = {"c": cdoc, "d": ddoc}
        if replace_abcd:
            a_src = Image.open(args.a_png).convert("RGB")
            g = np.asarray(a_src.convert("L"))
            m = g < 245
            cols, rows = m.sum(axis=0) >= 3, m.sum(axis=1) >= 3
            ax0, ax1 = int(np.argmax(cols)), len(cols) - int(np.argmax(cols[::-1]))
            ay0, ay1 = int(np.argmax(rows)), len(rows) - int(np.argmax(rows[::-1]))
            a_image = args.out_dir / "_ED_Fig_01_panel_a_crop.png"
            a_src.crop((ax0, ay0, ax1, ay1)).save(a_image)
            placements["a"] = fit_rect(A_SLOT, (ax1 - ax0) / (ay1 - ay0))
            bdoc = fitz.open(args.b_pdf)
            bdoc[0].set_cropbox(fitz.Rect(*B_CROPBOX))
            br = bdoc[0].rect
            placements["b"] = fit_rect(B_SLOT, br.width / br.height)
            pdf_pages["b"] = bdoc
        labels = list(LABELS)
    else:
        slots = dict(NEW_PANELS)
        if args.move_e_as_d:
            slots["d"] = E_DEST_SLOT
            sources = (("c", args.new_c), ("d", e_crop))
        elif args.drop_d:
            sources = (("c", args.new_c),)
        else:
            sources = (("c", args.new_c), ("d", args.new_d))
        for key, path in sources:
            im = Image.open(path)
            placements[key] = fit_rect(slots[key], im.width / im.height)
        labels = []
        for letter, text_bbox, box in LABELS:
            if args.drop_d and letter == "d":
                continue
            if args.move_e_as_d and letter == "e":
                continue
            labels.append((letter, text_bbox, box))

    offset = (0.0, 0.0)
    crop_rect = None
    if args.compact > 0:
        gray = np.asarray(raster.convert("L"))
        mask = gray < 245
        cols = mask.sum(axis=0) >= 3
        rows = mask.sum(axis=1) >= 3
        if not cols.any() or not rows.any():
            raise RuntimeError("background raster has no content to crop to")
        rx0, rx1 = int(np.argmax(cols)), len(cols) - int(np.argmax(cols[::-1]))
        ry0, ry1 = int(np.argmax(rows)), len(rows) - int(np.argmax(rows[::-1]))
        x0, y0 = px_to_pt(rx0, ry0)
        x1, y1 = px_to_pt(rx1, ry1)
        for rect in placements.values():
            x0, y0 = min(x0, rect[0]), min(y0, rect[1])
            x1, y1 = max(x1, rect[2]), max(y1, rect[3])
        for _, text_bbox, box in labels:
            x0, y0 = min(x0, text_bbox[0], box[0]), min(y0, text_bbox[1], box[1])
            x1, y1 = max(x1, text_bbox[2], box[2]), max(y1, text_bbox[3], box[3])
        pad = args.compact
        crop_rect = (x0 - pad, y0 - pad, x1 + pad, y1 + pad)
        offset = (crop_rect[0], crop_rect[1])
        new_size = (crop_rect[2] - crop_rect[0], crop_rect[3] - crop_rect[1])
    else:
        new_size = page_size
    ox, oy = offset

    def shift(rect):
        return fitz.Rect(rect[0] - ox, rect[1] - oy, rect[2] - ox, rect[3] - oy)

    out = fitz.open()
    new_page = out.new_page(width=new_size[0], height=new_size[1])
    new_page.insert_image(shift(RASTER_RECT), filename=str(bg), keep_proportion=True)
    if replace_cd:
        if replace_abcd:
            new_page.insert_image(shift(placements["a"]), filename=str(a_image),
                                  keep_proportion=True)
            new_page.show_pdf_page(shift(placements["b"]), pdf_pages["b"], 0)
        new_page.show_pdf_page(shift(placements["c"]), pdf_pages["c"], 0)
        new_page.show_pdf_page(shift(placements["d"]), pdf_pages["d"], 0)
    else:
        for key, rect in placements.items():
            new_page.insert_image(shift(rect), filename=str(sources[0 if key == "c" else 1][1]),
                                  keep_proportion=True)
    new_page.insert_font(fontname="dejavu", fontfile=str(args.font))
    for letter, text_bbox, box in labels:
        new_page.draw_rect(shift(box), color=None, fill=(1, 1, 1))
        baseline = text_bbox[3] - LABEL_DESCENDER * LABEL_FS
        new_page.insert_text((text_bbox[0] - ox, baseline - oy), letter, fontsize=LABEL_FS,
                             fontname="dejavu", color=(0, 0, 0))
    placements = {k: fitz.Rect(v[0] - ox, v[1] - oy, v[2] - ox, v[3] - oy)
                  for k, v in placements.items()}
    if replace_abcd:
        panel_desc = ("a/b/c/d replaced (a: high-res PNG, b/c: vector PDFs, d: IGV montage PDF); "
                      "e kept from original raster")
    elif replace_cd:
        panel_desc = "a/b/c/d/e (c,d replaced from PDFs; e kept from original raster)"
    elif args.move_e_as_d:
        panel_desc = "a/b/c/d (e moved under b, relabelled d, d-slot image dropped)"
    elif args.drop_d:
        panel_desc = "a/b/c/e (d dropped)"
    else:
        panel_desc = "a/b/c/d/e"
    subject = (f"a <- {args.a_png.name}, b <- {args.b_pdf.name} (vector), c <- {args.c_pdf.name} "
               f"(vector), d <- {args.d_pdf.name}; e from the 300 dpi page raster "
               f"({eff_dpi} dpi equivalent after upscale); page {new_size[0]:.1f}x{new_size[1]:.1f}pt"
               if replace_abcd else
               f"c <- {args.c_pdf.name} (vector), d <- {args.d_pdf.name}; a/b/e from the "
               f"300 dpi page raster ({eff_dpi} dpi equivalent after upscale); "
               f"page {new_size[0]:.1f}x{new_size[1]:.1f}pt" if replace_cd else
               f"Panels a/b/e from original raster ({eff_dpi} dpi equivalent after upscale); "
               f"c replaced with the 1688 dpi original; compact page {new_size[0]:.1f}x{new_size[1]:.1f}pt")
    out.set_metadata({
        "title": f"Extended Data Fig 1 (rebuilt, panels {panel_desc})",
        "subject": subject,
        "creator": "PyMuPDF rebuild 2026-09-21",
    })
    pdf = args.out_dir / f"{args.out_stem}.pdf"
    png = args.out_dir / f"{args.out_stem}.png"
    out.save(pdf, garbage=4, deflate=True)
    out.close()
    src.close()

    check = fitz.open(pdf)
    text = check[0].get_text()
    if replace_cd:
        expect_letters, expect_images = "abcde", None
    elif args.move_e_as_d:
        expect_letters, expect_images = "abcd", 3
    elif args.drop_d:
        expect_letters, expect_images = "abce", 2
    else:
        expect_letters, expect_images = "abcde", 3
    letters = [l for l in expect_letters if l in text]
    images = len(check[0].get_images())
    if letters != list(expect_letters) or (expect_images is not None and images != expect_images):
        raise RuntimeError(f"post-check failed images={images} letters={letters}")
    check[0].get_pixmap(dpi=args.dpi, alpha=False).save(png)
    check.close()

    with (args.out_dir / args.record).open("w", encoding="utf-8") as fh:
        extra = []
        if replace_cd:
            c_src, c_sha = str(args.c_pdf), sha256(args.c_pdf)
            d_src, d_sha = str(args.d_pdf), sha256(args.d_pdf)
            c_dpi = "vector (infinite)"
            note = ("c and d embedded as vector/PDF form XObjects; a/b/e Lanczos-upscaled from the "
                    "300 dpi page raster (no new detail); 'Fillled gaps' typo fixed in the embedded "
                    "c copy" if args.fix_typo else
                    "c and d embedded as vector/PDF form XObjects; a/b/e Lanczos-upscaled from the "
                    "300 dpi page raster (no new detail)")
            if replace_abcd:
                a_w = placements["a"].width if hasattr(placements["a"], "width") else 0
                extra = [
                    ("new_panel_a", f"{args.a_png} (cropped to content)"),
                    ("new_panel_a_sha256", sha256(args.a_png)),
                    ("new_panel_b", str(args.b_pdf)),
                    ("new_panel_b_sha256", sha256(args.b_pdf)),
                    ("panel_a_rect_pt", " ".join(f"{v:.2f}" for v in placements["a"])),
                    ("panel_b_rect_pt", " ".join(f"{v:.2f}" for v in placements["b"])),
                    ("panel_b_cropbox_pt", " ".join(f"{v:.1f}" for v in B_CROPBOX)),
                ]
                note = ("a from high-res PNG (cropped to content); b and c embedded as vector PDF "
                        "form XObjects; d embedded as PDF (IGV screenshots at ~1600 dpi effective); "
                        "e upscaled from the 300 dpi page raster (no new detail); 'Fillled gaps' typo "
                        "fixed in the embedded c copy")
        elif args.drop_d:
            c_src, c_sha = str(args.new_c), sha256(args.new_c)
            d_src, d_sha = "dropped", "none"
            c_dpi = f"{Image.open(args.new_c).width / (placements['c'].width / 72):.0f}"
            note = "panels a/b/e unchanged pixels from the 300 dpi page raster"
        elif args.move_e_as_d:
            c_src, c_sha = str(args.new_c), sha256(args.new_c)
            d_src, d_sha = f"{e_crop} (panel e cropped from original raster)", sha256(e_crop)
            c_dpi = f"{Image.open(args.new_c).width / (placements['c'].width / 72):.0f}"
            note = ("panels a/b and the moved e(d) are Lanczos-upscaled from the 300 dpi page "
                    "raster (no new detail); panel c is native ~1688 dpi")
        else:
            c_src, c_sha = str(args.new_c), sha256(args.new_c)
            d_src, d_sha = str(args.new_d), sha256(args.new_d)
            c_dpi = f"{Image.open(args.new_c).width / (placements['c'].width / 72):.0f}"
            note = ("panels a/b/e unchanged pixels from the 300 dpi page raster; panel c is the "
                    "~1688 dpi original")
        untouched = ("panels a/b unchanged; panel e kept in place from original raster"
                     if replace_cd else
                     ("panels a/b unchanged; panel e cropped from original raster and placed under b"
                      if args.move_e_as_d else
                      "panels a/b/e unchanged pixels from original page raster"))
        fh.write("item\tvalue\n")
        for key, val in (
            ("status", "passed"),
            ("original_pdf", str(args.original)),
            ("original_pdf_sha256", sha256(args.original)),
            ("panels_in_output", panel_desc),
            ("new_panel_c", c_src),
            ("new_panel_c_sha256", c_sha),
            ("new_panel_d", d_src),
            ("new_panel_d_sha256", d_sha),
            ("panel_e_source_rect_pt", " ".join(f"{v:.2f}" for v in E_SRC_RECT)
             if args.move_e_as_d else "n/a"),
            ("page_original_pt", f"{page_size[0]}x{page_size[1]}"),
            ("page_output_pt", f"{new_size[0]:.2f}x{new_size[1]:.2f}"),
            ("crop_rect_pt", " ".join(f"{v:.2f}" for v in crop_rect) if crop_rect else "none"),
            ("bg_upscale", f"{args.bg_upscale}x (Lanczos{', unsharp' if args.bg_sharpen else ''})"),
            ("panel_a_b_effective_dpi", str(eff_dpi)),
            ("panel_c_effective_dpi", c_dpi),
            ("render_dpi_png", str(args.dpi)),
            ("panel_c_rect_pt", " ".join(f"{v:.2f}" for v in placements["c"])),
            ("panel_d_rect_pt", " ".join(f"{v:.2f}" for v in placements["d"])
             if "d" in placements else "dropped"),
            ("panels_a_b_e", untouched),
            ("note", note),
            ("output_pdf", str(pdf)),
            ("output_pdf_sha256", sha256(pdf)),
            ("output_png", str(png)),
            ("output_png_sha256", sha256(png)),
        ):
            fh.write(f"{key}\t{val}\n")
    print(f"PASS {pdf}")
    print(f"PASS {png}")
    print(f"panels: {panel_desc}; page {new_size[0]:.1f}x{new_size[1]:.1f}pt; "
          f"c rect {placements['c']}"
          + (f" | d/e rect {placements['d']}" if "d" in placements else ""))


if __name__ == "__main__":
    main()
