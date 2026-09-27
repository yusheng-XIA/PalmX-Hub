#!/usr/bin/env python3
"""Figure 5: vector text corrections and height reduction.

Step 1 (src.pdf -> edited.pdf): delete selected text runs with text-only redactions
(images and line art untouched) and rewrite them at the original baseline, size and
colour with Arial / Arial Bold. Step 2 (edited.pdf -> Figure5.pdf): keep panel a (row 1,
4.5-pt labels) at 100% and place rows b–i, clipped from the same vector page, at 91.2%
so that the page fits 183 × 247 mm with all rescaled text >= 5 pt.
"""
import re
from pathlib import Path

import fitz

HERE = Path(__file__).resolve().parent
ARIAL = str(HERE / "ArialFix.ttf")  # Arial with U+00AD/U+00A0 dropped from cmap, so "-" and " " extract correctly
ARIAL_B = str(HERE / "ArialFixBold.ttf")
OUT = HERE.parent / "deliver/Main_Figures_revised"
OUT.mkdir(parents=True, exist_ok=True)

HAP = {"Dura": "TK", "dura": "TK", "Pisifera": "NS", "pisifera": "NS", "Nigerian": "Nigerian", "TN": "TN"}


def new_text(t, x0, y0):
    """Return replacement text, '' to delete, or None to keep."""
    m = re.fullmatch(r"(Dura|Pisifera|Nigerian|TN)_hap([12])", t)
    if m:                                            # 5h legend and 5i rows (ST24 naming)
        return f"{HAP[m.group(1)]}-Hap{m.group(2)}"
    m = re.fullmatch(r"EG-(\d+)", t)
    if m:                                            # 5a rows -> EG_008 form used in ST18/22/24
        return f"EG_{int(m.group(1)):03d}"
    if t == "EG-houke":
        return "EG_houke"
    if t == "Donor haplotype (n=35)":
        return "Donor haplotype (n = 35)"
    if t == "P < 2.2e-16":
        return "P < 2.2 × 10^−16"
    if t == "Number of SVs (x 10³)":
        return "Number of SVs (×10³)"
    if 179 < y0 < 182 and t == "0" and x0 > 40:     # 5a x axis: keep only chr01's "0"
        return ""
    if 179 < y0 < 182 and t == "0 25":
        return "25"
    return None


def main():
    doc = fitz.open(HERE / "src.pdf")
    page = doc[0]
    edits = []
    for b in page.get_text("rawdict")["blocks"]:
        for line in b.get("lines", []):
            for s in line["spans"]:
                t = "".join(c["c"] for c in s["chars"])
                nt = new_text(t, s["bbox"][0], s["bbox"][1])
                if nt is None or nt == t:
                    continue
                if t == "0 25":                      # keep the x position of "25"
                    two = [c for c in s["chars"] if c["c"] == "2"][0]
                    origin = (two["origin"][0], s["origin"][1])
                else:
                    origin = s["origin"]
                cen = [((c["bbox"][0] + c["bbox"][2]) / 2, c["origin"][1] - s["size"] * 0.33)
                       if line["dir"][0] == 1 else (c["origin"][0] - s["size"] * 0.33, (c["bbox"][1] + c["bbox"][3]) / 2)
                       for c in s["chars"] if c["c"].strip()]
                if t == "0 25":
                    cen = cen[:1]                    # only the leading "0" goes
                edits.append(dict(old=t, new=nt, cen=cen, bbox=fitz.Rect(s["bbox"]), origin=origin, size=s["size"],
                                  bold="Bold" in s["font"], color=s["color"], dir=line["dir"]))
    # text-only redactions: a 0.6-pt square at each glyph centre, so neighbouring glyphs survive
    for e in edits:
        for cx, cy in e["cen"]:
            page.add_redact_annot(fitz.Rect(cx - 0.3, cy - 0.3, cx + 0.3, cy + 0.3))
    page.apply_redactions(images=fitz.PDF_REDACT_IMAGE_NONE, graphics=fitz.PDF_REDACT_LINE_ART_NONE)

    fr, fb = fitz.Font(fontfile=ARIAL), fitz.Font(fontfile=ARIAL_B)
    log = []
    for e in edits:
        log.append((e["old"], e["new"], round(e["bbox"].x0, 1), round(e["bbox"].y0, 1)))
        if not e["new"]:
            continue
        col = tuple(((e["color"] >> k) & 255) / 255 for k in (16, 8, 0))
        font = fb if e["bold"] else fr
        size = e["size"]
        x, y = e["origin"]
        tw = fitz.TextWriter(page.rect, color=col)
        if e["dir"][0] != 1:                         # rotated y-axis title, reads bottom to top
            tw.append((x, y), e["new"], font=font, fontsize=size)
            tw.write_text(page, morph=(fitz.Point(x, y), fitz.Matrix(90)))
            continue
        if "^" in e["new"]:                          # P value with superscript exponent
            base, sup = e["new"].split("^")
            tw.append((x, y), base, font=font, fontsize=size)
            w = font.text_length(base, fontsize=size)
            tw.append((x + w + 0.2, y - size * 0.38), sup, font=font, fontsize=size * 0.7)
            tw.write_text(page)
            continue
        if (34.5 < e["bbox"].x1 < 35.6) or (280.5 < e["bbox"].x1 < 281.5):   # right-aligned row labels
            x = e["bbox"].x1 - font.text_length(e["new"], fontsize=size)
        tw.append((x, y), e["new"], font=font, fontsize=size)
        tw.write_text(page)
    doc.save(HERE / "edited.pdf", garbage=3, deflate=True)
    with open(HERE / "edit_log.tsv", "w") as fh:
        fh.write("old\tnew\tx0\ty0\n")
        for r in log:
            fh.write("\t".join(map(str, r)) + "\n")
    print(len(edits), "edits")

    # ---- recomposition -----------------------------------------------------------
    # 5g: replace the dSV ring (2,657 -> 1,480 dSVs); see replace_5g_dsv.py
    import replace_5g_dsv
    replace_5g_dsv.main(HERE / "edited.pdf", HERE / "edited_5g.pdf")
    src = fitz.open(HERE / "edited_5g.pdf")
    W, H = src[0].rect.width, src[0].rect.height
    cut = 187.6                                     # blank band between panel a and row b–d
    scale = 0.912
    lower_h = (H - cut) * scale
    out = fitz.open()
    pg = out.new_page(width=W, height=cut + lower_h)
    pg.show_pdf_page(fitz.Rect(0, 0, W, cut), src, 0, clip=fitz.Rect(0, 0, W, cut))
    lw = W * scale
    x0 = (W - lw) / 2
    pg.show_pdf_page(fitz.Rect(x0, cut, x0 + lw, cut + lower_h), src, 0, clip=fitz.Rect(0, cut, W, H))
    out.save(OUT / "Figure5.pdf", garbage=3, deflate=True)
    print("page mm", round(W / 72 * 25.4, 1), "x", round((cut + lower_h) / 72 * 25.4, 1))
    # final pass: true minus signs (5f axis) and FL/TN-Hap names; prevents regressions on re-run
    import minus_fix
    minus_fix.fix(OUT / "Figure5.pdf")
    # 5i colour-bar title: ln(1 + dSV + dSNP), per author Source Data 4 (African35 workbook)
    import colorbar_title_5i
    colorbar_title_5i.add_title(OUT / "Figure5.pdf")
    # 5i: in-figure notes moved to the legend (author request 2026-09-23)
    import remove_5i_notes
    remove_5i_notes.remove(OUT / "Figure5.pdf")


if __name__ == "__main__":
    main()
