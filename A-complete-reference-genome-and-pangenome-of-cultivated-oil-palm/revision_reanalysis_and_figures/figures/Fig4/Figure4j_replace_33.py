#!/usr/bin/env python3
"""Replace Fig. 4j (34 haplotype-level tips) with the 33-material redraw (f4j/panel_j_33.pdf).
The panel label 'j' and everything outside x 0-224.5 / y 347-620.8 are untouched.
Previous file: _superseded/Figure4_before_4j33.pdf."""
import fitz
from PIL import Image
SRC = "_superseded/Figure4_before_4j33.pdf"
PANEL = "../../f4j/panel_j_33.pdf"
doc = fitz.open(SRC); pg = doc[0]
for r in (fitz.Rect(11.5, 347, 224.5, 620.8), fitz.Rect(0, 367.5, 11.5, 620.8)):
    pg.add_redact_annot(r, fill=None)
pg.apply_redactions(images=fitz.PDF_REDACT_IMAGE_NONE, graphics=fitz.PDF_REDACT_LINE_ART_REMOVE_IF_COVERED,
                    text=fitz.PDF_REDACT_TEXT_REMOVE)
left = [d["rect"] for d in pg.get_drawings() if d["rect"].x1 <= 224.5 and d["rect"].y0 >= 347]
left_txt = [w[4] for w in pg.get_text("words", clip=fitz.Rect(0, 347, 224.5, 620.8))]
print("leftover drawings", len(left), "leftover words", left_txt)
pg.show_pdf_page(fitz.Rect(0, 347, 224.5, 620.8), fitz.open(PANEL), 0)
doc.save("Figure4.pdf", garbage=3, deflate=True)
pg = fitz.open("Figure4.pdf")[0]
pix = pg.get_pixmap(dpi=600, alpha=False)
Image.frombytes("RGB", (pix.width, pix.height), pix.samples).save("Figure4.tif", compression="tiff_lzw", dpi=(600, 600))
pg.get_pixmap(dpi=150, alpha=False).save("Figure4_preview.png")
pg.get_pixmap(dpi=220, alpha=False, clip=fitz.Rect(0, 330, 518.7, 620.8)).save("../../f4j/Figure4_jk_after.png")
