#!/usr/bin/env python3
"""Fig. 1b: TN-Hap1 BUSCO 99.6 -> 99.7 to match ST1/ST3 (odb12 run; V19 marks the TN value 99.56 'table-only').
Reversible: the previous PDF is kept as _superseded/Figure1_before_TN_BUSCO.pdf."""
import fitz
from PIL import Image
ARIAL = "/System/Library/Fonts/Supplemental/Arial.ttf"
src = "_superseded/Figure1_before_TN_BUSCO.pdf"
doc = fitz.open(src); pg = doc[0]
hit = [w for w in pg.get_text("words") if w[4] == "99.6" and abs(w[1] - 308.68) < 0.5]
assert len(hit) == 1, hit
r = fitz.Rect(hit[0][:4])
pg.add_redact_annot(r, fill=None)
pg.apply_redactions(images=fitz.PDF_REDACT_IMAGE_NONE, graphics=fitz.PDF_REDACT_LINE_ART_NONE)
pg.insert_font(fontname="ArialMT", fontfile=ARIAL)
pg.insert_text((375.4245910644531, 314.2841796875), "99.7", fontname="ArialMT", fontsize=5.568599700927734,
               color=(0x22/255,)*3)
doc.save("Figure1.pdf", garbage=3, deflate=True)
doc = fitz.open("Figure1.pdf"); pg = doc[0]
print([w[4] for w in pg.get_text("words") if 375 < w[0] < 376 and 300 < w[1] < 330])
pix = pg.get_pixmap(dpi=600, alpha=False)
Image.frombytes("RGB", (pix.width, pix.height), pix.samples).save("Figure1.tif", compression="tiff_lzw", dpi=(600, 600))
pg.get_pixmap(dpi=150, alpha=False).save("Figure1_preview.png")
pg.get_pixmap(dpi=600, alpha=False, clip=fitz.Rect(340, 190, 400, 390)).save("../../tmp_fig1b_zoom.png")
