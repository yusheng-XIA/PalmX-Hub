#!/usr/bin/env python3
"""3a/3b legend: 'Hap-specific' -> 'No 1:1 allele' (genes without a one-to-one allele in the synteny + GMAP
pairing; author decision 2026-09-24). Previous file: _superseded/Figure3_before_3b_label.pdf"""
import fitz
from PIL import Image
ARIAL = "/System/Library/Fonts/Supplemental/Arial.ttf"
doc = fitz.open("_superseded/Figure3_before_3b_label.pdf"); pg = doc[0]
hit = [w for w in pg.get_text("words") if w[4] == "Hap-specific"]
assert len(hit) == 1, hit
sp = [s for b in pg.get_text("dict")["blocks"] for l in b.get("lines", []) for s in l["spans"] if s["text"].strip() == "Hap-specific"][0]
pg.add_redact_annot(fitz.Rect(hit[0][:4]), fill=None)
pg.apply_redactions(images=fitz.PDF_REDACT_IMAGE_NONE, graphics=fitz.PDF_REDACT_LINE_ART_NONE)
pg.insert_font(fontname="ArialMT", fontfile=ARIAL)
c = sp["color"]; col = ((c >> 16 & 255) / 255, (c >> 8 & 255) / 255, (c & 255) / 255)
pg.insert_text(sp["origin"], "No 1:1 allele", fontname="ArialMT", fontsize=sp["size"], color=col)
doc.save("Figure3.pdf", garbage=3, deflate=True)
pg = fitz.open("Figure3.pdf")[0]
pix = pg.get_pixmap(dpi=600, alpha=False)
Image.frombytes("RGB", (pix.width, pix.height), pix.samples).save("Figure3.tif", compression="tiff_lzw", dpi=(600, 600))
pg.get_pixmap(dpi=150, alpha=False).save("Figure3_preview.png")
pg.get_pixmap(dpi=500, alpha=False, clip=fitz.Rect(0, 20, 200, 130)).save("../../fix/fig3k/fig3ab_label.png")
