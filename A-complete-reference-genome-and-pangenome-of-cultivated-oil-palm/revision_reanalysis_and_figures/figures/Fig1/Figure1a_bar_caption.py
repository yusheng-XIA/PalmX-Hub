#!/usr/bin/env python3
"""1a bar caption: 'Share of the selected major-crop / oil totals' -> 'Share of nine major oil crops'
(same font, size, colour and baseline; left edge kept)."""
import fitz
from PIL import Image
ARIAL = "/System/Library/Fonts/Supplemental/Arial.ttf"
doc = fitz.open("_superseded/Figure1_fao_shares_step1.pdf"); pg = doc[0]
sp = [s for b in pg.get_text("dict")["blocks"] for l in b.get("lines", []) for s in l["spans"]
      if "selected major-crop" in s["text"]]
assert len(sp) == 1, sp
s = sp[0]; print(s["font"], s["size"], hex(s["color"]), s["bbox"])
pg.add_redact_annot(fitz.Rect(s["bbox"]), fill=None)
pg.apply_redactions(images=fitz.PDF_REDACT_IMAGE_NONE, graphics=fitz.PDF_REDACT_LINE_ART_NONE)
pg.insert_font(fontname="ArialMT", fontfile=ARIAL)
c = s["color"]; col = ((c >> 16 & 255) / 255, (c >> 8 & 255) / 255, (c & 255) / 255)
pg.insert_text(s["origin"], "Share of nine major oil crops", fontname="ArialMT", fontsize=s["size"], color=col)
doc.save("Figure1.pdf", garbage=3, deflate=True)
pg = fitz.open("Figure1.pdf")[0]
pix = pg.get_pixmap(dpi=600, alpha=False)
Image.frombytes("RGB", (pix.width, pix.height), pix.samples).save("Figure1.tif", compression="tiff_lzw", dpi=(600, 600))
pg.get_pixmap(dpi=150, alpha=False).save("Figure1_preview.png")
pg.get_pixmap(dpi=500, alpha=False, clip=fitz.Rect(320, 20, 518, 90)).save("../../fix/fig1a_shares/bars_final.png")
