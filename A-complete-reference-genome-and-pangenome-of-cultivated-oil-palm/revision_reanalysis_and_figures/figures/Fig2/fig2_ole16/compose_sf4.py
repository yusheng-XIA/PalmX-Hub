#!/usr/bin/env python3
"""Supplementary Fig. 4 candidate: current SF4 (a-c, unchanged, vector) + new panel d (former Fig. 2c heatmap)."""
from pathlib import Path
import fitz
HERE = Path(__file__).resolve().parent
S = HERE.parents[1]
src = fitz.open(S / "deliver/Supplementary_Information/Supplementary_Figures_PDF/Supplementary_Fig_04.pdf")
pan = fitz.open(HERE / "panels/sf4_metab.pdf")
w, h = src[0].rect.width, src[0].rect.height
pr = pan[0].rect
out = fitz.open()
pg = out.new_page(width=w, height=h + 6 + pr.height)
pg.show_pdf_page(fitz.Rect(0, 0, w, h), src, 0)
pg.show_pdf_page(fitz.Rect(0, h + 6, pr.width, h + 6 + pr.height), pan, 0)
out.save(HERE / "SF4_candidate.pdf", garbage=4, deflate=True)
pg.get_pixmap(dpi=300).save(HERE / "SF4_candidate.png")
print(f"SF4 candidate {w/72*25.4:.0f} x {(h+6+pr.height)/72*25.4:.0f} mm")
