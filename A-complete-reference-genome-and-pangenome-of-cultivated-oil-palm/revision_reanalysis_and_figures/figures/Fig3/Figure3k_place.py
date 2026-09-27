#!/usr/bin/env python3
"""Redact the visible Fig. 3k area of Figure3.pdf and place a redrawn panel at 1:1.
usage: place_fig3k.py panel.pdf out.pdf
The 'k' panel letter (x < 281.2) and the neighbouring 'Stage' label above (y < 380.7)
are outside the redaction rectangle and remain untouched."""
import sys
import fitz
sys.path.insert(0, __import__('os').path.dirname(__file__))
from render_fig3k import X0, X1, Y0, Y1

SRC = '${WORK_DIR}/deliver/Main_Figures_revised/Figure3.pdf'
panel, out = sys.argv[1], sys.argv[2]
doc = fitz.open(SRC); page = doc[0]
R = fitz.Rect(X0, Y0, X1, Y1)
page.add_redact_annot(R, fill=(1, 1, 1))
page.apply_redactions(images=fitz.PDF_REDACT_IMAGE_NONE,
                      graphics=fitz.PDF_REDACT_LINE_ART_REMOVE_IF_COVERED,
                      text=fitz.PDF_REDACT_TEXT_REMOVE)
src = fitz.open(panel)
page.show_pdf_page(R, src, 0, keep_proportion=False)
doc.save(out, garbage=4, deflate=True)
