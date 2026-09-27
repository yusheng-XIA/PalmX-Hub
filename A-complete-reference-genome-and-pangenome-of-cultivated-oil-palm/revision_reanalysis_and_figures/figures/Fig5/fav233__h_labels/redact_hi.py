"""Remove panels h,i content (keeping the panel letters) from a copy of Figure5.pdf. Used by edit_figure5hi_mask.py."""
import fitz
RECTS = [fitz.Rect(33.0, 454.5, 250.0, 466.5), fitz.Rect(20.0, 466.5, 250.0, 694.9),
         fitz.Rect(256.0, 454.5, 518.7, 694.9), fitz.Rect(250.0, 466.5, 256.0, 694.9)]
def redact(page):
    for r in RECTS: page.add_redact_annot(r, fill=False)
    page.apply_redactions(images=fitz.PDF_REDACT_IMAGE_NONE, graphics=fitz.PDF_REDACT_LINE_ART_REMOVE_IF_COVERED,
                          text=fitz.PDF_REDACT_TEXT_REMOVE)
