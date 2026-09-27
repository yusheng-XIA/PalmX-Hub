"""Compose Figure 4 from the pre-repair Figure 4 (deliver/Main_Figures_revised/_superseded/Figure4_before_snprepair.pdf):
only the data areas of panels a-c are replaced by out/Fig4{a,b,c}_snp_repair.pdf; panel letters and all other panels untouched."""
import fitz
from pathlib import Path
H = Path(__file__).resolve().parent; S = H.parents[2]
SRC = S / "deliver/Main_Figures_revised/_superseded/Figure4_before_snprepair.pdf"
doc = fitz.open(SRC); pg = doc[0]
RED = [fitz.Rect(10, 12, 170, 114.8), fitz.Rect(172, 12, 342.5, 112.5), fitz.Rect(343.5, 12, 518.7, 112.5)]
PUT = {"a": fitz.Rect(0, 12, 170, 112.5), "b": fitz.Rect(172, 12, 342.5, 112.5), "c": fitz.Rect(343.53, 12, 518.53, 112.5)}
for r in RED: pg.add_redact_annot(r, fill=None)
pg.apply_redactions(images=fitz.PDF_REDACT_IMAGE_NONE, graphics=fitz.PDF_REDACT_LINE_ART_REMOVE_IF_COVERED, text=fitz.PDF_REDACT_TEXT_REMOVE)
for p, r in PUT.items(): pg.show_pdf_page(r, fitz.open(H / f"out/Fig4{p}_snp_repair.pdf"), 0)
doc.save(H / "out/Figure4_snp_repair.pdf", garbage=3, deflate=True)
print("saved", H / "out/Figure4_snp_repair.pdf")
