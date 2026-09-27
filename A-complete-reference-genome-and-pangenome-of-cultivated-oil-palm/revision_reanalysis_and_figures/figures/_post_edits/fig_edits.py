"""extfix in-figure edits (vector text only; graphics untouched). Writes candidates to fix/extfix/figs/ and, with --install,
moves the current files to _superseded/*_before_extfix.* and installs PDF, 600-dpi PNG and 2400-px docx image."""
import sys, shutil
from pathlib import Path
from collections import Counter
S = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(S / "fix/r1"))
import figtext as FT, fitz
OUT = S / "fix/extfix/figs"; OUT.mkdir(exist_ok=True)
D = S / "deliver"
ED = D / "Extended_Data_Figures"; SFP = D / "Supplementary_Information/Supplementary_Figures_PDF"; SFG = D / "Supplementary_Information/figures_png"
DX = S / "docx_img"

def at(txt, x, y):
    return lambda t, l: t == txt and abs(l["bbox"][0] - x) < 1.5 and abs(l["bbox"][1] - y) < 1.5

JOBS = {
  # new ED1 = old ED1: panel order a,b,c,d(bottom),e(right of c) -> swap d/e so the reading order is a-e
  ("ED", "Extended_Data_Fig_01"): [
      (at("e", 322.8, 246.2), [("d", "b")], "left", {"d": 1, "e": 1}),
      (at("d", 0.0, 402.1), [("e", "b")], "left", None)],
  # new ED6 = old ED4: e (CV, top right) -> c; c (tree) -> d; d (ancestry bars) -> e
  ("ED", "Extended_Data_Fig_04"): [
      (at("e", 392.5, 4.5), [("c", "b")], "left", None),
      (at("c", 15.3, 157.7), [("d", "b")], "left", None),
      (at("d", 346.6, 157.7), [("e", "b")], "left", None)],
  # new ED4 = old SF13: 'Hybrid' -> 'TN'
  ("SF", "Supplementary_Fig_13"): [
      (lambda t, l: t == "Hybrid", [("TN", "r")], "center", None),
      (lambda t, l: t == "Hybrid 50%", [("TN 50%", "r")], "left", None)],
  # new ED8 = old SF11 d: separate the Pearson n from the Wilcoxon unit
  ("SF", "Supplementary_Fig_11"): [
      (lambda t, l: t == "n = 560 haplotype–chromosome", [("Pearson: ", "r"), ("n", "i"), (" = 560 haplotype–chromosome", "r")], "left", None),
      (lambda t, l: t == "combinations", [("combinations; Wilcoxon: dSV–carrier pairs", "r")], "left", None)],
}

def paths(kind, stem):
    if kind == "ED":
        return ED / f"{stem}.pdf", ED / f"{stem}.png", DX / "Extended_Data_Figures" / f"{stem}.png"
    return SFP / f"{stem}.pdf", SFG / f"{stem}.png", DX / "Supplementary_Figures" / f"{stem}.png"

install = "--install" in sys.argv
for (kind, stem), edits in JOBS.items():
    pdf, png, dx = paths(kind, stem)
    out_pdf = OUT / f"{stem}.pdf"
    doc = fitz.open(pdf); p = doc[0]
    log = []
    for k, (sel, segs, align, _) in enumerate(edits):
        log.append(FT.edit_line(p, sel, segs, align=align, tag=f"{stem[-2:]}{k}"))
    doc.save(out_pdf, garbage=3, deflate=True)
    FT.export(pdf, out_pdf, OUT / f"{stem}.png", OUT / f"{stem}_docx2400.png")
    a = Counter(FT.text_words(pdf)); b = Counter(FT.text_words(out_pdf))
    print(stem, log, "| removed", dict(a - b), "| added", dict(b - a))
    if install:
        sup = pdf.parent.parent / "_superseded" if kind == "SF" else ED / "_superseded"
        sup.mkdir(exist_ok=True)
        for f, suf in ((pdf, ".pdf"), (png, ".png")):
            shutil.copy2(f, sup / f"{stem}_before_extfix{suf}")
        shutil.copy2(dx, sup / f"{stem}_docx_before_extfix.png")
        shutil.copy2(out_pdf, pdf); shutil.copy2(OUT / f"{stem}.png", png); shutil.copy2(OUT / f"{stem}_docx2400.png", dx)
        print("  installed", pdf.name)
