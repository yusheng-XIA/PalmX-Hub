"""Round-2 in-figure text edits (vector; ../r1/figtext.py): gene symbols in italics in Supplementary
Fig. 5e (legend) and 5i (y axis).  Outputs to fix/r2/figs; install_figs_r2.sh copies them into deliver/."""
import sys, collections
from pathlib import Path
import fitz
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "r1"))
import figtext as T

S = Path(__file__).resolve().parents[2]
O = Path(__file__).parent / "figs"
SF = S / "deliver/Supplementary_Information/Supplementary_Figures_PDF"
log = []


def run(src, name, edits, crops):
    d = fitz.open(src); p = d[0]
    changes = [T.edit_line(p, sel, segs, align=al, tag=f"r2{name[-2:]}{k}") for k, (sel, segs, al) in enumerate(edits)]
    out = O / f"{name}.pdf"
    d.save(out, garbage=3, deflate=True)
    T.export(src, out, png_out=O / f"{name}.png", docx_out=O / "docx_2400" / f"{name}.png")
    for tag, r in crops.items():
        T.before_after(src, out, r, O / f"{name}_{tag}_before_after.png")
    old_w = collections.Counter(T.text_words(src)); new_w = collections.Counter(T.text_words(out))
    assert old_w == new_w, (name, old_w - new_w, new_w - old_w)     # italics only: same words
    for o, n in changes:
        log.append((name, o.strip(), n.strip()))
    print(name, "ok", len(changes))


eq = lambda s: (lambda t, l: t.strip() == s)
sf5 = [(eq(f"OLE16{g} {m}"), [(f"OLE16{g}", "i"), (f" {m}", "r")], "left") for g in "ab" for m in ("FL", "TN")]
sf5.append((eq("Nuclei with OLE16a"), [("Nuclei with ", "r"), ("OLE16a", "i")], "center"))
run(SF / "Supplementary_Fig_05.pdf", "Supplementary_Fig_05", sf5, {"e": (240, 270, 350, 296), "i": (355, 550, 380, 620)})
with open(O / "fig_text_changes.tsv", "w") as f:
    f.write("figure_old_name\told\tnew\n")
    for r in log:
        f.write("\t".join(r) + "\n")
