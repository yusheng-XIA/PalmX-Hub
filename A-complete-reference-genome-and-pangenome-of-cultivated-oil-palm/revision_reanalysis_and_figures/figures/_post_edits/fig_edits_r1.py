"""Round-1 in-figure text edits (vector; see figtext.py).  Outputs to fix/r1/figs/<name>.pdf/.png(.tif),
docx_2400/<name>.png and before/after crops; install_figs_r1.sh copies them into deliver/ and docx_img/.
Figure numbers are the OLD file names (Supplementary_Fig_11 = Extended Data Fig. 8, etc.)."""
import sys, collections
from pathlib import Path
import fitz
sys.path.insert(0, str(Path(__file__).parent))
import figtext as T

S = Path(__file__).resolve().parents[2]
O = Path(__file__).parent / "figs"; (O / "docx_2400").mkdir(parents=True, exist_ok=True)
MAIN = S / "deliver/Main_Figures_revised"
SF = S / "deliver/Supplementary_Information/Supplementary_Figures_PDF"
ED = S / "deliver/Extended_Data_Figures"
log = []


def run(src, name, edits, crops, kind):
    d = fitz.open(src); p = d[0]
    changes = []
    for k, e in enumerate(edits):
        sel, segs, align = e[:3]
        changes.append(T.edit_line(p, sel, segs, align=align, tag=f"{name[:3]}{k}", size=(e[3] if len(e) > 3 else None)))
    out = O / f"{name}.pdf"
    d.save(out, garbage=3, deflate=True)
    if kind == "main":
        T.export(src, out, png_out=O / f"{name}.tif")
    else:
        T.export(src, out, png_out=O / f"{name}.png", docx_out=O / "docx_2400" / f"{name}.png")
    for tag, r in crops.items():
        T.before_after(src, out, r, O / f"{name}_{tag}_before_after.png")
    # only the targeted strings may differ
    old_w = collections.Counter(T.text_words(src)); new_w = collections.Counter(T.text_words(out))
    exp_rm = collections.Counter(w for o, n in changes for w in T.norm(o).split())
    exp_add = collections.Counter(w for o, n in changes for w in T.norm(n).split())
    assert old_w - new_w == exp_rm - exp_add, (name, old_w - new_w, exp_rm - exp_add)
    assert new_w - old_w == exp_add - exp_rm, (name, new_w - old_w, exp_add - exp_rm)
    for o, n in changes:
        log.append((name, o.strip(), n.strip()))
    print(name, "ok", len(changes))


eq = lambda s: (lambda t, l: t.strip() == s)
eqx = lambda s, xmin=None, xmax=None: (lambda t, l: t.strip() == s and (xmin is None or l["bbox"][0] >= xmin)
                                       and (xmax is None or l["bbox"][0] <= xmax))

# ---- Figure 2: plural fission/fusion labels (2a); italic gene symbol on y axes (2g, 2i)
fig2 = []
for t in ("16 Fission", "20 Fission", "28 Fission", "32 Fission", "33 Fission"):
    n = len([1 for _ in [0]])
fis = {"16 Fission": 2, "20 Fission": 1, "28 Fission": 1, "32 Fission": 1, "33 Fission": 1}
d0 = fitz.open(MAIN / "Figure2.pdf")[0]
labs = [l for l in T.lines(d0) if T.norm("".join(s["text"] for s in l["spans"])).strip().endswith(("Fission", "Fusion"))]
assert len(labs) == 12, len(labs)
for l in labs:
    t = T.norm("".join(s["text"] for s in l["spans"])).strip()
    x0 = l["bbox"][0]
    new = t.replace("Fission", "fissions").replace("Fusion", "fusions")
    new = " ".join(new.split())
    fig2.append((lambda tt, ll, t=t, x0=x0: tt.strip() == t and abs(ll["bbox"][0] - x0) < 0.5, [(new, "b")], "right"))
fig2.append((eq("OLE16a RNA"), [("OLE16a", "i"), (" RNA", "r")], "center"))
fig2.append((eq("OLE16a (RPM)"), [("OLE16a", "i"), (" (RPM)", "r")], "center"))
run(MAIN / "Figure2.pdf", "Figure2", fig2,
    {"2a": (140, 50, 510, 76), "2g": (140, 595, 170, 660), "2i": (350, 600, 372, 660)}, "main")

# ---- Figure 4: spelling consistent with the text (4e 'Pangenome', 4g 'Biallelic')
d4 = fitz.open(MAIN / "Figure4.pdf")[0]
pg = [l for l in T.lines(d4) if T.norm("".join(s["text"] for s in l["spans"])).strip() == "Pan-genome"]
assert len(pg) == 2
fig4 = [(lambda tt, ll, x0=l["bbox"][0]: tt.strip() == "Pan-genome" and abs(ll["bbox"][0] - x0) < 0.5,
         [("Pangenome", "r")], "left") for l in pg]
fig4.append((eq("Bi-allelic genes"), [("Biallelic genes", "r")], "left"))
run(MAIN / "Figure4.pdf", "Figure4", fig4, {"4e": (200, 118, 340, 145), "4g": (25, 232, 120, 250)}, "main")

# ---- Extended Data Fig. 8 (old SF11): panel f uses focal dSVs, not haplotype-chromosome combinations
sf11 = [
    (lambda t, l: t.strip() == "All38; n = 608" and l["bbox"][0] > 300, [("All38; ", "r"), ("n", "i"), (" = 1,480 dSVs", "r")], "center"),
    (lambda t, l: t.strip() == "Definition (c), 1/35; n = 560" and l["bbox"][0] > 300,
     [("Definition (c), 1/35; ", "r"), ("n", "i"), (" = 550 dSVs", "r")], "center"),
    (lambda t, l: t.startswith("Rare derived variants carried by the dSV haplotypes"),
     [("Rare variants carried by the dSV haplotypes (same carrier-frequency filters; date palm polarization unless stated)", "r")], "center"),
]
run(SF / "Supplementary_Fig_11.pdf", "Supplementary_Fig_11", sf11, {"ef": (90, 330, 500, 370)}, "sf")

# ---- post-harvest spelling (SF3, SF6, SF7)
run(SF / "Supplementary_Fig_03.pdf", "Supplementary_Fig_03", [
    (lambda t, l: t.strip() == "Postharvest" and abs(l["dir"][0]) < 0.99, [("Post-harvest", "r")], "right"),
    (lambda t, l: t.strip() == "Postharvest" and abs(l["dir"][0]) > 0.99, [("Post-harvest", "r")], "center"),
    (eq("Sampling stage (shaded: postharvest)"), [("Sampling stage (shaded: post-harvest)", "r")], "center"),
], {"c": (80, 300, 130, 340), "d": (200, 160, 260, 180), "g": (190, 505, 320, 525)}, "sf")
run(SF / "Supplementary_Fig_06.pdf", "Supplementary_Fig_06", [
    (eq("Postharvest"), [("Post-harvest", "r")], "center")], {"c": (255, 165, 315, 182)}, "sf")
run(SF / "Supplementary_Fig_07.pdf", "Supplementary_Fig_07", [
    (eq("Postharvest"), [("Post-harvest", "r")], "center"),
    (eq("Postharvest (6 stages)"), [("Post-harvest (6 stages)", "r")], "center")],
    {"a": (345, 4, 405, 20), "b": (300, 318, 395, 334)}, "sf")

# ---- Supplementary Fig. 4c: legend entries clear of the 30% cut-off line (5.8 -> 5.2 pt, left-aligned as before)
run(SF / "Supplementary_Fig_04.pdf", "Supplementary_Fig_04", [
    (eq("Positive (n = 12,742; median 8.4%)"), [("Positive (", "r"), ("n", "i"), (" = 12,742; median 8.4%)", "r")], "left", 5.2),
    (eq("Negative (n = 12,234; median 7.7%)"), [("Negative (", "r"), ("n", "i"), (" = 12,234; median 7.7%)", "r")], "left", 5.2),
], {"c": (110, 420, 260, 470)}, "sf")

# ---- Extended Data Fig. 2: gene symbols in italics
run(ED / "Extended_Data_Fig_02.pdf", "Extended_Data_Fig_02", [
    (eq("PDAT1 (evm.TU.chr03B.814)"), [("PDAT1", "i"), (" (evm.TU.chr03B.814)", "r")], "left"),
    (eq("HACD (evm.TU.chr05B.150)"), [("HACD", "i"), (" (evm.TU.chr05B.150)", "r")], "left"),
], {"b": (265, 14, 365, 30), "e": (12, 327, 110, 343)}, "sf")

(O / "fig_text_changes.tsv").write_text("figure\told\tnew\n" + "\n".join("\t".join(r) for r in log) + "\n", encoding="utf-8")
