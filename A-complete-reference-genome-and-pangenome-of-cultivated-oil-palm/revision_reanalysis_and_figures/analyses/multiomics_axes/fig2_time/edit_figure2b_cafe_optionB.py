#!/usr/bin/env python3
"""Option B for Fig. 2b: replace every CAFE5 number with the re-run on the final MCMCTree tree
(runs/gamma_newMCMC; values in sd/optionB/Fig2b_gene_families.tsv). Input = step 1 output:
  edit_figure2_time.py deliver/Main_Figures_revised/Figure2.pdf work/Figure2_optionB_step1.pdf --2b-geometry \
      --node-circles sd/optionB/Fig2b_gene_families.tsv
(2a labels, 2b geometry on the final run, node circles by the plot_palm_v6.R rule from the re-run nets). Every replacement keeps the digit count, so glyph advance and position are
unchanged (Arial digits share one advance width); left edge, baseline, size, colour and weight are kept.

usage: edit_figure2b_cafe_optionB.py work/Figure2_optionB_step1.pdf optionB/Figure2_candidate_optionB.pdf [--log LOG.tsv]
"""
import csv, re, sys
from pathlib import Path
import fitz

H = Path(__file__).resolve().parent
SRC, OUT = sys.argv[1], sys.argv[2]
LOG = sys.argv[sys.argv.index("--log") + 1] if "--log" in sys.argv else None
AR = "/System/Library/Fonts/Supplemental/Arial.ttf"
M = "−"
old = {int(r["node_id"]): (int(r["expanded_families"]), int(r["contracted_families"]))
       for r in csv.DictReader(open(H / "../plot_scripts/Fig2/whole_figure_ours/Fig2A_CAFE5_branch_turnover_source.tsv"), delimiter="\t")}
new = {int(r["node_id"]): (int(r["expanded_families"]), int(r["contracted_families"]))
       for r in csv.DictReader(open(H / "sd/optionB/Fig2b_gene_families.tsv"), delimiter="\t")}
# internal nodes: (author node id, y of the '+' span, y of the '−' span) in the current figure
INTERNAL = {2: (207.1, 207.8), 4: (215.9, 224.1), 6: (228.3, 236.4), 8: (239.6, 248.9), 10: (251.7, 261.3),
            14: (265.1, 273.5), 15: (298.3, 306.4), 18: (281.4, 289.2), 19: (323.4, 332.7), 20: (302.4, 311.6)}
TIPS = [0, 1, 3, 5, 7, 9, 11, 12, 13, 16, 17, 21]

doc = fitz.open(SRC); pg = doc[0]
spans = [s for b in pg.get_text("dict")["blocks"] for l in b.get("lines", []) for s in l["spans"]
         if s["bbox"][0] < 262 and 199 < s["bbox"][1] < 361]
edits, log = [], []

def one(text, y):
    hit = [s for s in spans if s["text"].strip() == text and abs(s["bbox"][1] - y) < 1.2]
    if len(hit) != 1:
        raise SystemExit(f"{text!r} at y~{y}: {len(hit)} hits")
    return hit[0]

def rewrite(s, t, node):
    if len(t) != len(s["text"].strip()):
        raise SystemExit(f"length change {s['text']!r} -> {t!r}")
    ox, oy = s["origin"]
    lead = len(s["text"]) - len(s["text"].lstrip())   # a merged leading space (e.g. ' −349') keeps its advance
    ox += fitz.Font(fontfile=AR).text_length(" " * lead, fontsize=s["size"]); yc = oy - 0.33 * s["size"]
    pg.add_redact_annot(fitz.Rect(s["bbox"][0] + 0.3, yc - 0.12, s["bbox"][2] - 0.3, yc + 0.12), fill=False)
    c = s["color"]; rgb = ((c >> 16) & 255) / 255, ((c >> 8) & 255) / 255, (c & 255) / 255
    edits.append((ox, oy, t, s["size"], rgb)); log.append((node, s["text"], t, f"{ox:.2f},{oy:.2f}"))

for n, (yp, ym) in INTERNAL.items():
    (e0, c0), (e1, c1) = old[n], new[n]
    rewrite(one(f"+{e0}", yp), f"+{e1}", n)
    rewrite(one(f"{M}{c0}", ym), f"{M}{c1}", n)
for n in TIPS:
    (e0, c0), (e1, c1) = old[n], new[n]
    rewrite(one(f"+{e0}\xa0/\xa0{M}{c0}", next(s["bbox"][1] for s in spans if s["text"].strip() == f"+{e0}\xa0/\xa0{M}{c0}")),
            f"+{e1}\xa0/\xa0{M}{c1}", n)
pg.apply_redactions(images=fitz.PDF_REDACT_IMAGE_NONE, graphics=fitz.PDF_REDACT_LINE_ART_NONE, text=fitz.PDF_REDACT_TEXT_REMOVE)
for x, y, t, sz, rgb in edits:
    pg.insert_text((x, y), t, fontname="ArialRC", fontfile=AR, fontsize=sz, color=rgb)
doc.subset_fonts()
doc.save(OUT, garbage=4, deflate=True)
if LOG:
    with open(LOG, "w") as f:
        f.write("cafe_node\told\tnew\torigin\n")
        for r in log: f.write("\t".join(map(str, r)) + "\n")
print(len(log), "labels rewritten")
