#!/usr/bin/env python3
"""Fig. 4d: replace only the two Nigerian (nrly_hap1/nrly_hap2) contig-Nx curves with the final 16-chromosome
assembly + unplaced contigs (split at N runs; fix/nrly_qc/nrly_final_contig_Nx_all_variants.tsv, variant
final_plus_unplaced) and recompute the 39-assembly median (dashed line, np.percentile(...,50)).
Adapted from deliver/Main_Figures_revised/Figure4_trace_fix_scripts/01_fig4d_contig_nx.py (same axis mapping,
same x offsets, same path/operator syntax); input = the delivered Figure4.pdf (already carries the contig-Nx
curves of the 8 RagTag haplotypes). Writes Figure4_candidate.pdf and Fig4d_Nx.tsv (complete sheet)."""
import csv, re
from pathlib import Path
import numpy as np, openpyxl, fitz

HERE = Path(__file__).resolve().parent
S = HERE.parents[1]
SRC = S / "deliver/Main_Figures_revised/Figure4.pdf"
DST = HERE / "Figure4_candidate.pdf"
SD_XLSX = S / "deliver/Source_Data_split/Source_Data_Fig4.xlsx"
NX = S / "fix/nrly_qc/nrly_final_contig_Nx_all_variants.tsv"
MAIN_FORM, MED_FORM = 78, 2002
K, Y0 = 0.283465, 407.9044
XS = [0, 14.844, 29.691, 44.539, 59.387, 74.23, 89.078, 103.926, 118.773, 133.617]
NRLY_COL = ".161 .686 .831"
SRCPATH = "${CLUSTER_WORK}/nrly_qc/ext/%s.final_plus_unplaced.fa"


def num(v, nd=3):
    s = f"{v:.{nd}f}".rstrip("0").rstrip(".")
    if s.startswith("0."): s = s[1:]
    elif s.startswith("-0."): s = "-" + s[2:]
    return "0" if s in ("", "-", "-0") else s


def path(vals_mb):
    y = [Y0 + K * v for v in vals_mb]
    return num(y[0], 4), " ".join(f"{num(XS[i])} {num(y[i] - y[0])} {'m' if i == 0 else 'l'}" for i in range(10))


new = {r["Sample"]: r for r in csv.DictReader(open(NX), delimiter="\t") if r["Variant"] == "final_plus_unplaced"}
assert sorted(new) == ["nrly_hap1", "nrly_hap2"]
ws = openpyxl.load_workbook(SD_XLSX)["Fig.4d_Nx"]
rows = [list(r) for r in ws.iter_rows(values_only=True)]
hdr, body = rows[0], rows[1:]
med_row = body[-1]; body = body[:-1]
assert med_row[0] == "Median_39" and len(body) == 39
old_Y = np.array([r[3:13] for r in body], float) / 1e6
med_old = np.percentile(old_Y, 50, axis=0)
assert np.abs(np.array(med_row[3:13], float) / 1e6 - med_old).max() < 1e-6, "sheet median != percentile of sheet"
changed = []
for r in body:
    if r[0] in new:
        n = new[r[0]]
        old = list(r)
        r[3:13] = [int(n[f"N{k}"]) for k in range(10, 101, 10)]
        r[13], r[14], r[15] = int(n["Total"]), int(n["N_contigs"]), int(n["Largest"])
        r[16] = "contig (final 16-chromosome assembly + unplaced contigs, split at N runs)"
        r[17] = "yes"
        r[18] = SRCPATH % r[0]
        changed.append((old, list(r)))
Y = np.array([r[3:13] for r in body], float) / 1e6
med_new = np.percentile(Y, 50, axis=0)

doc = fitz.open(SRC)
s = doc.xref_stream(MAIN_FORM).decode("latin1")
pat = re.compile(r"q 1 0 0 1 32\.168 ([\d.]+) cm " + re.escape(NRLY_COL) + r" RG ((?:0 J )?1 j 1 w )(0 0 m .*? l) S Q")
hits = list(pat.finditer(s))
assert len(hits) == 2, len(hits)
dash_at = s.index("[4 1.6]0 d", s.index("q 1 0 0 1 32.168 459.0024 cm"))
assert hits[0].start() < dash_at < hits[1].start()          # solid = Hap1, dashed = Hap2 (as in 01_fig4d)
cur = {r[0]: np.array(r[3:13], float) / 1e6 for (r, _) in [(o, 0) for o, _ in changed]}
out, pos = [], 0
for m, name in zip(hits, ["nrly_hap1", "nrly_hap2"]):
    pts = [float(m.group(1)) + float(v) for v in re.findall(r"-?[\d.]+ (-?[\d.]+) [ml]", m.group(3))]
    assert np.abs(np.array(pts) - (Y0 + K * cur[name])).max() < 0.01, (name, "curve != current Source Data")
    y0, b = path([float(new[name][f"N{k}"]) / 1e6 for k in range(10, 101, 10)])
    out += [s[pos:m.start()], f"q 1 0 0 1 32.168 {y0} cm {NRLY_COL} RG {m.group(2)}{b} S Q"]; pos = m.end()
out.append(s[pos:])
doc.update_stream(MAIN_FORM, "".join(out).encode("latin1"))
ms = doc.xref_stream(MED_FORM).decode("latin1")
m = re.search(r"q 1 0 0 1 32\.168 ([\d.]+) cm (0 0 m .*? l) S Q", ms)
old_pts = [float(m.group(1)) + float(v) for v in re.findall(r"-?[\d.]+ (-?[\d.]+) [ml]", m.group(2))]
assert np.abs(np.array(old_pts) - (Y0 + K * med_old)).max() < 0.01, "median line != 39-assembly median"
y0, b = path(med_new)
doc.update_stream(MED_FORM, (ms[:m.start()] + f"q 1 0 0 1 32.168 {y0} cm {b} S Q" + ms[m.end():]).encode("latin1"))
doc.save(DST, garbage=3, deflate=True)

with open(HERE / "Fig4d_Nx.tsv", "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t", lineterminator="\n")
    w.writerow(hdr)
    for r in body: w.writerow(["" if v is None else v for v in r])
    mr = list(med_row); mr[3:13] = [int(round(v * 1e6)) for v in med_new]
    w.writerow(["" if v is None else v for v in mr])
for o, n in changed:
    print(o[0], "N10..N100/Total/N_contigs old:", o[3:16]); print(o[0], "                       new:", n[3:16])
print("median old:", [int(v) for v in med_row[3:13]])
print("median new:", [int(round(v * 1e6)) for v in med_new])
print("wrote", DST)
