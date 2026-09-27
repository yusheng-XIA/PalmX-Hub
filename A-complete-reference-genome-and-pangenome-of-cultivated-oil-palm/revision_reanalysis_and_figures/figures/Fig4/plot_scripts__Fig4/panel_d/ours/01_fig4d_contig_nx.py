#!/usr/bin/env python3
"""Fig. 4d: replace the 8 RagTag haplotype curves (scaffold Nx) with contig Nx
(RagTag scaffolds split at every N run), and recompute the 39-assembly median
(dashed dark-grey line), which the source script draws as np.percentile(all 39, 50).
Only path coordinates change; colour, dash, width, cap/join operators are untouched.
Axis mapping (fitted on all 12 highlighted curves, max error 0.003 pt):
  x: N10..N100 -> 32.168 + (i * 14.8463) ; y(form) = 407.9044 + 0.283465 * Mb
in: steps/step0.pdf  out: steps/step1_4d.pdf  +  sd/Fig4d_Nx.tsv
"""
import csv, re, sys
import numpy as np, openpyxl, fitz
from common import *

SRC, DST = FIX / "steps/step0.pdf", FIX / "steps/step1_4d.pdf"
import sys as _s
if len(_s.argv) > 2: SRC, DST = Path(_s.argv[1]), Path(_s.argv[2])
K, Y0, X0 = 0.283465, 407.9044, 32.168
XS = [0, 14.844, 29.691, 44.539, 59.387, 74.23, 89.078, 103.926, 118.773, 133.617]   # original x offsets
NX = FIX / "src/nx_ragtag_full.tsv"
RAGTAG = ["dura_hap1", "dura_hap2", "pisifera_hap1", "pisifera_hap2",
          "nrly_hap1", "nrly_hap2", "meizhou4_hap1", "meizhou4_hap2"]
DISPLAY = {"dura": "TK", "pisifera": "NS", "nrly": "Nigerian", "meizhou4": "E. oleifera",
           "bk": "TN", "American": "FL", "Africa": "FL"}
COLOR = {"dura": ".275 .404 .682", "pisifera": ".949 .718 .416",
         "nrly": ".161 .686 .831", "meizhou4": ".655 .851 .867"}

# ---- data -------------------------------------------------------------------------
new = {}
for r in csv.DictReader(open(NX), delimiter="\t"):
    if r["level"] == "contig_splitN":
        new[r["Sample"]] = r
wb = openpyxl.load_workbook(SD_XLSX, read_only=True)
sd_rows = list(wb["Fig.4d_Nx"].iter_rows(values_only=True))
hdr, sd_rows = sd_rows[0], sd_rows[1:]
assert len(sd_rows) == 39
table = []
for r in sd_rows:
    s = str(r[0])
    if s in new:
        n = new[s]
        vals = [int(n[f"N{k}"]) for k in range(10, 101, 10)]
        table.append(dict(Sample=s, vals=vals, Total=int(n["Total"]), N_contigs=int(n["N_seqs"]),
                          Largest=int(n["Largest"]), Level="contig (RagTag scaffold split at N runs)",
                          Source="03_V3/01_figure1/00_final_input_resources/04_ragtag_6haps/fasta/%s.ragtag.fasta" % s,
                          Changed="yes"))
    else:
        table.append(dict(Sample=s, vals=list(r[1:11]), Total=r[11], N_contigs=r[12], Largest=r[13],
                          Level="contig (unchanged, contig_Nx_stats_full.tsv)", Source="Source Data Fig.4d_Nx (unchanged)",
                          Changed="no"))
assert sum(t["Changed"] == "yes" for t in table) == 8
Y = np.array([t["vals"] for t in table], float) / 1e6
med_old = np.percentile(np.array([r[1:11] for r in sd_rows], float) / 1e6, 50, axis=0)
med_new = np.percentile(Y, 50, axis=0)


def path(vals_mb):
    y = [Y0 + K * v for v in vals_mb]
    body = " ".join(f"{num(XS[i])} {num(y[i] - y[0])} {'m' if i == 0 else 'l'}" for i in range(10))
    return num(y[0], 4), body


# ---- edit ---------------------------------------------------------------------------
doc = fitz.open(SRC)
s = get_stream(doc)
a = s.index("q 1 0 0 1 32.168 459.0024 cm")                       # first highlighted curve
b = s.index("Q q .333 .333 .333 RG 2 J 0 j .5 w 29.195 487.275", a)
blk = s[a:b]
dash_at = blk.index("[4 1.6]0 d")
pat = re.compile(r"q 1 0 0 1 32\.168 ([\d.]+) cm ([\d. ]+) RG ((?:0 J )?1 j 1 w )(0 0 m .*? l) S Q")
out, pos, done = [], 0, []
for m in pat.finditer(blk):
    col = m.group(2)
    mat = next((k for k, c in COLOR.items() if c == col), None)
    out.append(blk[pos:m.start()])
    if mat is None:                     # FL / TN curves stay
        out.append(m.group(0))
    else:
        hap = "hap1" if m.start() < dash_at else "hap2"
        name = f"{mat}_{hap}"
        vals = [float(new[name][f"N{k}"]) / 1e6 for k in range(10, 101, 10)]
        y0, body = path(vals)
        out.append(f"q 1 0 0 1 32.168 {y0} cm {col} RG {m.group(3)}{body} S Q")
        done.append(name)
    pos = m.end()
out.append(blk[pos:])
assert sorted(done) == sorted(RAGTAG), done
s = s[:a] + "".join(out) + s[b:]
put_stream(doc, s)

# median line lives in its own form (xref 2002)
ms = get_stream(doc, 2002)
m = re.search(r"q 1 0 0 1 32\.168 ([\d.]+) cm (0 0 m .*? l) S Q", ms)
old_pts = [float(m.group(1)) + float(v) for v in re.findall(r"-?[\d.]+ (-?[\d.]+) [ml]", m.group(2))]
assert np.abs(np.array(old_pts) - (Y0 + K * med_old)).max() < 0.01, "median line is not the 39-assembly median"
y0, body = path(med_new)
ms = ms[:m.start()] + f"q 1 0 0 1 32.168 {y0} cm {body} S Q" + ms[m.end():]
bbox_top = float(doc.xref_get_key(2002, "BBox")[1].strip("[]").split()[1])
assert Y0 + K * med_new.max() < bbox_top
put_stream(doc, ms, 2002)
doc.save(DST, garbage=0, deflate=True)

# ---- Source Data table ---------------------------------------------------------------
with open(FIX / "sd/Fig4d_Nx.tsv", "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["Sample", "Display_name", "Haplotype", *[f"N{k}" for k in range(10, 101, 10)],
                "Total", "N_contigs", "Largest", "Level", "Changed_vs_submitted", "Source"])
    for t in table:
        s0 = t["Sample"]
        key = s0.split("_")[0]
        disp = DISPLAY.get(key, "EG_houke" if s0.startswith("houke") else f"EG_{int(s0):03d}" if s0.isdigit() else s0)
        hap = s0.split("_")[1].replace("hap", "Hap") if "_hap" in s0 else ""
        w.writerow([s0, disp, hap, *t["vals"], t["Total"], t["N_contigs"], t["Largest"], t["Level"], t["Changed"], t["Source"]])
    w.writerow(["Median_39", "Median of 39 assemblies (dashed grey line)", "", *[int(round(v * 1e6)) for v in med_new],
                "", "", "", "np.percentile(..., 50) over the 39 rows above", "yes (recomputed)", "derived"])
print("old median (Mb)", med_old.round(2))
print("new median (Mb)", med_new.round(2))
print("max new contig value (Mb)", Y.max().round(2), "y-axis top ~", round((487.275 - Y0) / K, 1))
