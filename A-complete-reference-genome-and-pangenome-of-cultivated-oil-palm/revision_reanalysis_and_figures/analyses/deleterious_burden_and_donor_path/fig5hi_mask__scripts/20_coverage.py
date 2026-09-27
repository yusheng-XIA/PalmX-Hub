#!/usr/bin/env python3
"""Per-500-kb-window reference coverage of one donor assembly, from the PAF that produced its dSNP calls.
Alignments are filtered exactly as paftools.js call does before calling (-l 1000: aligned length >= 1000;
MAPQ >= 5; primary only: tp:A:S / tp:A:i dropped, s1 without s2 dropped); covered = union of target intervals.
usage: 20_coverage.py LABEL SORTED_PAF OUT_TSV
"""
import sys, re
from collections import defaultdict
SNP_FAI = "${ANALYSIS_DIR}/14_pan_genome/04_SNP_calling/00_renamed_genomes/Africa_hap2.fasta.fai"
WIN = 500_000
label, paf, out = sys.argv[1:4]
L = {l.split()[0]: int(l.split()[1]) for l in open(SNP_FAI)}
iv = defaultdict(list); n = 0; kept = 0
for line in open(paf):
    t = line.split("\t"); n += 1
    if len(t) < 12 or t[5] == "*": continue
    if int(t[10]) < 1000 or int(t[11]) < 5: continue
    tags = line.rstrip("\n").split("\t")[12:]
    tp = next((x[5:] for x in tags if x.startswith("tp:A:")), None)
    s1 = any(x.startswith("s1:i:") for x in tags); s2 = any(x.startswith("s2:i:") for x in tags)
    if s1 and not s2: continue
    if tp in ("S", "i"): continue
    iv[t[5]].append((int(t[7]), int(t[8]))); kept += 1
with open(out, "w") as fh:
    fh.write("Donor\tChrom\tWindow_Index\tWindow_bp\tCovered_bp\tCoverage\n")
    for c in sorted(L, key=lambda c: int(c[3:5])):
        nw = (L[c] + WIN - 1) // WIN; cov = [0] * nw
        s = sorted(iv.get(c, [])); cur_s = cur_e = None; merged = []
        for a, b in s:
            if cur_e is None or a > cur_e: 
                if cur_e is not None: merged.append((cur_s, cur_e))
                cur_s, cur_e = a, b
            else: cur_e = max(cur_e, b)
        if cur_e is not None: merged.append((cur_s, cur_e))
        for a, b in merged:
            w = a // WIN
            while w * WIN < b and w < nw:
                x = min(b, (w + 1) * WIN) - max(a, w * WIN)
                if x > 0: cov[w] += x
                w += 1
        for w in range(nw):
            wb = min(WIN, L[c] - w * WIN)
            fh.write(f"{label}\t{c}B\t{w}\t{wb}\t{cov[w]}\t{cov[w]/wb:.6f}\n")
print(label, n, kept)
