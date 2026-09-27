#!/usr/bin/env python3
"""stdin: bcftools query -H -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' (per-group AN/AC computed here).
Counts, for each test group X vs reference R, biallelic sites where one allele has frequency >= 0.05 in X and
< 0.01 in R (and == 0 in R), requiring called-allele fraction >= MINCALL in both groups (denominator: samples of the
group not listed as absent for this chromosome).
argv: tag chrom groups.txt absent.tsv(sample chrom; '-' none) mincall out_prefix"""
import sys, gzip, collections
tag, chrom, gfile, afile, mincall, outp = sys.argv[1:7]
mincall = float(mincall)
GROUPS = ["ALL", "COMM", "NONC", "AFR", "HHG", "IDB", "SAEG", "SEAA", "SEAB", "K4P1", "K4P2", "K4P3", "K4P4"]
absent = set()
if afile != "-":
    for l in open(afile):
        f = l.split()
        if len(f) >= 2 and f[1] == chrom: absent.add(f[0])
npres = collections.Counter()
for l in open(gfile):
    s, g = l.split()
    if s in absent: continue
    for x in g.split(","): npres[x] += 1
maxan = {g: 2 * npres[g] for g in GROUPS}
PAIRS = [("NONC", "COMM"), ("AFR", "COMM"), ("IDB", "COMM"), ("SEAB", "COMM"), ("HHG", "COMM"),
         ("K4P1", "K4P2"), ("K4P3", "K4P2"), ("K4P4", "K4P2"), ("COMM", "NONC")]
idx = {g: i for i, g in enumerate(GROUPS)}
C = collections.Counter()
lenbins = collections.Counter()
out = gzip.open(outp + ".sites.tsv.gz", "wt")
out.write("chrom\tpos\tref\talt\tallele\tfreq_NONC\tfreq_COMM\tfreq_AFR\tfreq_IDB\tfreq_SEAB\tfreq_HHG\n")
nsite = 0
import numpy as np
memb = {}
for l in open(gfile):
    s_, g_ = l.split()
    memb[s_] = set(g_.split(","))
hdr = sys.stdin.readline().rstrip("\n").split("\t")
names = [h.split("]", 1)[1].split(":")[0] for h in hdr[4:]]
M = np.array([[1 if g in memb[n] else 0 for g in GROUPS] for n in names], dtype=np.int64)
CODE = {}
for a in ("0", "1", "."):
    for b in ("0", "1", "."):
        for sep in ("/", "|"):
            k = a + sep + b
            CODE[k] = ((a != ".") + (b != "."), (a == "1") + (b == "1"))
CODE["0"] = (1, 0); CODE["1"] = (1, 1); CODE["."] = (0, 0)
for line in sys.stdin:
    f = line.rstrip("\n").split("\t")
    pos, ref, alt = f[1], f[2], f[3]
    cc = [CODE[x] for x in f[4:]]
    A = np.array(cc, dtype=np.int64)
    an = (A[:, 0] @ M).tolist(); ac = (A[:, 1] @ M).tolist()
    nsite += 1
    ok = [maxan[g] > 0 and an[i] >= mincall * maxan[g] for i, g in enumerate(GROUPS)]
    fr = [ac[i] / an[i] if an[i] else float("nan") for i in range(len(GROUPS))]
    if ok[0]:
        C[("callable_ALL",)] += 1
    for X, R in PAIRS:
        i, j = idx[X], idx[R]
        if not (ok[i] and ok[j]): continue
        C[(X, R, "callable")] += 1
        fx, fr_ = fr[i], fr[j]
        if min(fx, 1 - fx) >= 0.05: C[(X, R, "X_MAF>=0.05")] += 1
        for allele, ax, ar in (("ALT", fx, fr_), ("REF", 1 - fx, 1 - fr_)):
            if ax >= 0.05 and ar < 0.01:
                C[(X, R, "X>=0.05_R<0.01")] += 1
                if ar == 0: C[(X, R, "X>=0.05_R==0")] += 1
                if ax >= 0.20: C[(X, R, "X>=0.20_R<0.01")] += 1
                if tag == "sv":
                    L = abs(len(ref) - len(alt)); t = "DEL" if len(ref) > len(alt) else "INS"
                    if allele == "REF": t = {"DEL": "INS", "INS": "DEL"}[t] + "(ref-allele)"
                    lenbins[(X, R, t)] += 1
                if X == "NONC" and R == "COMM":
                    g = lambda k: "%.4f" % ((fr[idx[k]] if allele == "ALT" else 1 - fr[idx[k]]) if an[idx[k]] else float("nan"))
                    out.write("\t".join([chrom, pos, ref if tag != "sv" else str(len(ref)), alt if tag != "sv" else str(len(alt)),
                                         allele, g("NONC"), g("COMM"), g("AFR"), g("IDB"), g("SEAB"), g("HHG")]) + "\n")
out.close()
with open(outp + ".counts.tsv", "w") as h:
    h.write("chrom\ttest\tref\tmetric\tn\n")
    h.write(f"{chrom}\t-\t-\trecords\t{nsite}\n")
    for k, n in sorted(C.items()):
        if len(k) == 1: h.write(f"{chrom}\tALL\t-\t{k[0]}\t{n}\n")
        else: h.write(f"{chrom}\t{k[0]}\t{k[1]}\t{k[2]}\t{n}\n")
    for k, n in sorted(lenbins.items()):
        h.write(f"{chrom}\t{k[0]}\t{k[1]}\ttype_{k[2]}\t{n}\n")
    h.write(f"{chrom}\t-\t-\tnpresent\t" + ";".join(f"{g}={npres[g]}" for g in GROUPS) + "\n")
