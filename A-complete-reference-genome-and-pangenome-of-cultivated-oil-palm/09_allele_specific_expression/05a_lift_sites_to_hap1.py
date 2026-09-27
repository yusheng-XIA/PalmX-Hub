#!/usr/bin/env python3
"""Lift FL diagnostic ASE sites from FL-Hap2 (allele A, backbone) to FL-Hap1 (allele B) coordinates.

Each site's 201-bp FL-Hap2 flank, with the centre base replaced by the FL-Hap1 allele, is aligned to FL-Hap1
with BWA-MEM; the aligned position of the centre base gives the FL-Hap1 coordinate (MAPQ >= 20 used downstream).

usage: 05a_lift_sites_to_hap1.py SITES.tsv.gz FL_Hap2.fa FL_Hap1.fa OUTDIR
  SITES: diagnostic sites with columns chrom, pos, ref, alt, gene_id, retained (from the marker QC step)
"""
import subprocess, re, sys, os
import pandas as pd
MQ, FA_A, REF_B, W = sys.argv[1:5]
os.makedirs(W, exist_ok=True)
BIN = ""  # samtools / bwa on PATH
s = pd.read_csv(MQ, sep="\t", usecols=["chrom","pos","ref","alt","gene_id","retained"])
s = s.sort_values(["chrom","pos"]).reset_index(drop=True)
s["sid"] = s.chrom + ":" + s.pos.astype(str)
s[["chrom","pos"]].to_csv(f"{W}/FL_A_positions.tsv", sep="\t", header=False, index=False)
s.to_csv(f"{W}/FL_sites.tsv.gz", sep="\t", index=False, compression="gzip")
F = 100
with open(f"{W}/regions.txt", "w") as f:
    for c, p in zip(s.chrom, s.pos): f.write(f"{c}:{max(1,p-F)}-{p+F}\n")
fa = subprocess.run([BIN+"samtools", "faidx", FA_A, "-r", f"{W}/regions.txt"], capture_output=True, text=True, check=True).stdout
seqs = []; cur = None; buf = []
for line in fa.splitlines():
    if line.startswith(">"):
        if cur: seqs.append((cur, "".join(buf)))
        cur = line[1:]; buf = []
    else: buf.append(line.strip())
seqs.append((cur, "".join(buf)))
assert len(seqs) == len(s)
bad = 0
with open(f"{W}/flanks_Ballele.fa", "w") as f:
    for (name, seq), r in zip(seqs, s.itertuples()):
        off = r.pos - max(1, r.pos - F)
        if seq[off].upper() != r.ref: bad += 1
        q = seq[:off] + r.alt + seq[off+1:]
        f.write(f">{r.sid}|{off}\n{q}\n")
print("ref-base mismatches vs FL-Hap2:", bad, "of", len(s))
cmd = f"{BIN}bwa mem -t 12 -v 1 {REF_B} {W}/flanks_Ballele.fa | {BIN}samtools view -F 0x900 - > {W}/flanks_on_B.sam"
subprocess.run(cmd, shell=True, check=True)
comp = str.maketrans("ACGT", "TGCA")
out = []
for line in open(f"{W}/flanks_on_B.sam"):
    t = line.split("\t"); flag = int(t[1])
    if flag & 4: continue
    sid, off = t[0].split("|"); off = int(off); mapq = int(t[4]); rpos = int(t[3]); cig = t[5]; seq = t[9]
    nm = next((int(x[5:]) for x in t[11:] if x.startswith("NM:i:")), -1)
    rev = bool(flag & 16)
    qoff = (len(seq) - 1 - off) if rev else off     # SAM SEQ is revcomp for reverse; offset in SAM seq coords
    qi = 0; ri = rpos; hit = None
    for n, op in re.findall(r"(\d+)([MIDNSHP=X])", cig):
        n = int(n)
        if op in "M=X":
            if qi <= qoff < qi + n: hit = ri + (qoff - qi); break
            qi += n; ri += n
        elif op in "IS": 
            if qi <= qoff < qi + n: break
            qi += n
        elif op in "DN": ri += n
    if hit is None: continue
    out.append((sid, t[2], hit, "-" if rev else "+", mapq, nm, seq[qoff]))
b = pd.DataFrame(out, columns=["sid","chrom_B","pos_B","strand_B","mapq_B","nm_B","base_in_query"])
m = s.merge(b, on="sid", how="left")
# expected B-reference base at pos_B: alt (forward) or complement(alt) (reverse)
m["alt_onB"] = [a if st == "+" else a.translate(comp) for a, st in zip(m.alt, m.strand_B.fillna("+"))]
m["ref_onB"] = [a if st == "+" else a.translate(comp) for a, st in zip(m.ref, m.strand_B.fillna("+"))]
m.to_csv(f"{W}/FL_sites_liftB.tsv.gz", sep="\t", index=False, compression="gzip")
ok = m[(m.mapq_B >= 20)]
ok[["chrom_B","pos_B"]].astype({"pos_B": int}).drop_duplicates().sort_values(["chrom_B","pos_B"]).to_csv(f"{W}/FL_B_positions.tsv", sep="\t", header=False, index=False)
print("sites", len(s), "lifted", m.pos_B.notna().sum(), "mapq>=20", len(ok), "chrom A->B same number", (ok.chrom.str[:5] == ok.chrom_B.str[:5]).mean())
