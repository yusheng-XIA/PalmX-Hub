#!/usr/bin/env python3
"""Whole-chromosome tidk-style scan (10-kb windows anchored at 0); prints windows with
max(TTTAGGG, CCCTAAA) >= 20, to locate telomere arrays that are internal or >50 kb from an end.
usage: telo_scan.py manifest.tsv out.tsv"""
import sys
sys.path.insert(0, '.')
from telo_ends import read_fai, fetch, cnt, CHR_RE, W

man, outp = sys.argv[1], sys.argv[2]
out = open(outp, 'w')
out.write('label\tversion\tseq\tlength\twin_start\twin_end\tF\tR\tdist_to_nearest_end\n')
for line in open(man):
    if not line.strip() or line.startswith('#'):
        continue
    label, ver, fa = line.rstrip('\n').split('\t')
    chroms = [e for e in read_fai(fa) if CHR_RE.match(e[0]) and e[1] > 20_000_000][:16]
    with open(fa, 'rb') as fh:
        for ent in chroms:
            L = ent[1]
            seq = fetch(fh, ent, 0, L)
            for s in range(0, L, W):
                f, r = cnt(seq[s:s + W])
                if max(f, r) >= 20:
                    e = min(s + W, L)
                    out.write(f'{label}\t{ver}\t{ent[0]}\t{L}\t{s}\t{e}\t{f}\t{r}\t{min(s, L - e)}\n')
            del seq
    out.flush()
    print(label, ver, flush=True)
