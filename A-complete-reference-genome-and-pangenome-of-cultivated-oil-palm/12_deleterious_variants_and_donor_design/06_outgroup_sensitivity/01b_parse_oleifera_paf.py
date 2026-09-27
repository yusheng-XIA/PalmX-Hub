#!/usr/bin/env python3
"""Convert minimap2 asm20 --cs PAF (E. oleifera haplotype vs Africa_hap2, MAPQ >= 20, the
same preset/threshold as the Phoenix dSNP polarization) into the pickle format of
enhB_01: merged aligned target blocks, SNP bases (only at dSNP candidate sites),
query insertions (target pos, len) and target bases deleted in the query."""
import csv, pickle, re, sys, os
import numpy as np
from collections import defaultdict

W = 'work/outgroup_sensitivity'
hap, paf = sys.argv[1], sys.argv[2]
os.makedirs(W + '/mz4_mm2', exist_ok=True)
TOK = re.compile(r'(:\d+|=[A-Za-z]+|\*[A-Za-z][A-Za-z]|\+[A-Za-z]+|-[A-Za-z]+|~[A-Za-z]{2}\d+[A-Za-z]{2})')
want = defaultdict(set)
with open(W + '/dsnp_candidates_eg35.syri.tsv') as fh:
    for r in csv.DictReader(fh, delimiter='\t'):
        want[r['Chrom']].add(int(r['Pos']))


def merge(iv):
    iv.sort(); out = []
    for s, e in iv:
        if out and s <= out[-1][1]:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([s, e])
    return np.array(out, dtype=np.int64).reshape(-1, 2)


blocks = defaultdict(list); dels = defaultdict(list); ins = defaultdict(list); snp = defaultdict(dict)
n = used = 0
with open(paf) as fh:
    for line in fh:
        f = line.rstrip('\n').split('\t'); n += 1
        if int(f[11]) < 20:
            continue
        ch = f[5]; ts = int(f[7]); te = int(f[8]); used += 1
        blocks[ch].append([ts, te])
        cs = next((x[5:] for x in f[12:] if x.startswith('cs:Z:')), None)
        if cs is None:
            continue
        t = ts; w = want.get(ch, ())
        for tok in TOK.findall(cs):
            op = tok[0]
            if op == ':':
                t += int(tok[1:])
            elif op == '=':
                t += len(tok) - 1
            elif op == '*':
                if (t + 1) in w:
                    snp[ch][t + 1] = tok[2].upper()
                t += 1
            elif op == '-':
                L = len(tok) - 1; dels[ch].append([t, t + L]); t += L
            elif op == '+':
                ins[ch].append((t, len(tok) - 1))
            elif op == '~':
                m = re.match(r'~[A-Za-z]{2}(\d+)[A-Za-z]{2}', tok); t += int(m.group(1))
res = {'blocks': {c: merge(v) for c, v in blocks.items()}, 'snp': dict(snp),
       'ins': {c: np.array(sorted(v), dtype=np.int64).reshape(-1, 2) for c, v in ins.items()},
       'del': {c: merge(v) for c, v in dels.items()}}
tot = sum(int((b[:, 1] - b[:, 0]).sum()) for b in res['blocks'].values())
print(hap, 'paf lines', n, 'mapq20', used, 'aligned_bp(merged)', tot, flush=True)
pickle.dump(res, open(f'{W}/mz4_mm2/{hap}.pkl', 'wb'), protocol=4)
