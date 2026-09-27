#!/usr/bin/env python3
"""Parse SyRI output of the two E. oleifera haplotypes (MZ4 hap1/hap2 vs Africa_hap2)
into reference-coordinate aligned blocks, SNP calls, insertion and deletion records.
Output: one pickle per haplotype in W/mz4/."""
import pickle, sys, os
import numpy as np
from collections import defaultdict

W = 'work/outgroup_sensitivity'
M = 'sv_calls/syri'
ALN = {'SYNAL', 'INVAL', 'TRANSAL', 'INVTRAL', 'DUPAL', 'INVDPAL'}
os.makedirs(W + '/mz4', exist_ok=True)


def chrom_map(c):
    # SyRI reference chromosome 'chr01' -> 'chr01B'
    return c if c.endswith('B') else c + 'B'


def merge(iv):
    iv.sort()
    out = []
    for s, e in iv:
        if out and s <= out[-1][1]:
            if e > out[-1][1]:
                out[-1][1] = e
        else:
            out.append([s, e])
    return np.array(out, dtype=np.int64).reshape(-1, 2)


for hap in sys.argv[1:] or ['meizhou4_hap1', 'meizhou4_hap2']:
    blocks = defaultdict(list)          # 0-based half-open on reference
    snp = defaultdict(dict)             # pos(1-based) -> alt base
    ins = defaultdict(list)             # (pos, len)
    dels = defaultdict(list)            # (start0, end0)
    n = 0
    with open(f'{M}/{hap}/syri/syri.out') as fh:
        for line in fh:
            f = line.rstrip('\n').split('\t')
            t = f[10]
            if t in ALN:
                c = chrom_map(f[0]); s = int(f[1]) - 1; e = int(f[2])
                blocks[c].append((min(s, e), max(s, e)))
            elif t == 'SNP':
                snp[chrom_map(f[0])][int(f[1])] = f[4].upper()
            elif t == 'INS':
                ins[chrom_map(f[0])].append((int(f[1]), len(f[4]) - len(f[3])))
            elif t == 'DEL':
                c = chrom_map(f[0]); s = int(f[1]); e = int(f[2])
                # SyRI DEL: ref positions f[1]..f[2]; first base is the anchor
                dels[c].append((s, e))
            n += 1
    res = {'blocks': {c: merge(v) for c, v in blocks.items()},
           'snp': dict(snp),
           'ins': {c: np.array(sorted(v), dtype=np.int64).reshape(-1, 2) for c, v in ins.items()},
           'del': {c: merge([list(x) for x in v]) for c, v in dels.items()}}
    tot = sum(int((b[:, 1] - b[:, 0]).sum()) for b in res['blocks'].values())
    print(hap, 'lines', n, 'aligned_bp(merged)', tot, 'snps', sum(len(v) for v in snp.values()),
          'ins', sum(len(v) for v in ins.values()), 'del', sum(len(v) for v in dels.values()), flush=True)
    with open(f'{W}/mz4/{hap}.pkl', 'wb') as out:
        pickle.dump(res, out, protocol=4)
