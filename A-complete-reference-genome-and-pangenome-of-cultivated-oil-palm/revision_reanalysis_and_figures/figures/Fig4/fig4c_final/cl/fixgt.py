#!/usr/bin/env python3
"""stdin VCF (GT-only) -> stdout; haploid missing '.' -> './.'; counts written to argv[1]."""
import sys
nd = 0; nsite = 0; other_hap = 0
out = sys.stdout.buffer
for line in sys.stdin.buffer:
    if line[:1] == b'#':
        out.write(line); continue
    nsite += 1
    f = line.rstrip(b'\n').split(b'\t')
    g = f[9:]
    changed = False
    for i, x in enumerate(g):
        if len(x) == 1:
            if x == b'.':
                g[i] = b'./.'; nd += 1; changed = True
            else:
                other_hap += 1
    if changed:
        line = b'\t'.join(f[:9] + g) + b'\n'
    out.write(line)
with open(sys.argv[1], 'w') as h:
    h.write('sites\t%d\nhaploid_dot_converted\t%d\nother_haploid_calls\t%d\n' % (nsite, nd, other_hap))
