#!/usr/bin/env python3
"""Stream one byte range of the 349 GB oil_filter VCF and write GT-only per-chromosome body parts.
FIXED output: haploid missing GT '.' -> diploid missing './.' (the only change).
RAW output (validation chromosomes only): GT strings kept exactly as in the source.
Also counts per-sample haploid GTs per chromosome."""
import sys, os, collections
V, wid, nw, outdir = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), sys.argv[4]
RAWCHR = set(sys.argv[5].split(',')) if len(sys.argv) > 5 else set()
size = os.path.getsize(V)
start = size * wid // nw
end = size * (wid + 1) // nw
fh = open(V, 'rb', buffering=16 * 1024 * 1024)
fh.seek(start)
pos = start
if start > 0:
    fh.seek(start - 1)
    pos = start - 1
    l = fh.readline(); pos += len(l)   # skip to the first line that starts at >= start
outs = {}; raws = {}
nsite = collections.Counter(); hapc = collections.defaultdict(collections.Counter)
nonstd = collections.Counter()
while pos < end:
    line = fh.readline()
    if not line:
        break
    pos += len(line)
    if line[:1] == b'#':
        continue
    f = line.rstrip(b'\n').split(b'\t')
    c = f[0]
    gts = [x.split(b':', 1)[0] for x in f[9:]]
    nsite[c] += 1
    fixed = []
    for i, g in enumerate(gts):
        if len(g) == 1:                       # haploid call ('.', or '0'/'1')
            hapc[c][i] += 1
            if g == b'.':
                g = b'./.'
            else:
                nonstd[g] += 1
        fixed.append(g)
    base = b'\t'.join(f[:5]) + b'\t.\t.\t.\tGT\t'
    o = outs.get(c)
    if o is None:
        o = outs[c] = open('%s/part%02d_%s.body' % (outdir, wid, c.decode()), 'wb', buffering=8 * 1024 * 1024)
    o.write(base + b'\t'.join(fixed) + b'\n')
    if c.decode() in RAWCHR:
        r = raws.get(c)
        if r is None:
            r = raws[c] = open('%s/raw%02d_%s.body' % (outdir, wid, c.decode()), 'wb', buffering=8 * 1024 * 1024)
        r.write(base + b'\t'.join(gts) + b'\n')
for o in list(outs.values()) + list(raws.values()):
    o.close()
with open('%s/stats%02d.tsv' % (outdir, wid), 'w') as s:
    for c in nsite:
        s.write('SITES\t%s\t%d\n' % (c.decode(), nsite[c]))
    for c in hapc:
        for i, n in hapc[c].items():
            s.write('HAP\t%s\t%d\t%d\n' % (c.decode(), i, n))
    for g, n in nonstd.items():
        s.write('NONSTD\t%s\t%d\n' % (g.decode(), n))
    s.write('DONE\t%d\t%d\n' % (start, end))
