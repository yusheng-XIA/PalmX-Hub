#!/usr/bin/env python3
# Convert minimap2 SAM (EG11 ref, FL query) to PAF-like blocks and compare with Source Data synteny sheet
import re, csv, sys, collections
SAM = sys.argv[1]
W = '${CLUSTER_WORK}/trace/Fig1c'
def norm(c):
    m = re.search(r'chr(\d+)', c); return 'chr%02d' % int(m.group(1)) if m else c
rows = list(csv.reader(open(W + '/sd/Fig1c_synteny.tsv'), delimiter='\t'))[1:]
sdset = set((r[0], int(float(r[1])), int(float(r[2])), r[3], int(float(r[4])), int(float(r[5])), r[6], int(float(r[7])), int(float(r[8]))) for r in rows)
qlen = {}
cig = re.compile(r'(\d+)([MIDNSHP=X])')
out = open(W + '/out/sam_blocks.tsv', 'w')
allb = set(); n = 0
for l in open(SAM):
    if l[0] == '@':
        continue
    f = l.split('\t', 11)
    flag = int(f[1])
    if flag & 4: continue
    rname, pos, mq, cg = f[2], int(f[3]) - 1, int(f[4]), f[5]
    ops = cig.findall(cg)
    lead = 0; trail = 0; qa = 0; ra = 0; blk = 0
    for i, (k, o) in enumerate(ops):
        k = int(k)
        if o in 'SH':
            if i == 0: lead = k
            else: trail = k
        elif o in 'M=X': qa += k; ra += k; blk += k
        elif o == 'I': qa += k; blk += k
        elif o in 'DN': ra += k; blk += k
    qtot = lead + qa + trail
    strand = '-' if flag & 16 else '+'
    if strand == '+': qs = lead
    else: qs = trail
    qe = qs + qa
    rec = (norm(rname), pos, pos + ra, norm(f[0]), qs, qe, strand, blk, mq)
    n += 1
    out.write('\t'.join(map(str, rec)) + '\n')
    if rec[0] == rec[3] and mq >= 20 and blk >= 50000: allb.add(rec)
out.close()
print('records', n, 'filtered', len(allb), 'SD', len(sdset), 'exact', len(allb & sdset))
print('only_SD', sorted(sdset - allb)[:5]); print('only_mine', sorted(allb - sdset)[:5])
