#!/usr/bin/env python3
"""stdin: 308-sample GT-only VCF already filtered QUAL>=30 & MAF>=0.05 (MAF on called alleles of the 308).
Converts haploid '.' -> './.'. Writes:
  strict  : F_MISSING (over all 308) <= 0.2   == original vcftools --max-missing 0.8
  present : F_MISSING computed only over samples NOT listed as 'absent' for this chromosome <= 0.2
            (absent = legacy sample/chromosome pairs whose genotypes are entirely '.' in the final callset)
argv: chrom absent_list(sample<TAB>chrom) out_strict out_present stats_file"""
import sys, subprocess
chrom, absfile, o1, o2, statf = sys.argv[1:6]
absent = set()
for l in open(absfile):
    f = l.split()
    if len(f) >= 2 and f[1] == chrom: absent.add(f[0])
BGZ = 'bgzip'
p1 = subprocess.Popen([BGZ, '-c'], stdin=subprocess.PIPE, stdout=open(o1, 'wb'))
p2 = subprocess.Popen([BGZ, '-c'], stdin=subprocess.PIPE, stdout=open(o2, 'wb'))
w1, w2 = p1.stdin, p2.stdin
idx_present = None; n = n1 = n2 = 0; ndot = 0; dot_abs = 0; nonabs_called_on_abs = 0
for line in sys.stdin.buffer:
    if line[:1] == b'#':
        if line.startswith(b'#CHROM'):
            s = line.rstrip(b'\n').split(b'\t')[9:]
            names = [x.decode() for x in s]
            idx_present = [i for i, x in enumerate(names) if x not in absent]
            idx_abs = [i for i, x in enumerate(names) if x in absent]
            nabs_found = len(idx_abs)
        w1.write(line); w2.write(line); continue
    n += 1
    f = line.rstrip(b'\n').split(b'\t')
    g = f[9:]
    miss = [False] * len(g)
    for i, x in enumerate(g):
        if x[:1] == b'.':
            miss[i] = True
            if len(x) == 1:
                g[i] = b'./.'; ndot += 1
    for i in idx_abs:
        if not miss[i]: nonabs_called_on_abs += 1
    m_all = sum(miss) / len(g)
    m_pr = sum(miss[i] for i in idx_present) / len(idx_present)
    out = b'\t'.join(f[:9] + g) + b'\n'
    if m_all <= 0.2 + 1e-12:
        w1.write(out); n1 += 1
    if m_pr <= 0.2 + 1e-12:
        w2.write(out); n2 += 1
w1.close(); w2.close(); p1.wait(); p2.wait()
with open(statf, 'w') as h:
    h.write('chrom\t%s\nabsent_samples_in_308\t%d\nsites_maf_qual\t%d\nsites_strict\t%d\nsites_present\t%d\nhaploid_dot_converted\t%d\ncalled_GT_in_absent_blocks\t%d\n'
            % (chrom, nabs_found, n, n1, n2, ndot, nonabs_called_on_abs))
