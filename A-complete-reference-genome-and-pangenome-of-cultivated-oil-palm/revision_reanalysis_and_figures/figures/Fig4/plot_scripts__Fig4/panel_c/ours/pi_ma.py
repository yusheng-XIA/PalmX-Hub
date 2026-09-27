#!/usr/bin/env python3
"""Window pi per K4 population, 100-kb windows / 50-kb step (vcftools bins: 1, 50001, ...).
Writes, per window with >=1 site polymorphic in the population:
 CHROM BIN_START BIN_END N_VARIANTS PI_VCFTOOLS PI_MISSING_AWARE N_SITES_CALLED MEAN_CALLED_ALLELES
 PI_VCFTOOLS    = sum_site k(n-k)*2 / (L*N(N-1) - sum_polysite (N(N-1)-n(n-1)))  [reproduces vcftools 0.1.16 --window-pi exactly]
 PI_MISSING_AWARE = sum_site 2k(n-k)/(n(n-1)) / L   [per-site estimator on called alleles; == vcftools when no missing GT]
usage: pi_ma.py in.vcf.gz groups_dir out_prefix"""
import sys, gzip, collections
vcf, gdir, outp = sys.argv[1:4]
L = 100000; STEP = 50000
pops = ['K4_Pop1', 'K4_Pop2', 'K4_Pop3', 'K4_Pop4']
mem = {p: set(l.strip() for l in open('%s/%s.txt' % (gdir, p)) if l.strip()) for p in pops}
acc = {p: collections.defaultdict(lambda: [0, 0.0, 0.0, 0.0, 0, 0]) for p in pops}  # npoly, mism, adj, pi_ma, ncalled_sites, sum_n
chrom = None
with (gzip.open(vcf, 'rt') if vcf.endswith('.gz') else open(vcf)) as h:
    for line in h:
        if line[0] == '#':
            if line.startswith('#CHROM'):
                names = line.rstrip('\n').split('\t')[9:]
                idx = {p: [i + 9 for i, s in enumerate(names) if s in mem[p]] for p in pops}
                NN = {p: 2 * len(idx[p]) for p in pops}
            continue
        f = line.rstrip('\n').split('\t')
        chrom = f[0]; pos = int(f[1])
        j = (pos - 1) // STEP
        bins = [b for b in (j, j - 1) if b >= 0]
        for p in pops:
            g = ''.join([f[i] for i in idx[p]])
            k = g.count('1'); r = g.count('0'); n = k + r
            if n < 2: continue
            N = NN[p]
            for b in bins:
                a = acc[p][b]; a[4] += 1; a[5] += n
            if k == 0 or r == 0: continue
            mm = 2.0 * k * r
            for b in bins:
                a = acc[p][b]
                a[0] += 1; a[1] += mm; a[2] += N * (N - 1) - n * (n - 1); a[3] += mm / (n * (n - 1))
for p in pops:
    N = NN[p]
    with open('%s__%s_100kb.pima.tsv' % (outp, p), 'w') as o:
        o.write('CHROM\tBIN_START\tBIN_END\tN_VARIANTS\tPI_VCFTOOLS\tPI_MISSING_AWARE\tN_SITES_CALLED\tMEAN_CALLED_ALLELES\n')
        for b in sorted(acc[p]):
            a = acc[p][b]
            if a[0] == 0: continue
            s = b * STEP + 1
            o.write('%s\t%d\t%d\t%d\t%.6g\t%.6g\t%d\t%.2f\n' % (chrom, s, s + L - 1, a[0], a[1] / (L * N * (N - 1) - a[2]), a[3] / L, a[4], a[5] / a[4]))
