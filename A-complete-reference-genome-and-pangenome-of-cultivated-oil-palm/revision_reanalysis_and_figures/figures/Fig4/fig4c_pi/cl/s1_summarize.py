#!/usr/bin/env python3
"""Merge per-chromosome vcftools outputs and recompute genome-wide pi / FST exactly as 03_summarize_and_plot.py:
pi = mean of window PI weighted by N_VARIANTS+N_MONOMORPHIC; FST = mean of window MEAN_FST weighted by N_VARIANTS."""
import sys, glob, os, collections
base = sys.argv[1]          # dir with orig/, out/fixed, out/raw
ORIG = os.path.join(base, 'orig')
POPS = ['K4_Pop1', 'K4_Pop2', 'K4_Pop3', 'K4_Pop4']
PAIRS = [('K4_Pop1','K4_Pop2'),('K4_Pop1','K4_Pop3'),('K4_Pop1','K4_Pop4'),('K4_Pop2','K4_Pop3'),('K4_Pop2','K4_Pop4'),('K4_Pop3','K4_Pop4')]

def read(path):
    rows = []
    with open(path) as h:
        hdr = h.readline().split()
        for l in h:
            f = l.split()
            if f: rows.append(dict(zip(hdr, f)))
    return rows

def pi_stats(rows, chroms=None):
    num = den = 0.0; n = 0
    for r in rows:
        if chroms and r['CHROM'] not in chroms: continue
        w = (float(r['N_VARIANTS']) + float(r['N_MONOMORPHIC'])) if 'N_MONOMORPHIC' in r else (int(r['BIN_END']) - int(r['BIN_START']) + 1)
        v = r['PI']
        if v in ('nan','-nan','NA') or w <= 0: continue
        num += float(v) * w; den += w; n += 1
    return (num/den if den else float('nan')), n

def fst_stats(rows, chroms=None):
    num = den = 0.0; n = 0
    for r in rows:
        if chroms and r['CHROM'] not in chroms: continue
        w = float(r['N_VARIANTS']); v = r['MEAN_FST']
        if v in ('nan','-nan','NA') or w <= 0: continue
        num += float(v) * w; den += w; n += 1
    return (num/den if den else float('nan')), n

def load(setname, name, kind):
    rows = []
    suf = '.windowed.pi' if kind == 'pi' else '.windowed.weir.fst'
    tag = name + ('_100kb.pi' if kind == 'pi' else '_100kb_fst')
    for p in sorted(glob.glob(os.path.join(base, 'out', setname, '*__' + tag + suf))):
        rows += read(p)
    return rows

def orig(name, kind):
    return read(os.path.join(ORIG, name + ('_100kb.pi.windowed.pi' if kind=='pi' else '_100kb_fst.windowed.weir.fst')))

def chrom_cov(rows):
    c = collections.Counter(r['CHROM'] for r in rows)
    return c

out = []
# 1. validation: raw (unfixed) re-run vs original windows on chr14B, chr16B
print('## Validation: unfixed GT-only rerun (vcftools 0.1.16) vs original windows (0.1.17)')
for kind, names in (('pi', POPS), ('fst', ['%s_%s' % p for p in PAIRS])):
    for nm in names:
        o = orig(nm, kind); r = load('raw', nm, kind)
        for c in ('chr14B', 'chr16B'):
            oc = {(x['BIN_START']): x for x in o if x['CHROM']==c}
            rc = {(x['BIN_START']): x for x in r if x['CHROM']==c}
            col = 'PI' if kind=='pi' else 'MEAN_FST'
            same = sum(1 for k in oc if k in rc and abs(float(oc[k][col]) - float(rc[k][col])) < 1e-9*max(1,abs(float(oc[k][col]))) + 1e-12 and oc[k]['N_VARIANTS']==rc[k]['N_VARIANTS'])
            print('%s\t%s\t%s\torig_windows=%d\trerun_windows=%d\tidentical=%d' % (kind, nm, c, len(oc), len(rc), same))

print('\n## Genome-wide values')
print('metric\tgroup\tfigure_value(orig)\torig_nwin\tfixed_value\tfixed_nwin\tabs_diff\trel_diff_%\tfixed_same_chroms_as_orig\tfixed_14common_chroms')
common = None
for p in POPS:
    s = set(r['CHROM'] for r in orig(p, 'pi'))
    common = s if common is None else common & s
for kind, names in (('pi', POPS), ('fst', ['%s_%s' % p for p in PAIRS])):
    f = pi_stats if kind=='pi' else fst_stats
    for nm in names:
        o = orig(nm, kind); x = load('fixed', nm, kind)
        ov, on = f(o); xv, xn = f(x)
        oc = set(r['CHROM'] for r in o)
        sv, _ = f(x, oc); cv, _ = f(x, common)
        print('%s\t%s\t%.6g\t%d\t%.6g\t%d\t%+.6g\t%+.2f\t%.6g\t%.6g' % (kind, nm, ov, on, xv, xn, xv-ov, 100*(xv-ov)/ov, sv, cv))

print('\n## Chromosome coverage (windows; max BIN_END Mb) fixed vs orig')
for kind, names in (('pi', POPS), ('fst', ['%s_%s' % p for p in PAIRS])):
    for nm in names:
        o = chrom_cov(orig(nm, kind)); x = load('fixed', nm, kind); xc = chrom_cov(x)
        mx = collections.defaultdict(int)
        for r in x: mx[r['CHROM']] = max(mx[r['CHROM']], int(r['BIN_END']))
        diffs = ['%s:%d->%d' % (c, o.get(c,0), xc.get(c,0)) for c in sorted(set(o)|set(xc)) if o.get(c,0)!=xc.get(c,0)]
        print('%s\t%s\tchroms_orig=%d\tchroms_fixed=%d\tchr05B_end=%.1fMb\tchr10B_end=%.1fMb\tchanged: %s' % (kind, nm, len(o), len(xc), mx['chr05B']/1e6, mx['chr10B']/1e6, ' '.join(diffs)))
