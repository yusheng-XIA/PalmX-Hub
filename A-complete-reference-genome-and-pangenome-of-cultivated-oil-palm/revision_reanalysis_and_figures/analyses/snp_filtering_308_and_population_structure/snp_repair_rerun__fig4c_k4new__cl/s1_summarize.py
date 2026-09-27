#!/usr/bin/env python3
"""Summarize per-chromosome vcftools windows exactly as 03_summarize_and_plot.py:
pi  = mean of window PI weighted by N_VARIANTS+N_MONOMORPHIC (vcftools reports N_MONO only with invariant sites;
      when absent the original script's weights are None -> fall back to equal weights, i.e. simple mean; we report both)
FST = mean of window MEAN_FST weighted by N_VARIANTS (W&C mean), plus reference N_VARIANTS-weighted WEIGHTED_FST."""
import glob, os, collections, csv, sys
B = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(B, 'out'); ORIG = os.path.join(B, '..', 'fig4c_pi', 'orig'); PREV = os.path.join(B, '..', 'fig4c_pi', 'out', 'fixed')
POPS = ['K4_Pop1', 'K4_Pop2', 'K4_Pop3', 'K4_Pop4']
PAIRS = [('K4_Pop1','K4_Pop2'),('K4_Pop1','K4_Pop3'),('K4_Pop1','K4_Pop4'),('K4_Pop2','K4_Pop3'),('K4_Pop2','K4_Pop4'),('K4_Pop3','K4_Pop4')]
CHR = [l.split()[0] for l in open(os.path.join(B, 'meta', 'chrom_len.tsv'))]
LEN = {l.split()[0]: int(l.split()[1]) for l in open(os.path.join(B, 'meta', 'chrom_len.tsv'))}
FIG = {'K4_Pop1':0.0056650237166,'K4_Pop2':0.00470491658359,'K4_Pop3':0.00577306769722,'K4_Pop4':0.00597151287293,
 'K4_Pop1_K4_Pop2':0.126357843235,'K4_Pop1_K4_Pop3':0.0657473265601,'K4_Pop1_K4_Pop4':0.0436513270168,'K4_Pop2_K4_Pop3':0.0928941177597,'K4_Pop2_K4_Pop4':0.129326048602,'K4_Pop3_K4_Pop4':0.0819804232249}
def rd(p):
    with open(p) as h:
        hd = h.readline().split(); return [dict(zip(hd, l.split())) for l in h if l.strip()]
def fl(x):
    try:
        v = float(x); return None if v != v else v
    except: return None
def load(d, nm, kind):
    suf = '_100kb.pi.windowed.pi' if kind == 'pi' else '_100kb_fst.windowed.weir.fst'
    rows = []
    for c in CHR:
        p = os.path.join(d, '%s__%s%s' % (c, nm, suf))
        if os.path.exists(p): rows += rd(p)
    return rows
def pi_sum(rows):
    num = den = 0.0; n = 0; nv = 0
    for r in rows:
        v = fl(r['PI']);
        if v is None: continue
        w = (float(r['N_VARIANTS']) + float(r['N_MONOMORPHIC'])) if 'N_MONOMORPHIC' in r else 1.0
        num += v * w; den += w; n += 1
        if (int(r['BIN_START']) - 1) % 100000 == 0: nv += int(r['N_VARIANTS'])
    return num / den, n, nv
def fst_sum(rows):
    num = den = numw = 0.0; n = 0; nv = 0
    for r in rows:
        v = fl(r['MEAN_FST']); w = float(r['N_VARIANTS'])
        if v is None or w <= 0: continue
        num += v * w; den += w; n += 1; numw += (fl(r['WEIGHTED_FST']) or 0) * w
        if (int(r['BIN_START']) - 1) % 100000 == 0: nv += int(r['N_VARIANTS'])
    return num / den, n, nv, numw / den
def prev_val(nm, kind):
    rows = []
    suf = '_100kb.pi.windowed.pi' if kind == 'pi' else '_100kb_fst.windowed.weir.fst'
    for p in glob.glob(os.path.join(PREV, '*__' + nm + suf)): rows += rd(p)
    if not rows: return None
    return (pi_sum(rows) if kind == 'pi' else fst_sum(rows))[0]
sites = {}
for s in ('strict', 'present'):
    sites[s] = {}
    for c in CHR:
        p = os.path.join(B, 'logs', 'split_%s.txt' % c)
        d = dict(l.rstrip('\n').split('\t') for l in open(p))
        sites[s][c] = int(d['sites_' + s])
res = []
for s in ('strict', 'present'):
    d = os.path.join(OUT, s)
    for nm in POPS:
        v, n, nv = pi_sum(load(d, nm, 'pi'))
        res.append([s, 'pi', nm, v, n, nv, sum(sites[s].values()), '', FIG[nm], prev_val(nm, 'pi')])
    for a, b in PAIRS:
        nm = a + '_' + b
        v, n, nv, wf = fst_sum(load(d, nm, 'fst'))
        res.append([s, 'fst_mean', nm, v, n, nv, sum(sites[s].values()), wf, FIG[nm], prev_val(nm, 'fst')])
with open(os.path.join(B, 'pi_fst_genomewide.tsv'), 'w') as h:
    h.write('snp_set\tmetric\tgroup_or_pair\tvalue\tn_windows\tn_SNPs_in_windows(non-overlapping 100kb bins)\tn_SNPs_in_set_16chr\tref_WEIGHTED_FST_nvar_weighted\tfigure_value\tprevious_fixed_value(old VCF, dot->./.)\tdiff_vs_figure_%\tdiff_vs_previous_%\n')
    for r in res:
        dv = 100 * (r[3] / r[8] - 1); dp = 100 * (r[3] / r[9] - 1) if r[9] else float('nan')
        h.write('\t'.join([r[0], r[1], r[2], '%.6g' % r[3], str(r[4]), str(r[5]), str(r[6]), ('%.6g' % r[7]) if r[7] != '' else '', '%.6g' % r[8], '%.6g' % r[9] if r[9] else 'NA', '%+.2f' % dv, '%+.2f' % dp]) + '\n')
# coverage per chromosome
with open(os.path.join(B, 'coverage_by_chrom.tsv'), 'w') as h:
    h.write('snp_set\tgroup_or_pair\tchrom\tchrom_len_bp\tsnps_in_set\tn_windows\texpected_windows\tfirst_bin_start\tlast_bin_end\tlast_end_to_chrom_end_bp\tmax_gap_between_windows_bp\tchrom_value\n')
    for s in ('strict', 'present'):
        d = os.path.join(OUT, s)
        for nm, kind in [(p, 'pi') for p in POPS] + [(a + '_' + b, 'fst') for a, b in PAIRS]:
            rows = load(d, nm, kind)
            by = collections.defaultdict(list)
            for r in rows: by[r['CHROM']].append(r)
            for c in CHR:
                rr = sorted(by.get(c, []), key=lambda r: int(r['BIN_START']))
                exp = (LEN[c] - 1) // 50000 + 1 - 1
                if rr:
                    st = [int(r['BIN_START']) for r in rr]; gap = max([b - a for a, b in zip(st, st[1:])] + [0])
                    val = (pi_sum(rr) if kind == 'pi' else fst_sum(rr))[0]
                    h.write('\t'.join(map(str, [s, nm, c, LEN[c], sites[s][c], len(rr), exp, st[0], rr[-1]['BIN_END'], LEN[c] - int(rr[-1]['BIN_END']), gap, '%.6g' % val])) + '\n')
                else:
                    h.write('\t'.join(map(str, [s, nm, c, LEN[c], sites[s][c], 0, exp, 'NA', 'NA', 'NA', 'NA', 'NA'])) + '\n')
print(open(os.path.join(B, 'pi_fst_genomewide.tsv')).read())
# ---- same-chromosome comparison on the 11 chromosomes without absent sample blocks ----
CLEAN = ['chr01B','chr02B','chr03B','chr06B','chr08B','chr09B','chr11B','chr12B','chr13B','chr14B','chr15B']
def orig_rows(nm, kind):
    return rd(os.path.join(ORIG, nm + ('_100kb.pi.windowed.pi' if kind == 'pi' else '_100kb_fst.windowed.weir.fst')))
def prev_rows(nm, kind):
    suf = '_100kb.pi.windowed.pi' if kind == 'pi' else '_100kb_fst.windowed.weir.fst'
    r = []
    for p in glob.glob(os.path.join(PREV, '*__' + nm + suf)): r += rd(p)
    return r
with open(os.path.join(B, 'compare_11_clean_chroms.tsv'), 'w') as h:
    h.write('metric\tgroup_or_pair\torig_figure_VCF\tprev_fixed_oldVCF\tfinal_strict\tfinal_present\tfinal_strict_vs_orig_%\tnwin_orig\tnwin_final\n')
    for nm, kind in [(p, 'pi') for p in POPS] + [(a + '_' + b, 'fst') for a, b in PAIRS]:
        f = pi_sum if kind == 'pi' else fst_sum
        sel = lambda rows: [r for r in rows if r['CHROM'] in CLEAN]
        o = f(sel(orig_rows(nm, kind))); pv = f(sel(prev_rows(nm, kind)))
        st = f(sel(load(os.path.join(OUT, 'strict'), nm, kind))); pr = f(sel(load(os.path.join(OUT, 'present'), nm, kind)))
        h.write('%s\t%s\t%.6g\t%.6g\t%.6g\t%.6g\t%+.2f\t%d\t%d\n' % (kind, nm, o[0], pv[0], st[0], pr[0], 100 * (st[0] / o[0] - 1), o[1], st[1]))
print(open(os.path.join(B, 'compare_11_clean_chroms.tsv')).read())
