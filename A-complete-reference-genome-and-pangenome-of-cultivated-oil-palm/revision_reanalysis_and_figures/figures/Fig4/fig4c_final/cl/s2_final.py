#!/usr/bin/env python3
"""Final tables for Fig. 4c from the final (corrected joint call → 308) callset.
Recommended = 'present' SNP set (all 16 chr), pi = missing-aware per-site estimator, FST = vcftools W&C MEAN_FST N_VARIANTS-weighted."""
import glob, os, collections, csv
B = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
POPS = ['K4_Pop1', 'K4_Pop2', 'K4_Pop3', 'K4_Pop4']
PAIRS = [('K4_Pop1','K4_Pop2'),('K4_Pop1','K4_Pop3'),('K4_Pop1','K4_Pop4'),('K4_Pop2','K4_Pop3'),('K4_Pop2','K4_Pop4'),('K4_Pop3','K4_Pop4')]
CL = [l.split() for l in open(os.path.join(B, 'meta', 'chrom_len.tsv'))]
CHR = [c for c, _ in CL]; LEN = {c: int(n) for c, n in CL}
CLEAN = ['chr01B','chr02B','chr03B','chr06B','chr08B','chr09B','chr11B','chr12B','chr13B','chr14B','chr15B']
FIG = {'K4_Pop1':0.0056650237166,'K4_Pop2':0.00470491658359,'K4_Pop3':0.00577306769722,'K4_Pop4':0.00597151287293,
 'K4_Pop1_K4_Pop2':0.126357843235,'K4_Pop1_K4_Pop3':0.0657473265601,'K4_Pop1_K4_Pop4':0.0436513270168,'K4_Pop2_K4_Pop3':0.0928941177597,'K4_Pop2_K4_Pop4':0.129326048602,'K4_Pop3_K4_Pop4':0.0819804232249}
PREV = {'K4_Pop1':0.00554711,'K4_Pop2':0.00456683,'K4_Pop3':0.00569249,'K4_Pop4':0.00586287,
 'K4_Pop1_K4_Pop2':0.122687,'K4_Pop1_K4_Pop3':0.0633396,'K4_Pop1_K4_Pop4':0.042167,'K4_Pop2_K4_Pop3':0.0929228,'K4_Pop2_K4_Pop4':0.129178,'K4_Pop3_K4_Pop4':0.0807861}
def rd(p):
    with open(p) as h:
        hd = h.readline().split(); return [dict(zip(hd, l.split())) for l in h if l.strip()]
def pima(sub, nm, chroms=None):
    rows = []
    for c in CHR:
        p = os.path.join(B, 'out', sub, '%s__%s_100kb.pima.tsv' % (c, nm))
        if os.path.exists(p): rows += rd(p)
    if chroms: rows = [r for r in rows if r['CHROM'] in chroms]
    return rows
def fstrows(sub, nm, chroms=None):
    rows = []
    for c in CHR:
        p = os.path.join(B, 'out', sub, '%s__%s_100kb_fst.windowed.weir.fst' % (c, nm))
        if os.path.exists(p): rows += rd(p)
    if chroms: rows = [r for r in rows if r['CHROM'] in chroms]
    return rows
def mean(rows, col):
    v = [float(r[col]) for r in rows]; return sum(v) / len(v), len(v)
def fst(rows):
    num = den = numw = 0.0; n = 0
    for r in rows:
        if r['MEAN_FST'] in ('nan', '-nan'): continue
        w = float(r['N_VARIANTS']); num += float(r['MEAN_FST']) * w; numw += float(r['WEIGHTED_FST']) * w; den += w; n += 1
    return num / den, n, numw / den
def nsnp_nonoverlap(rows):
    return sum(int(r['N_VARIANTS']) for r in rows if (int(r['BIN_START']) - 1) % 100000 == 0)
sites = {s: {} for s in ('strict', 'present')}
for c in CHR:
    d = dict(l.rstrip('\n').split('\t') for l in open(os.path.join(B, 'logs', 'split_%s.txt' % c)))
    for s in sites: sites[s][c] = int(d['sites_' + s])
T = {s: sum(v.values()) for s, v in sites.items()}
# ---------- genome-wide table ----------
rows = []
for nm in POPS:
    for s in ('present', 'strict'):
        r = pima('pima_' + s, nm)
        ma, n = mean(r, 'PI_MISSING_AWARE'); vt, _ = mean(r, 'PI_VCFTOOLS')
        rows.append(dict(snp_set=s, metric='pi', estimator='missing-aware per-site (recommended)' , group=nm, value=ma, n_windows=n, n_chrom=len(set(x['CHROM'] for x in r)), snps_pop=nsnp_nonoverlap(r), snps_set=T[s]))
        rows.append(dict(snp_set=s, metric='pi', estimator='vcftools --window-pi (original method)', group=nm, value=vt, n_windows=n, n_chrom=len(set(x['CHROM'] for x in r)), snps_pop=nsnp_nonoverlap(r), snps_set=T[s]))
for a, b in PAIRS:
    nm = a + '_' + b
    for s in ('present', 'strict'):
        r = fstrows(s, nm); v, n, w = fst(r)
        rows.append(dict(snp_set=s, metric='fst', estimator='W&C mean FST, window MEAN_FST weighted by N_VARIANTS (original method)', group=nm, value=v, n_windows=n, n_chrom=len(set(x['CHROM'] for x in r)), snps_pop=nsnp_nonoverlap(r), snps_set=T[s], wfst=w))
with open(os.path.join(B, 'pi_fst_genomewide.tsv'), 'w') as h:
    h.write('metric\tgroup_or_pair\tsnp_set\testimator\tvalue\tn_windows\tn_chromosomes\tn_SNPs_used(polymorphic_in_pop/pair; non-overlapping bins)\tn_SNPs_in_set\tref_weighted_FST(not plotted)\tfigure_value\tprevious_fixed_value(old VCF)\tdiff_vs_figure_%\tdiff_vs_previous_%\n')
    for r in rows:
        f = FIG[r['group']]; p = PREV[r['group']]
        h.write('\t'.join(map(str, [r['metric'], r['group'], r['snp_set'], r['estimator'], '%.6g' % r['value'], r['n_windows'], r['n_chrom'], r['snps_pop'], r['snps_set'],
                '%.6g' % r['wfst'] if 'wfst' in r else '', '%.6g' % f, '%.6g' % p, '%+.2f' % (100 * (r['value'] / f - 1)), '%+.2f' % (100 * (r['value'] / p - 1))])) + '\n')
# ---------- coverage per population x chromosome (recommended set) ----------
with open(os.path.join(B, 'coverage_by_pop_chrom.tsv'), 'w') as h:
    h.write('group_or_pair\tchrom\tchrom_len_bp\tSNPs_in_present_set\tSNPs_in_strict_set\tabsent_samples_in_group(all-missing on this chrom)\tn_windows\tn_bins_expected\tfirst_bin_start\tlast_bin_end\tlast_bin_end_minus_chrom_len\tlargest_gap_between_window_starts_bp\tmean_called_alleles_per_site\tpi_missing_aware\tpi_vcftools\tFST_mean\n')
    absent = collections.defaultdict(set)
    for l in open(os.path.join(B, 'meta', 'absent_308.tsv')):
        s, c = l.split(); absent[c].add(s)
    mem = {p: set(l.strip() for l in open(os.path.join(B, '..', 'fig4c_pi', 'groups', p + '.txt')) if l.strip()) for p in POPS}
    for nm in POPS + ['%s_%s' % pr for pr in PAIRS]:
        grp = set().union(*[mem[p] for p in POPS if p in nm])
        for c in CHR:
            if nm in POPS:
                r = pima('pima_present', nm, [c])
            else:
                r = fstrows('present', nm, [c])
            st = sorted(int(x['BIN_START']) for x in r)
            gap = max([b - a for a, b in zip(st, st[1:])] + [0]) if st else 'NA'
            nb = (LEN[c] - 1) // 50000
            if nm in POPS:
                vals = ['%.2f' % (sum(float(x['MEAN_CALLED_ALLELES']) for x in r) / len(r)), '%.6g' % mean(r, 'PI_MISSING_AWARE')[0], '%.6g' % mean(r, 'PI_VCFTOOLS')[0], '']
            else:
                vals = ['', '', '', '%.6g' % fst(r)[0]]
            h.write('\t'.join(map(str, [nm, c, LEN[c], sites['present'][c], sites['strict'][c], len(absent[c] & grp), len(r), nb, st[0] if st else 'NA',
                    max(int(x['BIN_END']) for x in r) if r else 'NA', (max(int(x['BIN_END']) for x in r) - LEN[c]) if r else 'NA', gap] + vals)) + '\n')
# ---------- 11 clean chromosomes: callset comparison (old VCF vs final) ----------
with open(os.path.join(B, 'compare_callsets_11clean_chroms.tsv'), 'w') as h:
    h.write('metric\tgroup\told_VCF_vcftools\tfinal_strict_vcftools\told_VCF_missing_aware\tfinal_strict_missing_aware\tmean_called_alleles_old\tmean_called_alleles_final\n')
    for nm in POPS:
        o = pima('pima_oldvcf', nm, CLEAN) if os.path.isdir(os.path.join(B, 'out', 'pima_oldvcf')) else []
        f = pima('pima_strict', nm, CLEAN)
        mc = lambda r: sum(float(x['MEAN_CALLED_ALLELES']) for x in r) / len(r) if r else float('nan')
        h.write('pi\t%s\t%s\t%.6g\t%s\t%.6g\t%.1f\t%.1f\n' % (nm, '%.6g' % mean(o, 'PI_VCFTOOLS')[0] if o else 'NA', mean(f, 'PI_VCFTOOLS')[0], '%.6g' % mean(o, 'PI_MISSING_AWARE')[0] if o else 'NA', mean(f, 'PI_MISSING_AWARE')[0], mc(o), mc(f)))
for fn in ('pi_fst_genomewide.tsv', 'compare_callsets_11clean_chroms.tsv'):
    print(open(os.path.join(B, fn)).read())
