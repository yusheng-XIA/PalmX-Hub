#!/usr/bin/env python3
"""Confounding controls for Fig. 5e,f (dSV-dSNP association).

Sets (dSV set / derived-SNV rule / panel):
  A608  : formal 1,480 dSVs (All38)             / n38 == 1, Phoenix ALT-derived      / 38 haplotypes x 16 chr = 608
  S560  : strict African35 dSVs (c, 1/35; 550)  / n35 == 1, Phoenix + both MZ4 REF   / 35 x 16 = 560
  V560  : 924 African35 dSVs x 70,269 dSNPs     / (as used in the previous validation) / 560
Derived-SNV classes (same frequency and polarity rule within a set):
  dSNP (formal functional classes), SYN (synonymous, outside conserved proxy; disjoint from dSNP),
  SYNall (all synonymous, including those inside conserved intervals, which are also dSNPs),
  NCNC (intron/intergenic/UTR/+-2 kb, outside conserved proxy); NCunpol = NCNC without the polarity
  (and E. oleifera) requirement, same frequency rule (mostly sites without date palm alignment).
  For SYN/SYNall (far fewer sites than dSNP) the 5e r is also computed for dSNP randomly thinned to the same
  number of sites (mean of 200 draws), so that the comparison is not driven by count sparsity.
(a) Fig. 5e: raw Pearson r; haplotype-demeaned; two-way (haplotype + chromosome) fixed-effect residual r;
    density per Mb and per functional Mb (CDS U conserved proxy); FE partial r additionally controlling
    the SYN and NCNC counts; Meng-Rosenthal-Rubin test of r(dSV,dSNP) vs r(dSV,neutral) on FE residuals.
(b) Fig. 5f: script-50 counting (carrier-matched derived SNVs within +-1 Mb, 10-kb bins, 20 background
    centres per focal dSV) with three backgrounds:
      B0 uniform on the same chromosome (original; same seed and draw order -> reproduces Fig. 5f)
      B1 same chromosome, +-1 Mb functional bp (CDS U conserved) in the same within-chromosome decile
      B2 as B1 and the background centre itself inside CDS U conserved (as required of a dSV)
    The same background centres are used for dSNP, SYN and NCNC.
"""
import csv, gzip, importlib.util, sys, os
from collections import defaultdict
import numpy as np, pandas as pd
from scipy import stats

W0 = '${CLUSTER_WORK}/enh_B'
W = '${CLUSTER_WORK}/enh_dsv_confound'
B = '${ANALYSIS_DIR}/21_MS/06_result/dSVs'
R = B + '/results-8.9'
CONS = B + '/results/archive_old_versions/05_dsv_v3_multi_outgroup/multi_outgroup_conserved_regions/multi_outgroup_conserved_regions.bed'
OUT = W + '/res'; os.makedirs(OUT, exist_ok=True)
spec = importlib.util.spec_from_file_location('m50', W0 + '/50_build_dsv_dsnp_colocalization.py')
m50 = importlib.util.module_from_spec(spec); spec.loader.exec_module(m50)
CHR = [f'chr{i:02d}B' for i in range(1, 17)]
lens = m50.read_fai(B + '/input/Africa_hap2.fa.fai')
man = pd.read_csv(R + '/config/Sample_Manifest.tsv', sep='\t')
P38 = sorted(man.loc[man.Include_All38 == 'Yes', 'Sample_ID'])
P35 = sorted(man.loc[man.Include_African35 == 'Yes', 'Sample_ID'])
NB = (2 * m50.FLANK_BP) // m50.BIN_SIZE
NC = {'Intron', 'Intergenic', 'UTR', 'Upstream_2kb', 'Downstream_2kb'}
FUNC = {'Missense', 'Stop_Gained', 'Stop_Lost', 'Start_Lost', 'Splice_Region', 'Conserved_Region_Proxy'}

# ---------------------------------------------------------------- functional mask (CDS U conserved)
iv = defaultdict(list)
with open(B + '/input/Africa_hap2.EVM.gff3') as fh:
    for line in fh:
        if line.startswith('#'):
            continue
        f = line.split('\t')
        if len(f) > 4 and f[2] == 'CDS' and f[0] in CHR:
            iv[f[0]].append((int(f[3]) - 1, int(f[4])))
cds_bp = {}
for c in CHR:
    a = sorted(iv[c]); m = []
    for s, e in a:
        if m and s <= m[-1][1]:
            m[-1][1] = max(m[-1][1], e)
        else:
            m.append([s, e])
    cds_bp[c] = sum(e - s for s, e in m)
with open(CONS) as fh:
    for line in fh:
        f = line.split('\t')
        if f[0] in CHR:
            iv[f[0]].append((int(f[1]), int(f[2])))
MASK = {}
for c in CHR:
    a = sorted(iv[c]); m = []
    for s, e in a:
        if m and s <= m[-1][1]:
            m[-1][1] = max(m[-1][1], e)
        else:
            m.append([s, e])
    m = np.array(m, dtype=np.int64)
    cum = np.concatenate([[0], np.cumsum(m[:, 1] - m[:, 0])])
    MASK[c] = (m, cum)
FBP = {c: int(MASK[c][1][-1]) for c in CHR}
print('functional Mb per chr', {c: round(FBP[c] / 1e6, 2) for c in CHR}, 'total', sum(FBP.values()) / 1e6, flush=True)


def fcum(c, x):
    """functional bp in [0, x) (0-based half-open)"""
    m, cum = MASK[c]
    x = np.asarray(x, dtype=np.int64)
    i = np.searchsorted(m[:, 0], x, side='right') - 1
    ii = np.maximum(i, 0)
    val = cum[ii] + np.minimum(np.maximum(x - m[ii, 0], 0), m[ii, 1] - m[ii, 0])
    return np.where(i >= 0, val, 0)


def wfunc(c, centre):
    """functional bp in the +-1 Mb window used by script 50 (positions centre-1Mb .. centre+1Mb, 1-based)."""
    lo = np.maximum(np.asarray(centre) - m50.FLANK_BP - 1, 0); hi = np.minimum(np.asarray(centre) + m50.FLANK_BP, lens[c])
    return fcum(c, hi) - fcum(c, lo)


def infunc(c, p1):
    m, _ = MASK[c]; p0 = int(p1) - 1
    i = np.searchsorted(m[:, 0], p0, side='right') - 1
    return i >= 0 and m[i, 0] <= p0 < m[i, 1]


GRID = {}; DEC = {}
for c in CHR:
    g = np.arange(5000, lens[c], 10000)
    wf = wfunc(c, g)
    GRID[c] = (g, wf)
    DEC[c] = np.quantile(wf, np.linspace(0, 1, 11)[1:-1])


def decile(c, x):
    return np.searchsorted(DEC[c], x, side='right')


# ---------------------------------------------------------------- inputs
def rows_of(df, col):
    out = []
    for r in df.to_dict('records'):
        r = {k: ('' if (isinstance(v, float) and np.isnan(v)) else v) for k, v in r.items()}
        r['Samples'] = r[col]; out.append(r)
    return out


snv = pd.read_csv(W + '/snv_classes.tsv.gz', sep='\t', dtype={'Samples_African35': str}, low_memory=False)
snv['phx'] = snv.phx.astype(str) == 'True'; snv['ole'] = snv.ole.astype(str) == 'True'
snv['Samples_African35'] = snv.Samples_African35.fillna('')
snv['cls'] = np.where(snv.FC.isin(FUNC), 'dSNP', np.where(snv.FC == 'Synonymous', 'SYN', np.where(snv.FC.isin(NC), 'NCNC', 'other')))
a_rule = (snv.n38 == 1) & snv.phx
s_rule = (snv.n35 == 1) & snv.phx & snv.ole
d70 = pd.read_csv(R + '/06_dSNP_phoenix_polarity/polarity/dsnp_v2_phoenix_alt_derived_candidates.tsv', sep='\t',
                  usecols=['SNP_ID', 'Chrom', 'Pos', 'Samples'])
sn35 = pd.read_csv(W0 + '/dsnp_candidates_eg35.syri.tsv', sep='\t', dtype={'phx': str, 'ole': str})
sn35 = sn35[(sn35.phx == 'True') & (sn35.ole == 'True') & (sn35.n35 == 1)]
# consistency: the class table reproduces both published dSNP sets
k70 = set(zip(d70.Chrom, d70.Pos)); kA = set(zip(snv.Chrom[a_rule & (snv.cls == 'dSNP')], snv.Pos[a_rule & (snv.cls == 'dSNP')]))
k35 = set(zip(sn35.Chrom, sn35.Pos)); kS = set(zip(snv.Chrom[s_rule & (snv.cls == 'dSNP')], snv.Pos[s_rule & (snv.cls == 'dSNP')]))
print('dSNP reproduce All38', len(k70), len(kA), len(k70 ^ kA), '| strict', len(k35), len(kS), len(k35 ^ kS), flush=True)
assert len(k70 ^ kA) == 0 and len(k35 ^ kS) == 0

SETS = {
    'A608': dict(dsv=rows_of(pd.read_csv(R + '/01_core_dsv/dsv_hap38_candidates.tsv', sep='\t', low_memory=False), 'Samples'),
                 panel=P38, rule=a_rule, unpol=(snv.n38 == 1), col='Samples', dsnp=rows_of(d70, 'Samples')),
    'S560': dict(dsv=rows_of(pd.read_csv(W0 + '/dsv_syri/dsv_c_Phoenix_and_Oleifera_k1.tsv', sep='\t', low_memory=False), 'Samples_African35'),
                 panel=P35, rule=s_rule, unpol=(snv.n35 == 1), col='Samples_African35', dsnp=rows_of(sn35, 'Samples_African35')),
    'V560': dict(dsv=rows_of(pd.read_csv(R + '/04_species_qc/african35_reidentified_candidates.tsv', sep='\t', low_memory=False), 'Samples_African35'),
                 panel=P35, rule=a_rule, unpol=(snv.n38 == 1), col='Samples', dsnp=rows_of(d70, 'Samples')),
}


def carrier_index(rows, pset):
    d = defaultdict(list)
    for r in rows:
        c = r.get('Chrom')
        if c not in CHR:
            continue
        for smp in m50.parse_samples(r.get('Samples')):
            if smp in pset:
                d[(c, smp)].append(int(r['Pos']))
    return {k: np.array(sorted(v), dtype=np.int64) for k, v in d.items()}


EMPTY = np.zeros(0, dtype=np.int64)


def fast_bins(cidx, chrom, center, carriers):
    arr = np.zeros(NB)
    for smp in carriers:
        P = cidx.get((chrom, smp), EMPTY)
        lo = np.searchsorted(P, center - m50.FLANK_BP, side='left'); hi = np.searchsorted(P, center + m50.FLANK_BP, side='right')
        b = np.minimum((P[lo:hi] - center + m50.FLANK_BP) // m50.BIN_SIZE, NB - 1)
        arr += np.bincount(b, minlength=NB)
    return arr / len(carriers)


def fe_resid(df, v):
    x = df.pivot(index='Sample', columns='Chrom', values=v).values.astype(float)
    r = x - x.mean(1, keepdims=True) - x.mean(0, keepdims=True) + x.mean()
    return r.ravel()


def pr(x, y, df_res):
    r = float(np.corrcoef(x, y)[0, 1]); t = r * np.sqrt(df_res / max(1e-15, 1 - r * r))
    return r, float(2 * stats.t.sf(abs(t), df_res))


def resid_on(y, Z):
    beta, *_ = np.linalg.lstsq(Z, y, rcond=None); return y - Z @ beta


def mrr(r1, r2, r12, n):
    """Meng, Rosenthal & Rubin (1992) Z for two dependent correlations sharing one variable."""
    z1, z2 = np.arctanh(r1), np.arctanh(r2); rm2 = (r1 * r1 + r2 * r2) / 2
    f = min(1.0, (1 - r12) / (2 * (1 - rm2))); h = (1 - f * rm2) / (1 - rm2)
    z = (z1 - z2) * np.sqrt((n - 3) / (2 * (1 - r12) * h)); return float(z), float(2 * stats.norm.sf(abs(z)))


e_rows, f_rows = [], []
for S, cfg in SETS.items():
    panel = cfg['panel']; pset = set(panel)
    sub = snv[cfg['rule']]
    ys = {'dSNP': cfg['dsnp'],
          'SYN': rows_of(sub[sub.cls == 'SYN'][['Chrom', 'Pos', cfg['col']]], cfg['col']),
          'SYNall': rows_of(sub[sub.FE == 'Synonymous'][['Chrom', 'Pos', cfg['col']]], cfg['col']),
          'NCNC': rows_of(sub[sub.cls == 'NCNC'][['Chrom', 'Pos', cfg['col']]], cfg['col'])}
    un = snv[cfg['unpol'] & (snv.cls == 'NCNC')]
    ys['NCunpol'] = rows_of(un[['Chrom', 'Pos', cfg['col']]], cfg['col'])
    YS = ('dSNP', 'SYN', 'SYNall', 'NCNC', 'NCunpol')
    print(S, {k: len(v) for k, v in ys.items()}, flush=True)
    # ------------------------------------------------ (a) Fig. 5e
    base = pd.DataFrame(m50.build_sample_chrom_counts(cfg['dsv'], ys['dSNP'], CHR, panel, S)).rename(columns={'dSV_Count': 'dSV', 'dSNP_Count': 'dSNP'})
    for y in YS[1:]:
        t = pd.DataFrame(m50.build_sample_chrom_counts(cfg['dsv'], ys[y], CHR, panel, S))
        assert (t.Sample.values == base.Sample.values).all() and (t.Chrom.values == base.Chrom.values).all()
        base[y] = t.dSNP_Count.values
    base['Chrom_Mb'] = base.Chrom.map(lambda c: lens[c] / 1e6); base['Functional_Mb'] = base.Chrom.map(lambda c: FBP[c] / 1e6)
    base['CDS_Mb'] = base.Chrom.map(lambda c: cds_bp[c] / 1e6)
    base.drop(columns=['Analysis']).to_csv(f'{OUT}/e_counts_{S}.tsv', sep='\t', index=False)
    n = len(base); H = base.Sample.nunique(); C = base.Chrom.nunique(); dfe = n - H - C
    X = base.dSV.values.astype(float)
    Xfe = fe_resid(base, 'dSV')
    fe = {v: fe_resid(base, v) for v in YS}
    for y in YS:
        Y = base[y].values.astype(float)
        r, p = stats.pearsonr(X, Y); e_rows.append((S, n, y, '1_raw_count', r, p, n - 2))
        xc = X - base.groupby('Sample').dSV.transform('mean').values; yc = Y - base.groupby('Sample')[y].transform('mean').values
        r, p = pr(xc, yc, n - H - 1); e_rows.append((S, n, y, '2_haplotype_demeaned', r, p, n - H - 1))
        r, p = pr(Xfe, fe[y], dfe); e_rows.append((S, n, y, '3_two_way_FE', r, p, dfe))
        for dn, col in (('4_density_per_Mb', 'Chrom_Mb'), ('5_density_per_functional_Mb', 'Functional_Mb'), ('5b_density_per_CDS_Mb', 'CDS_Mb')):
            r, p = stats.pearsonr(X / base[col].values, Y / base[col].values); e_rows.append((S, n, y, dn, r, p, n - 2))
            tmp = base.assign(xd=X / base[col].values, yd=Y / base[col].values)
            xd = tmp.xd - tmp.groupby('Sample').xd.transform('mean'); yd = tmp.yd - tmp.groupby('Sample').yd.transform('mean')
            r, p = pr(xd, yd, n - H - 1); e_rows.append((S, n, y, dn + '_haplotype_demeaned', r, p, n - H - 1))
        r, p = stats.spearmanr(Xfe, fe[y]); e_rows.append((S, n, y, '3s_two_way_FE_spearman', r, p, dfe))
    Z = np.column_stack([fe['SYN'], fe['NCNC']])
    r, p = pr(resid_on(Xfe, Z), resid_on(fe['dSNP'], Z), dfe - 2)
    e_rows.append((S, n, 'dSNP', '6_two_way_FE_plus_neutral_counts', r, p, dfe - 2))
    # thinned dSNP (same number of sites as SYN / SYNall), mean r over 200 draws
    trng = np.random.default_rng(20260926)
    for y in ('SYN', 'SYNall'):
        k = int(base[y].sum()); rr_raw, rr_fe = [], []
        for _ in range(200):
            keep = trng.choice(len(ys['dSNP']), size=k, replace=False)
            t = pd.DataFrame(m50.build_sample_chrom_counts(cfg['dsv'], [ys['dSNP'][j] for j in keep], CHR, panel, S))
            yy = t.dSNP_Count.values.astype(float)
            rr_raw.append(np.corrcoef(X, yy)[0, 1])
            tt = base[['Sample', 'Chrom']].assign(v=yy)
            rr_fe.append(np.corrcoef(Xfe, fe_resid(tt, 'v'))[0, 1])
        e_rows.append((S, n, f'dSNP_thinned_to_{y}', '1_raw_count', float(np.mean(rr_raw)), np.nan, k))
        e_rows.append((S, n, f'dSNP_thinned_to_{y}', '3_two_way_FE', float(np.mean(rr_fe)), np.nan, k))
        e_rows.append((S, n, f'dSNP_thinned_to_{y}', '3_two_way_FE_q025', float(np.quantile(rr_fe, 0.025)), np.nan, k))
        e_rows.append((S, n, f'dSNP_thinned_to_{y}', '3_two_way_FE_q975', float(np.quantile(rr_fe, 0.975)), np.nan, k))
    for y in ('SYN', 'SYNall', 'NCNC', 'NCunpol'):
        r1 = np.corrcoef(Xfe, fe['dSNP'])[0, 1]; r2 = np.corrcoef(Xfe, fe[y])[0, 1]; r12 = np.corrcoef(fe['dSNP'], fe[y])[0, 1]
        z, p = mrr(r1, r2, r12, dfe + 2)
        e_rows.append((S, n, f'dSNP_vs_{y}', '7_MRR_test_FE_r_difference', r1 - r2, p, dfe))
        r1 = np.corrcoef(X, base.dSNP)[0, 1]; r2 = np.corrcoef(X, base[y])[0, 1]; r12 = np.corrcoef(base.dSNP, base[y])[0, 1]
        z, p = mrr(r1, r2, r12, n)
        e_rows.append((S, n, f'dSNP_vs_{y}', '7_MRR_test_raw_r_difference', r1 - r2, p, n - 2))
    print(S, 'e done', flush=True)
    # ------------------------------------------------ (b) Fig. 5f
    foc = m50.prepare_focal_infos(cfg['dsv'], set(CHR), pset)
    rng0 = np.random.default_rng(20260630 + 31)
    rng1 = np.random.default_rng(20260924)
    rng2 = np.random.default_rng(20260925)
    cent = {'B0': [], 'B1': [], 'B2': []}; finfo = []
    for f in foc:
        c = f['Chrom']; L = lens[c]
        cent['B0'].append([m50.random_center(L, rng0) for _ in range(20)])
        wf = int(wfunc(c, f['Focal_Pos'])); d = int(decile(c, wf))
        got = []
        while len(got) < 20:
            x = rng1.integers(1, L + 1, size=200); ok = x[decile(c, wfunc(c, x)) == d]
            got.extend(ok[:20 - len(got)].tolist())
        cent['B1'].append(got)
        m, cum = MASK[c]; tot = cum[-1]; got = []; tries = 0
        while len(got) < 20 and tries < 200:
            u = rng2.integers(0, tot, size=200); i = np.searchsorted(cum, u, side='right') - 1
            x = m[i, 0] + (u - cum[i]) + 1
            ok = x[decile(c, wfunc(c, x)) == d]; got.extend(ok[:20 - len(got)].tolist()); tries += 1
        assert len(got) == 20, (c, d)
        cent['B2'].append(got)
        finfo.append((f['Focal_DSV_ID'], c, f['Focal_Pos'], ','.join(f['Focal_Carriers']), wf, d, infunc(c, f['Focal_Pos'])))
    pairs = pd.DataFrame(finfo, columns=['Focal_DSV_ID', 'Chrom', 'Focal_Pos', 'Carriers', 'Window_functional_bp', 'Decile', 'Focal_in_functional'])
    bidx = np.random.default_rng(20260630).integers(0, len(foc), size=(1000, len(foc)))
    TOT = {}
    for y in YS:
        cidx = carrier_index(ys[y], pset)
        O = np.zeros((len(foc), NB)); G = {k: np.zeros((len(foc), NB)) for k in cent}
        for i, f in enumerate(foc):
            O[i] = fast_bins(cidx, f['Chrom'], f['Focal_Pos'], f['Focal_Carriers'])
            for k in cent:
                for x in cent[k][i]:
                    G[k][i] += fast_bins(cidx, f['Chrom'], int(x), f['Focal_Carriers'])
                G[k][i] /= 20
        if y == 'dSNP':  # the fast counter equals script 50's count_carrier_bins (checked on the first 50 focal dSVs)
            idx = m50.prepare_dsnp_index(ys[y], CHR, pset); pos = {c: [it[0] for it in v] for c, v in idx.items()}
            for i, f in enumerate(foc[:50]):
                a = m50.count_carrier_bins(idx.get(f['Chrom'], []), f['Focal_Pos'], f['Focal_Carriers'], NB, pos.get(f['Chrom'], []))
                assert np.allclose(a, O[i]), ('fast counter mismatch', i)
                a = m50.count_carrier_bins(idx.get(f['Chrom'], []), cent['B0'][i][0], f['Focal_Carriers'], NB, pos.get(f['Chrom'], []))
                assert np.allclose(a, fast_bins(cidx, f['Chrom'], cent['B0'][i][0], f['Focal_Carriers']))
        o = O.sum(1); pairs[f'Obs_{y}'] = o
        prof = pd.DataFrame({'bin_mid_kb': (np.arange(NB) * m50.BIN_SIZE - m50.FLANK_BP + 5000) / 1000, 'Observed': O.sum(0)})
        for k in cent:
            g = G[k].sum(1); pairs[f'{k}_{y}'] = g; prof[k] = G[k].sum(0)
            w = stats.wilcoxon(o, g)
            bt = o[bidx].mean(1) / np.maximum(g[bidx].mean(1), 1e-12)
            TOT[(y, k)] = (o, g)
            mid = NB // 2
            f_rows.append((S, len(foc), y, k, o.mean(), g.mean(), o.mean() / g.mean(), *np.quantile(bt, [0.025, 0.975]),
                           (o > g).mean(), w.pvalue, O[:, mid - 5:mid + 5].sum() / max(G[k][:, mid - 5:mid + 5].sum(), 1e-9)))
        prof.to_csv(f'{OUT}/f_profile_{S}_{y}.tsv', sep='\t', index=False)
    for k in cent:  # enrichment of dSNP relative to each neutral class (ratio of ratios, paired bootstrap over focal dSVs)
        od, gd = TOT[('dSNP', k)]
        for y in ('SYN', 'SYNall', 'NCNC', 'NCunpol'):
            on, gn = TOT[(y, k)]
            rr = (od.mean() / gd.mean()) / (on.mean() / gn.mean())
            bt = (od[bidx].mean(1) / gd[bidx].mean(1)) / (on[bidx].mean(1) / gn[bidx].mean(1))
            pb = 2 * min((bt <= 1).mean(), (bt >= 1).mean()); pb = max(pb, 1 / 1000)
            f_rows.append((S, len(foc), f'dSNP_over_{y}', k, np.nan, np.nan, rr, *np.quantile(bt, [0.025, 0.975]), np.nan, pb, np.nan))
    pairs.to_csv(f'{OUT}/f_pairs_{S}.tsv', sep='\t', index=False)
    print(S, 'f done', flush=True)

E = pd.DataFrame(e_rows, columns=['Set', 'n', 'Y', 'Metric', 'r', 'P', 'df'])
E.to_csv(f'{OUT}/e_stats.tsv', sep='\t', index=False)
F = pd.DataFrame(f_rows, columns=['Set', 'Focal_dSVs', 'Y', 'Background', 'Obs_mean', 'Bg_mean', 'Ratio', 'Ratio_lo95', 'Ratio_hi95',
                                  'Frac_obs_gt_bg', 'P', 'Ratio_central100kb'])
F.to_csv(f'{OUT}/f_stats.tsv', sep='\t', index=False)
pd.set_option('display.width', 250); pd.set_option('display.max_rows', 500)
print(E.to_string()); print(F.to_string())
print('DONE')
