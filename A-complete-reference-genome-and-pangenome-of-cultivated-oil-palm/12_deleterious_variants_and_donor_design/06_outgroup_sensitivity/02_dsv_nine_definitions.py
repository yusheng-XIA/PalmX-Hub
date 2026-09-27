#!/usr/bin/env python3
"""Re-identify candidate dSVs inside E. guineensis only (African35 panel) under three
polarization schemes and three carrier thresholds.
  a  Phoenix outgroup (ALT_Polarity_v5 == ALT_Derived, unchanged)
  b  E. oleifera outgroup: both MZ4 haplotypes REF-like (non-carrier and aligned/spanned)
  c  a AND b
Carrier threshold k: 1 <= carriers among African35 <= k (k=1 equals freq <= 0.05 of 35).
Functional filter and SV types identical to the formal rule."""
import pickle, re
import numpy as np
import pandas as pd

W = 'work/outgroup_sensitivity'
import sys as _s, os as _o
TAG = _s.argv[1] if len(_s.argv) > 1 else 'syri'
OD = W + '/dsv_' + TAG
_o.makedirs(OD, exist_ok=True)
MZDIR = W + ('/mz4' if TAG == 'syri' else '/mz4_' + TAG)
R = 'dsv_analysis/results_hap38'
man = pd.read_csv(R + '/config/Sample_Manifest.tsv', sep='\t')
AF = set(man.loc[man.Include_African35 == 'Yes', 'Sample_ID'])
assert len(AF) == 35, len(AF)
MZ = ['meizhou4_hap1', 'meizhou4_hap2']
mz = {h: pickle.load(open(f'{MZDIR}/{h}.pkl', 'rb')) for h in MZ}

c = pd.read_csv(R + '/01_core_dsv/sv_catalog.dsv_hap38.tsv', sep='\t', low_memory=False)
print('catalog', len(c))


def split(v):
    v = '' if pd.isna(v) else str(v)
    return [x for x in re.split(r'[;,|]+', v) if x and x != 'NA']


c['car38'] = c.Samples.map(split)
c['car35'] = c.car38.map(lambda s: [x for x in s if x in AF])
c['n35'] = c.car35.map(len)
c['Samples_African35'] = c.car35.map(';'.join)
func = (c.CDS_Overlap_Count.fillna(0).astype(float) > 0) | (c.Conserved_Region_Overlap == 'Yes') | \
       (c.Conserved_Region_Covered_Bases.fillna(0).astype(float) > 0)
c['func'] = func
c['typeok'] = c.SVTYPE.isin(['DEL', 'INS', 'DUP', 'INV'])
c['phx'] = c.ALT_Polarity_v5 == 'ALT_Derived'

# reproduce formal All38 and the existing African35 set
formal = c.typeok & (c.Frequency <= 0.05) & c.phx & ((c.CDS_Overlap_Count > 0) | (c.Conserved_Region_Overlap == 'Yes'))
print('formal All38 recomputed', formal.sum(), 'flag', (c.dSV_v5_Flag == 'Yes').sum())


def cov_frac(blocks, s, e):
    if blocks is None or e <= s:
        return 0.0
    i = np.searchsorted(blocks[:, 1], s, side='right')
    tot = 0
    while i < len(blocks) and blocks[i, 0] < e:
        tot += min(e, blocks[i, 1]) - max(s, blocks[i, 0])
        i += 1
    return tot / (e - s)


def spans(blocks, s, e):
    if blocks is None:
        return False
    i = np.searchsorted(blocks[:, 1], e, side='left')  # first block with end >= e
    return i < len(blocks) and blocks[i, 0] <= s and blocks[i, 1] >= e


def del_frac(dels, s, e):
    return cov_frac(dels, s, e) if dels is not None else 0.0


def mz_state(row, h):
    if h in row.car38:
        return 'ALT_carrier'
    d = mz[h]
    ch = row.Chrom
    b = d['blocks'].get(ch)
    t = row.SVTYPE
    if t == 'INS':
        p = int(row.Pos)
        if not spans(b, p - 100, p + 100):
            return 'Unaligned'
        ins = d['ins'].get(ch)
        if ins is not None and len(ins):
            lo = np.searchsorted(ins[:, 0], p - 50); hi = np.searchsorted(ins[:, 0], p + 51)
            need = max(50, int(round(0.5 * abs(float(row.Abs_SVLEN)))))
            if lo < hi and ins[lo:hi, 1].max() >= need:
                return 'ALT_like_uncalled'
        return 'REF_like'
    s, e = int(row.Interval_Start0), int(row.Interval_End0)
    f = cov_frac(b, s, e)
    if t == 'DEL' and del_frac(d['del'].get(ch), s, e) >= 0.5:
        return 'ALT_like_uncalled'
    if f >= 0.8:
        return 'REF_like'
    if f <= 0.2:
        return 'Unaligned'
    return 'Partial'


# only need MZ4 state for rows that could qualify
cand = c.typeok & c.func & (c.n35 >= 1) & (c.n35 <= 3)
sub = c[cand].copy()
for h in MZ:
    sub['MZ_' + h] = [mz_state(r, h) for r in sub.itertuples()]
sub['ole'] = (sub['MZ_meizhou4_hap1'] == 'REF_like') & (sub['MZ_meizhou4_hap2'] == 'REF_like')
sub['formal1480'] = formal[cand].values
keep_cols = ['SV_Key', 'Chrom', 'Pos', 'Start', 'End', 'Interval_Start0', 'Interval_End0', 'SVTYPE', 'Abs_SVLEN',
             'Sample_Count', 'Samples', 'Samples_African35', 'n35', 'ALT_Polarity_v5', 'ALT_Polarity_v5_Evidence',
             'MZ_meizhou4_hap1', 'MZ_meizhou4_hap2', 'phx', 'ole', 'CDS_Overlap_Count', 'Conserved_Region_Overlap',
             'dSV_v5_Evidence_Class', 'Gene_Impact_Class', 'Gene_IDs', 'CDS_Gene_IDs', 'formal1480']
sub[keep_cols].to_csv(OD + '/dsv_candidates_eg35.tsv', sep='\t', index=False)

a35 = pd.read_csv(R + '/04_species_qc/african35_reidentified_candidates.tsv', sep='\t', usecols=['SV_Key'])
A924 = set(a35.SV_Key)
rows = []
for k in (1, 2, 3):
    for sch, m in (('a_Phoenix', sub.phx), ('b_Oleifera', sub.ole), ('c_Phoenix_and_Oleifera', sub.phx & sub.ole)):
        d = sub[m & (sub.n35 <= k)]
        keys = set(d.SV_Key)
        tc = d.SVTYPE.value_counts().to_dict()
        rows.append({'k_max_carriers35': k, 'scheme': sch, 'n_dSV': len(d),
                     'DEL': tc.get('DEL', 0), 'INS': tc.get('INS', 0), 'DUP': tc.get('DUP', 0), 'INV': tc.get('INV', 0),
                     'singleton_frac': round((d.n35 == 1).mean(), 4) if len(d) else np.nan,
                     'overlap_A924': len(keys & A924), 'overlap_formal1480': int(d.formal1480.sum()),
                     'carrier_haps': d.Samples_African35.str.split(';').explode().nunique()})
        d.to_csv(f'{OD}/dsv_{sch}_k{k}.tsv', sep='\t', index=False)
summ = pd.DataFrame(rows)
summ.to_csv(OD + '/dsv_scheme_summary.tsv', sep='\t', index=False)
print(summ.to_string())

# per-haplotype distribution
per = {}
for k in (1, 2, 3):
    for sch, m in (('a_Phoenix', sub.phx), ('b_Oleifera', sub.ole), ('c_Phoenix_and_Oleifera', sub.phx & sub.ole)):
        d = sub[m & (sub.n35 <= k)]
        per[f'{sch}_k{k}'] = d.Samples_African35.str.split(';').explode().value_counts()
per = pd.DataFrame(per).reindex(sorted(AF)).fillna(0).astype(int)
formal_ph = c[formal].car38.explode().value_counts()
per.insert(0, 'formal1480_All38', formal_ph.reindex(per.index).fillna(0).astype(int))
per.to_csv(OD + '/dsv_per_haplotype.tsv', sep='\t')
print(per.sort_values('b_Oleifera_k1', ascending=False).to_string())

# MZ4 state among the 924 (scheme a k=1): why oleifera polarization drops them
a1 = sub[sub.phx & (sub.n35 == 1)]
st = a1.groupby(['MZ_meizhou4_hap1', 'MZ_meizhou4_hap2']).size().rename('n').reset_index()
st.to_csv(OD + '/dsv_a1_mz4_state.tsv', sep='\t', index=False)
print(st.to_string())
print('a1 = 924 check', len(a1), 'in A924', len(set(a1.SV_Key) & A924))
print('A924 not in a1', len(A924 - set(a1.SV_Key)))
# of the 631 formal dSVs carried by non-African haplotypes, how many are shared with an African carrier
f = c[formal]
print('formal with any African carrier', (f.n35 > 0).sum(), 'formal only non-African', (f.n35 == 0).sum())
