#!/usr/bin/env python3
"""Fig. 5e/5f statistics for each dSV/dSNP definition, using the functions of the
original script 50 (same bins, flank, 20 random background centres per focal dSV,
seed 20260630 + 31) and a two-sided paired Wilcoxon signed-rank test on the per-focal
+-1 Mb totals (observed vs mean background)."""
import csv, importlib.util, sys
import numpy as np
import pandas as pd
from scipy import stats

W = '${CLUSTER_WORK}/enh_B'
TAG = sys.argv[1] if len(sys.argv) > 1 else 'syri'
B = '${ANALYSIS_DIR}/21_MS/06_result/dSVs'
R = B + '/results-8.9'
spec = importlib.util.spec_from_file_location('m50', W + '/50_build_dsv_dsnp_colocalization.py')
m50 = importlib.util.module_from_spec(spec); spec.loader.exec_module(m50)

man = pd.read_csv(R + '/config/Sample_Manifest.tsv', sep='\t')
P38 = sorted(man.loc[man.Include_All38 == 'Yes', 'Sample_ID'])
P35 = sorted(man.loc[man.Include_African35 == 'Yes', 'Sample_ID'])
lens = m50.read_fai(B + '/input/Africa_hap2.fa.fai')
CHR = [f'chr{i:02d}B' for i in range(1, 17)]


def rows_of(df, col):
    out = []
    for r in df.to_dict('records'):
        r = {k: ('' if (isinstance(v, float) and np.isnan(v)) else v) for k, v in r.items()}
        r['Samples'] = r[col]
        out.append(r)
    return out


def run(label, dsv_rows, dsnp_rows, panel):
    pset = set(panel)
    counts = m50.build_sample_chrom_counts(dsv_rows, dsnp_rows, CHR, panel, label)
    corr = m50.pearson_summary(counts, label)
    idx = m50.prepare_dsnp_index(dsnp_rows, CHR, pset)
    foc = m50.prepare_focal_infos(dsv_rows, set(CHR), pset)
    n_bins = (2 * m50.FLANK_BP) // m50.BIN_SIZE
    rng = np.random.default_rng(20260630 + 31)
    pos = {c: [it[0] for it in v] for c, v in idx.items()}
    obs, bg, prof_o, prof_b = [], [], np.zeros(n_bins), np.zeros(n_bins)
    for f in foc:
        items = idx.get(f['Chrom'], []); p = pos.get(f['Chrom'], [])
        o = m50.count_carrier_bins(items, f['Focal_Pos'], f['Focal_Carriers'], n_bins, p)
        b = np.zeros(n_bins)
        for _ in range(20):
            b += m50.count_carrier_bins(items, m50.random_center(lens.get(f['Chrom'], 0), rng), f['Focal_Carriers'], n_bins, p)
        b /= 20
        obs.append(o.sum()); bg.append(b.sum()); prof_o += o; prof_b += b
    obs = np.array(obs); bg = np.array(bg)
    w = stats.wilcoxon(obs, bg) if len(obs) > 10 else None
    # central +-50 kb enrichment for the profile
    mid = n_bins // 2
    res = {'label': label, 'panel_n': len(panel), 'n_dSV': len(dsv_rows), 'n_dSNP': len(dsnp_rows),
           'e_N': corr['N'], 'e_r': float(corr['Pearson_r']), 'e_P': float(corr['P_value']),
           'e_dSV_total': corr['dSV_Total'], 'e_dSNP_total': corr['dSNP_Total'],
           'f_pairs': len(obs), 'f_obs_mean': obs.mean() if len(obs) else np.nan,
           'f_bg_mean': bg.mean() if len(bg) else np.nan,
           'f_obs_median': np.median(obs) if len(obs) else np.nan, 'f_bg_median': np.median(bg) if len(bg) else np.nan,
           'f_frac_obs_gt_bg': (obs > bg).mean() if len(obs) else np.nan,
           'f_wilcoxon_P': w.pvalue if w else np.nan,
           'f_ratio_central100kb': prof_o[mid - 5:mid + 5].sum() / max(prof_b[mid - 5:mid + 5].sum(), 1e-9)}
    pd.DataFrame(counts).to_csv(f'{W}/coloc_{TAG}/{label}.fig_e_counts.tsv', sep='\t', index=False)
    pd.DataFrame({'Focal_DSV_ID': [f['Focal_DSV_ID'] for f in foc], 'Carriers': [','.join(f['Focal_Carriers']) for f in foc],
                  'Observed': obs, 'Background': bg}).to_csv(f'{W}/coloc_{TAG}/{label}.fig_f_pairs.tsv', sep='\t', index=False)
    pd.DataFrame({'bin_mid_kb': (np.arange(n_bins) * m50.BIN_SIZE - m50.FLANK_BP + 5000) / 1000,
                  'Observed': prof_o, 'Background': prof_b}).to_csv(f'{W}/coloc_{TAG}/{label}.fig_f_profile.tsv', sep='\t', index=False)
    print({k: (f'{v:.4g}' if isinstance(v, float) else v) for k, v in res.items()}, flush=True)
    return res


import os
os.makedirs(f'{W}/coloc_{TAG}', exist_ok=True)
out = []
d70 = pd.read_csv(R + '/06_dSNP_phoenix_polarity/polarity/dsnp_v2_phoenix_alt_derived_candidates.tsv', sep='\t',
                  usecols=['SNP_ID', 'Chrom', 'Pos', 'Samples'])
snp70 = rows_of(d70, 'Samples')
if TAG == 'syri':  # validation runs (independent of the oleifera callability source)
    f1480 = pd.read_csv(R + '/01_core_dsv/dsv_hap38_candidates.tsv', sep='\t', low_memory=False)
    out.append(run('V0_All38_formal_1480', rows_of(f1480, 'Samples'), snp70, P38))
    a924 = pd.read_csv(R + '/04_species_qc/african35_reidentified_candidates.tsv', sep='\t', low_memory=False)
    out.append(run('V1_African35_924_dSNP70269', rows_of(a924, 'Samples_African35'), snp70, P35))

snp = pd.read_csv(f'{W}/dsnp_candidates_eg35.{TAG}.tsv', sep='\t', dtype={'phx': str, 'ole': str})
snp['phx'] = snp.phx == 'True'; snp['ole'] = snp.ole == 'True'
for k in (1, 2, 3):
    for sch in ('a_Phoenix', 'b_Oleifera', 'c_Phoenix_and_Oleifera'):
        m = {'a_Phoenix': snp.phx, 'b_Oleifera': snp.ole, 'c_Phoenix_and_Oleifera': snp.phx & snp.ole}[sch]
        s = snp[m & (snp.n35 <= k)]
        dv = pd.read_csv(f'{W}/dsv_{TAG}/dsv_{sch}_k{k}.tsv', sep='\t', low_memory=False)
        out.append(run(f'EG35_{sch}_k{k}', rows_of(dv, 'Samples_African35'), rows_of(s, 'Samples_African35'), P35))
        if k == 1:  # same dSV set against the published dSNP set (restricted to African35 carriers)
            out.append(run(f'EG35_{sch}_k{k}_x_dSNP70269', rows_of(dv, 'Samples_African35'), snp70, P35))
pd.DataFrame(out).to_csv(f'{W}/coloc_{TAG}/coloc_summary.tsv', sep='\t', index=False)
