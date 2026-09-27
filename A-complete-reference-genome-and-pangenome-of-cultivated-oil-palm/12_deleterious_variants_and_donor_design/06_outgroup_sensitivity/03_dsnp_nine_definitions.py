#!/usr/bin/env python3
"""Rebuild dSNPs inside E. guineensis (African35) with the same rule as the formal set
(PASS, Mean_Qual >= 30, not within 50 bp of an SV breakpoint, functional class in
Missense/Stop/Start/Splice/Conserved-proxy) but with carriers counted among the 35
E. guineensis haplotypes (1 <= n35 <= 3 kept for sensitivity), then polarize with
Phoenix (script 41 logic, unchanged) and with the two E. oleifera haplotypes (SyRI)."""
import csv, importlib.util, pickle, sys
from collections import Counter, defaultdict
import numpy as np

W = 'work/outgroup_sensitivity'
import sys as _s, os as _o
import os as _os
SCRIPTS = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))  # 12_deleterious_variants_and_donor_design/
TAG = _s.argv[1] if len(_s.argv) > 1 else 'syri'
MZDIR = W + ('/mz4' if TAG == 'syri' else '/mz4_' + TAG)
B = 'dsv_analysis'
R = B + '/results_hap38'


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m); return m


m40 = load('m40', SCRIPTS + '/03_build_dsnp_catalog.py')
m41 = load('m41', SCRIPTS + '/04_polarize_dsnp_date_palm.py')

AF = set()
with open(R + '/config/Sample_Manifest.tsv') as fh:
    for r in csv.DictReader(fh, delimiter='\t'):
        if r['Include_African35'] == 'Yes':
            AF.add(r['Sample_ID'])
assert len(AF) == 35

# 1-2. load annotation, then stream the population catalogue and annotate on the fly
lengths = m40.read_fai(B + '/input/Africa_hap2.fa.fai')
seqs = m40.read_fasta(B + '/input/Africa_hap2.fa', set(lengths))
ann = m40.Annotation(B + '/input/Africa_hap2.EVM.gff3', seqs, 2000, 2, 100000)
print('annotation loaded', flush=True)
func = []
stat = Counter()
KEEP = ['SNP_ID', 'Chrom', 'Pos', 'Ref', 'Alt', 'Sample_Count', 'Samples']
with open(R + '/05_dSNP_minimap_hap38/snp_population_catalog.tsv') as fh:
    rd = csv.DictReader(fh, delimiter='\t')
    for r in rd:
        stat['all'] += 1
        s = r['Samples'].split(';')
        n38 = len(s)
        if n38 > 6:
            continue
        s35 = [x for x in s if x in AF]
        n35 = len(s35)
        if not (1 <= n35 <= 3 or n38 == 1):
            continue
        if not (r['Filter'] == 'PASS' and float(r['Mean_Qual']) >= 30 and r['Near_SV_Breakpoint'] == 'No'):
            continue
        stat['kept_basic'] += 1
        a = ann.annotate(r['Chrom'], int(r['Pos']), r['Ref'], r['Alt'], r['Conserved_Region_Proxy'] == 'Yes')
        if a['Functional_Class'] not in m40.FUNCTIONAL_EFFECTS:
            continue
        o = {k: r[k] for k in KEEP}
        o.update({'n35': n35, 'Samples_African35': ';'.join(s35), 'n38': n38, 'Functional_Class': a['Functional_Class']})
        func.append(o)
del seqs, ann
print('catalog stream', dict(stat), 'functional', len(func), flush=True)
print('check: n38==1 functional (expect 125899):', sum(1 for r in func if r['n38'] == 1), flush=True)

# 3. Phoenix polarization (script 41 functions, MAPQ >= 20)
positions_by_chrom = defaultdict(list); index_by_site = defaultdict(list)
for i, r in enumerate(func):
    positions_by_chrom[r['Chrom']].append(int(r['Pos'])); index_by_site[(r['Chrom'], int(r['Pos']))].append(i)
for ch in positions_by_chrom:
    positions_by_chrom[ch] = sorted(set(positions_by_chrom[ch]))
calls, pst = m41.parse_paf_calls(R + '/06_dSNP_phoenix_polarity/alignments/Phoenix_vs_Africa_hap2.asm20.cs.primary.paf',
                                 func, positions_by_chrom, index_by_site, 20)
for i, r in enumerate(func):
    r.update(m41.summarize_call(r, calls.get(i, [])))
print('phoenix derived n38==1 (expect 70269):',
      sum(1 for r in func if r['n38'] == 1 and r['ALT_Polarity_Phoenix'] == 'ALT_Derived'), flush=True)

# 4. E. oleifera state from SyRI
MZ = ['meizhou4_hap1', 'meizhou4_hap2']
mz = {h: pickle.load(open(f'{MZDIR}/{h}.pkl', 'rb')) for h in MZ}


def inblk(b, p0):
    if b is None or not len(b):
        return False
    i = np.searchsorted(b[:, 1], p0, side='right')
    return i < len(b) and b[i, 0] <= p0


for r in func:
    s38 = set(r['Samples'].split(';'))
    p = int(r['Pos']); ch = r['Chrom']
    for h in MZ:
        if h in s38:
            st = 'ALT_carrier'
        else:
            d = mz[h]
            if not inblk(d['blocks'].get(ch), p - 1):
                st = 'Unaligned'
            elif inblk(d['del'].get(ch), p - 1):
                st = 'Deleted'
            else:
                b = d['snp'].get(ch, {}).get(p)
                if b is None:
                    st = 'REF_like'
                elif b == r['Alt'].upper():
                    st = 'ALT_like_uncalled'
                else:
                    st = 'Third_base'
        r['MZ_' + h] = st
    r['phx'] = r['ALT_Polarity_Phoenix'] == 'ALT_Derived'
    r['ole'] = r['MZ_meizhou4_hap1'] == 'REF_like' and r['MZ_meizhou4_hap2'] == 'REF_like'

cols = ['SNP_ID', 'Chrom', 'Pos', 'Ref', 'Alt', 'Sample_Count', 'Samples', 'n38', 'n35', 'Samples_African35',
        'Functional_Class', 'Phoenix_Base', 'ALT_Polarity_Phoenix', 'MZ_meizhou4_hap1', 'MZ_meizhou4_hap2', 'phx', 'ole']
with open(W + f'/dsnp_candidates_eg35.{TAG}.tsv', 'w', newline='') as out:
    w = csv.DictWriter(out, fieldnames=cols, delimiter='\t', extrasaction='ignore', lineterminator='\n')
    w.writeheader(); w.writerows(func)
summ = Counter()
for r in func:
    for k in (1, 2, 3):
        if 1 <= r['n35'] <= k:
            summ[(k, 'a_Phoenix')] += r['phx']
            summ[(k, 'b_Oleifera')] += r['ole']
            summ[(k, 'c_Phoenix_and_Oleifera')] += r['phx'] and r['ole']
    if r['n38'] == 1 and r['phx']:
        summ[('All38_formal', 'Phoenix')] += 1
        if r['n35'] == 1:
            summ[('All38_formal_African_carrier', 'Phoenix')] += 1
with open(W + f'/dsnp_scheme_summary.{TAG}.tsv', 'w') as out:
    out.write('k_max_carriers35\tscheme\tn_dSNP\n')
    for (k, s), v in sorted(summ.items(), key=lambda x: str(x[0])):
        out.write(f'{k}\t{s}\t{v}\n'); print(k, s, v)
st = Counter((r['MZ_meizhou4_hap1'], r['MZ_meizhou4_hap2']) for r in func if r['n35'] == 1 and r['phx'])
print('MZ4 states among Phoenix-derived n35==1 dSNPs:', st.most_common())
