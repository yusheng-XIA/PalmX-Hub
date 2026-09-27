#!/usr/bin/env python3
"""Neutral-control SNVs for the Fig. 5e,f confounding checks.

Same stream, filters, annotation, Phoenix polarisation and E. oleifera (MZ4, SyRI) states as
enh_B/enhB_03_dsnp.py (itself the formal dSNP rule of scripts 40/41), but every functional
class is kept, not only the dSNP classes.  Classes used downstream:
  dSNP   : Missense/Stop/Start/Splice/Conserved_Region_Proxy  (formal dSNP classes)
  SYN    : Synonymous and outside the conserved-region proxy
  NCNC   : Intron/Intergenic/UTR/Upstream_2kb/Downstream_2kb and outside the conserved-region proxy
Rows kept: PASS, Mean_Qual >= 30, not within 50 bp of an SV breakpoint, n38 <= 6 and
(n38 == 1 or 1 <= n35 <= 3), exactly as enhB_03.
Output: snv_classes.tsv.gz  (Chrom Pos Ref Alt n38 n35 Samples Samples_African35 FC FE phx ole; FE = effect before the
        conserved-proxy relabelling, so FE == Synonymous also covers synonymous sites inside conserved intervals)
"""
import csv, gzip, importlib.util, pickle
from collections import Counter, defaultdict
import numpy as np

W0 = '${CLUSTER_WORK}/enh_B'
W = '${CLUSTER_WORK}/enh_dsv_confound'
B = '${ANALYSIS_DIR}/21_MS/06_result/dSVs'
R = B + '/results-8.9'


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m); return m


m40 = load('m40', W0 + '/40_build_dsnp_hap38_catalog.py')
m41 = load('m41', W0 + '/41_polarize_dsnp_hap38_phoenix.py')
AF = set()
with open(R + '/config/Sample_Manifest.tsv') as fh:
    for r in csv.DictReader(fh, delimiter='\t'):
        if r['Include_African35'] == 'Yes':
            AF.add(r['Sample_ID'])
assert len(AF) == 35

lengths = m40.read_fai(B + '/input/Africa_hap2.fa.fai')
seqs = m40.read_fasta(B + '/input/Africa_hap2.fa', set(lengths))
ann = m40.Annotation(B + '/input/Africa_hap2.EVM.gff3', seqs, 2000, 2, 100000)
print('annotation loaded', flush=True)
rows = []
stat = Counter()
with open(R + '/05_dSNP_minimap_hap38/snp_population_catalog.tsv') as fh:
    for r in csv.DictReader(fh, delimiter='\t'):
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
        rows.append({'Chrom': r['Chrom'], 'Pos': r['Pos'], 'Ref': r['Ref'], 'Alt': r['Alt'], 'Samples': r['Samples'],
                     'n38': n38, 'n35': n35, 'Samples_African35': ';'.join(s35), 'FC': a['Functional_Class'], 'FE': a['Functional_Effect']})
del seqs, ann
print('stream', dict(stat), flush=True)
fc = Counter(r['FC'] for r in rows if r['n38'] == 1)
print('n38==1 classes', fc.most_common(), flush=True)
print('check n38==1 functional (expect 125899):', sum(v for k, v in fc.items() if k in m40.FUNCTIONAL_EFFECTS), flush=True)

pbc = defaultdict(list); ibs = defaultdict(list)
for i, r in enumerate(rows):
    pbc[r['Chrom']].append(int(r['Pos'])); ibs[(r['Chrom'], int(r['Pos']))].append(i)
for ch in pbc:
    pbc[ch] = sorted(set(pbc[ch]))
calls, pst = m41.parse_paf_calls(R + '/06_dSNP_phoenix_polarity/alignments/Phoenix_vs_Africa_hap2.asm20.cs.primary.paf',
                                 rows, pbc, ibs, 20)
for i, r in enumerate(rows):
    r['phx'] = m41.summarize_call(r, calls.get(i, []))['ALT_Polarity_Phoenix'] == 'ALT_Derived'
del calls
print('check phoenix-derived functional n38==1 (expect 70269):',
      sum(1 for r in rows if r['n38'] == 1 and r['phx'] and r['FC'] in m40.FUNCTIONAL_EFFECTS), flush=True)

MZ = ['meizhou4_hap1', 'meizhou4_hap2']
mz = {h: pickle.load(open(f'{W0}/mz4/{h}.pkl', 'rb')) for h in MZ}


def inblk(b, p0):
    if b is None or not len(b):
        return False
    i = np.searchsorted(b[:, 1], p0, side='right')
    return i < len(b) and b[i, 0] <= p0


for r in rows:
    s38 = set(r['Samples'].split(';')); p = int(r['Pos']); ch = r['Chrom']; ok = True
    for h in MZ:
        if h in s38:
            ok = False; break
        d = mz[h]
        if not inblk(d['blocks'].get(ch), p - 1) or inblk(d['del'].get(ch), p - 1):
            ok = False; break
        if d['snp'].get(ch, {}).get(p) is not None:
            ok = False; break
    r['ole'] = ok
print('check strict functional n35==1 (expect 27918):',
      sum(1 for r in rows if r['n35'] == 1 and r['phx'] and r['ole'] and r['FC'] in m40.FUNCTIONAL_EFFECTS), flush=True)
cols = ['Chrom', 'Pos', 'Ref', 'Alt', 'n38', 'n35', 'Samples', 'Samples_African35', 'FC', 'FE', 'phx', 'ole']
with gzip.open(W + '/snv_classes.tsv.gz', 'wt') as out:
    w = csv.DictWriter(out, fieldnames=cols, delimiter='\t', extrasaction='ignore', lineterminator='\n')
    w.writeheader(); w.writerows(rows)
c2 = Counter()
for r in rows:
    for tag, ok in (('All38_n38eq1_phx', r['n38'] == 1 and r['phx']),
                    ('EG35_n35eq1_phx_ole', r['n35'] == 1 and r['phx'] and r['ole'])):
        if ok:
            c2[(tag, r['FC'])] += 1
for k, v in sorted(c2.items()):
    print('count', k[0], k[1], v)
print('DONE', flush=True)
