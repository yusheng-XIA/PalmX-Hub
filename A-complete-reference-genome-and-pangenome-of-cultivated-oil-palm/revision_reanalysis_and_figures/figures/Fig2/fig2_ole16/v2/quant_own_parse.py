#!/usr/bin/env python3
"""Parse own FL/TN libraries with exactly the ole16_cis/quant_parse.py counting rules
(alignment >=60 bp, NM <= 3% of aligned length, primary alignment assigned, RPM over sampled read-1 records)."""
import os, re, glob, collections, csv
C = '../../ole16_cis'
fam = {}; uid = {}
for l in open(f'{C}/target_family.tsv'):
    k, u, f = l.rstrip('\n').split('\t'); fam[k] = f; uid[k] = u
meta = {}
for l in open('own_libs.tsv'):
    if l.startswith('lib'): continue
    lib, mat, st = l.rstrip('\n').split('\t'); meta[lib] = (mat, st)
OLEa = {'African_hap2__evm.TU.chr11B.1497': 'Eg', 'American_hap1__evm.TU.chr11A.1059': 'Eo'}
rows = []
for sam in sorted(glob.glob('own/own_rq/*.sam')):
    r = os.path.basename(sam)[:-4]
    T = int(open(f'own/own_hits/{r}.count').read().strip())
    best = {}; alle = collections.defaultdict(dict)
    for l in open(sam):
        f = l.rstrip('\n').split('\t'); q, flag, t = f[0], int(f[1]), f[2]
        if t == '*' or flag & 4: continue
        alen = sum(int(n) for n, o in re.findall(r'(\d+)([MI=X])', f[5]))
        tags = dict((x[:2], x[5:]) for x in f[6:] if len(x) > 5)
        nm = int(tags.get('NM', 99)); AS = int(tags.get('AS', -999))
        if alen < 60 or nm > 0.03 * alen: continue
        if t in OLEa: alle[q][OLEa[t]] = max(AS, alle[q].get(OLEa[t], -999))
        if flag & 256 or flag & 2048: continue
        best[q] = t
    cnt = collections.Counter(uid.get(t, t) for t in best.values())
    famc = collections.Counter(fam.get(t, '?') for t in best.values())
    eo = sum(1 for q, d in alle.items() if 'Eo' in d and d.get('Eo', -999) > d.get('Eg', -999))
    eg = sum(1 for q, d in alle.items() if 'Eg' in d and d.get('Eg', -999) > d.get('Eo', -999))
    mat, st = meta[r]
    rows.append(dict(run=r, material=mat, stage=st, total_reads=T,
        OLE16a_rpm=round(cnt['UFTN008330'] / T * 1e6, 1), OLE16b_rpm=round(cnt['UFTN008332'] / T * 1e6, 1),
        LDAP_chr10_rpm=round(cnt['UFTN008529'] / T * 1e6, 1), ACTIN_rpm=round(cnt['ACTIN'] / T * 1e6, 1),
        oleosin_all_rpm=round(famc['oleosin'] / T * 1e6, 1), LDAP_all_rpm=round(famc['REF'] / T * 1e6, 1),
        OLE16a_Eo_allele_reads=eo, OLE16a_Eg_allele_reads=eg))
w = csv.DictWriter(open('own_quant.tsv', 'w'), fieldnames=list(rows[0]), delimiter='\t'); w.writeheader(); [w.writerow(x) for x in rows]
for x in rows: print('\t'.join(map(str, x.values())))
