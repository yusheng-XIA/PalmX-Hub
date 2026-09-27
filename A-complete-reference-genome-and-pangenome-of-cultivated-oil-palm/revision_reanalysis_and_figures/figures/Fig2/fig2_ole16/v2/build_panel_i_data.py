#!/usr/bin/env python3
"""Source table for Fig. 2i: OLE16a RPM per library (public + this study), identical counting pipeline."""
import csv, re
C = '../../ole16_cis'
out = []
meta_db = {r['run_accession']: r for r in csv.DictReader(open(f'{C}/prjdb1773.tsv'), delimiter='\t')}
for r in csv.DictReader(open(f'{C}/public_rna_quant.tsv'), delimiter='\t'):
    s = r['sample']
    if re.match(r'K[DPT]', s):
        g, tis, fr, stg, proj = 'EG kernel', 'kernel (endosperm)', {'D': 'dura', 'P': 'pisifera', 'T': 'tenera'}[s[1]], s[2:] + ' months after fertilization', 'PRJDB1773'
    elif re.match(r'M[DPT]', s):
        g, tis, fr, stg, proj = 'EG mesocarp', 'mesocarp', {'D': 'dura', 'P': 'pisifera', 'T': 'tenera'}[s[1]], s[2:] + ' months after fertilization', 'PRJDB1773'
    elif s.startswith('G1-M-'):
        g, tis, fr, stg, proj = 'EG mesocarp', 'mesocarp', 'E. guineensis', s.split('-')[-1] + ' DAP', 'PRJEB11097'
    elif s.startswith('Op-M-'):
        g, tis, fr, stg, proj = 'EO mesocarp', 'mesocarp', 'E. oleifera', s.split('-')[-1] + ' DAP', 'PRJEB11097'
    elif s.startswith('BC'):
        g, tis, fr, stg, proj = 'Backcross mesocarp', 'mesocarp', '(E. oleifera x E. guineensis) x E. guineensis', s.split('-')[-1] + ' (stage as labelled)', 'PRJEB11097'
    else:
        continue
    out.append(dict(Group=g, Material=fr, Tissue=tis, Sample=s, Stage=stg, Run=r['run'], Project=proj,
                    Reads_sampled=r['total_reads'], OLE16a_RPM=r['OLE16a_rpm'], OLE16b_RPM=r['OLE16b_rpm'],
                    LDAP_chr10_RPM=r['LDAP_chr10_rpm'], ACTIN_RPM=r['ACTIN_rpm']))
for r in csv.DictReader(open('own_quant.tsv'), delimiter='\t'):
    out.append(dict(Group=f"{r['material']} mesocarp", Material=r['material'], Tissue='mesocarp', Sample=r['run'],
                    Stage=r['stage'], Run=r['run'], Project='this study', Reads_sampled=r['total_reads'],
                    OLE16a_RPM=r['OLE16a_rpm'], OLE16b_RPM=r['OLE16b_rpm'], LDAP_chr10_RPM=r['LDAP_chr10_rpm'],
                    ACTIN_RPM=r['ACTIN_rpm']))
# backcross: per-individual maximum over sampled stages (as plotted)
bc = {}
for r in out:
    if r['Group'] == 'Backcross mesocarp':
        ind = r['Sample'].split('-')[0]
        bc[ind] = max(bc.get(ind, 0.0), float(r['OLE16a_RPM']))
for r in out:
    r['Plotted'] = 'yes' if r['Group'] != 'Backcross mesocarp' else ''
for ind, v in bc.items():
    for r in out:
        if r['Group'] == 'Backcross mesocarp' and r['Sample'].split('-')[0] == ind and float(r['OLE16a_RPM']) == v and r['Plotted'] == '':
            r['Plotted'] = 'yes (individual maximum)'; break
w = csv.DictWriter(open('Fig2i_OLE16a_public_RNA.tsv', 'w'), fieldnames=list(out[0]), delimiter='\t')
w.writeheader(); [w.writerow(x) for x in out]
import collections
d = collections.defaultdict(list)
for r in out:
    if r['Plotted'].startswith('yes'):
        d[r['Group']].append(float(r['OLE16a_RPM']))
for g, v in d.items():
    print(f'{g}\tn={len(v)}\tOLE16a RPM {min(v):g}-{max(v):g}')
