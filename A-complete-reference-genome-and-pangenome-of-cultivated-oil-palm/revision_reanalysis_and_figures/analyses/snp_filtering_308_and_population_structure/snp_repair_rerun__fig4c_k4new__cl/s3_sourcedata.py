#!/usr/bin/env python3
"""Source Data / ED4b replacement TSVs + figure values JSON from pi_fst_genomewide.tsv (recommended = present set)."""
import csv, json, os
B = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
R = list(csv.DictReader(open(os.path.join(B, 'pi_fst_genomewide.tsv')), delimiter='\t'))
POPS = ['K4_Pop1', 'K4_Pop2', 'K4_Pop3', 'K4_Pop4']
PAIRS = [('K4_Pop1','K4_Pop2'),('K4_Pop1','K4_Pop3'),('K4_Pop1','K4_Pop4'),('K4_Pop2','K4_Pop3'),('K4_Pop2','K4_Pop4'),('K4_Pop3','K4_Pop4')]
def get(metric, g, est_prefix):
    for r in R:
        if r['metric'] == metric and r['group_or_pair'] == g and r['snp_set'] == 'present' and r['estimator'].startswith(est_prefix): return r
    raise KeyError(g)
SD = os.path.join(B, 'source_data'); os.makedirs(SD, exist_ok=True)
pi = {p: get('pi', p, 'missing-aware') for p in POPS}
fs = {a + '_' + b: get('fst', a + '_' + b, 'W&C') for a, b in PAIRS}
w = lambda fn: csv.writer(open(os.path.join(SD, fn), 'w', newline=''), delimiter='\t', lineterminator='\n')
t = w('Fig.4c_pi.tsv'); t.writerow(['group', 'PI_weighted', 'n_used_weighted', 'n_windows'])
for p in POPS: t.writerow([p, pi[p]['value'], pi[p]['n_windows'], pi[p]['n_windows']])
t = w('Fig.4c_fst_long.tsv'); t.writerow(['pop1', 'pop2', 'FST_mean', 'n_windows_used', 'n_windows', 'for_reference_WEIGHTED_FST_not_plotted'])
for a, b in PAIRS:
    r = fs[a + '_' + b]; t.writerow([a, b, r['value'], r['n_windows'], r['n_windows'], r['ref_weighted_FST(not plotted)']])
t = w('Fig.4c_fst_matrix.tsv'); t.writerow([''] + POPS)
for a in POPS:
    row = [a]
    for b in POPS:
        row.append(0 if a == b else fs[a + '_' + b]['value'] if a + '_' + b in fs else fs[b + '_' + a]['value'])
    t.writerow(row)
t = w('ED4b_pi.tsv'); t.writerow(['group', 'pi', 'n_windows'])
for p in POPS: t.writerow([p, pi[p]['value'], pi[p]['n_windows']])
t = w('ED4b_fst.tsv'); t.writerow(['pop1', 'pop2', 'FST_mean', 'n_windows'])
for a, b in PAIRS: t.writerow([a, b, fs[a + '_' + b]['value'], fs[a + '_' + b]['n_windows']])
vals = {'pi': {p.replace('K4_', ''): float(pi[p]['value']) * 1e3 for p in POPS},
        'fst': {k.replace('K4_', ''): float(v['value']) for k, v in fs.items()}}
fmin = min(vals['fst'].values()); fmax = max(vals['fst'].values())
vals['ticks'] = [t for t in (0.05, 0.07, 0.09, 0.11, 0.13) if fmin <= t <= fmax]
json.dump(vals, open(os.path.join(B, 'figure4c_values.json'), 'w'), indent=1)
print(json.dumps(vals, indent=1))
