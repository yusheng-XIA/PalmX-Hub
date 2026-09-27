#!/usr/bin/env python3
import argparse
import csv
from collections import Counter
from pathlib import Path

ap = argparse.ArgumentParser()
ap.add_argument('--tier1', required=True)
ap.add_argument('--syri-large', required=True)
ap.add_argument('--fai', required=True)
ap.add_argument('--out-dir', required=True)
args = ap.parse_args()
outdir = Path(args.out_dir)
outdir.mkdir(parents=True, exist_ok=True)

def read_rows(path):
    with open(path) as h:
        return list(csv.DictReader(h, delimiter='\t'))

tier1 = read_rows(args.tier1)
syri = read_rows(args.syri_large)
strict = [r for r in tier1 if r['SVTYPE'] in {'DEL','INS','INV','DUP','TRA'}]
mixed = []
for r in tier1:
    if r['SVTYPE'] in {'DEL','INS'}:
        x = dict(r); x['Evidence_Layer'] = 'Tier1_cuteSV_AND_assembly'; mixed.append(x)
for r in syri:
    if r['SVTYPE'] in {'INV','DUP','TRA'}:
        x = dict(r); x['Evidence_Layer'] = 'SyRI_large_rearrangement'; mixed.append(x)

chrom_rank = {}
contigs = []
with open(args.fai) as h:
    for i, line in enumerate(h):
        f = line.rstrip().split('\t')
        chrom_rank[f[0]] = i; contigs.append((f[0], f[1]))
mixed.sort(key=lambda r: (chrom_rank.get(r['Chrom'], 999), int(r['Start']), r['SVTYPE']))

def write_tsv(path, rows, add_layer=False):
    if not rows: raise SystemExit(f'no rows for {path}')
    fields = list(rows[0].keys())
    if add_layer and 'Evidence_Layer' not in fields: fields.append('Evidence_Layer')
    with open(path, 'w', newline='') as out:
        w = csv.DictWriter(out, fieldnames=fields, delimiter='\t', extrasaction='ignore')
        w.writeheader(); w.writerows(rows)

write_tsv(outdir / 'oilpalm_hap38.tier1_strict_sv.tsv', strict)
write_tsv(outdir / 'oilpalm_hap38.final_highconfidence_sv.tsv', mixed, True)

with open(outdir / 'oilpalm_hap38.final_highconfidence_sv.vcf', 'w') as out:
    out.write('##fileformat=VCFv4.2\n')
    out.write('##source=oilpalm_hap39_incremental_multievidence_v1\n')
    for c, n in contigs: out.write(f'##contig=<ID={c},length={n}>\n')
    out.write('##INFO=<ID=SVTYPE,Number=1,Type=String,Description="SV type">\n')
    out.write('##INFO=<ID=END,Number=1,Type=Integer,Description="End position">\n')
    out.write('##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Median SV length">\n')
    out.write('##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of non-reference assemblies">\n')
    out.write('##INFO=<ID=FREQ,Number=1,Type=Float,Description="Assembly frequency among 38 non-reference assemblies">\n')
    out.write('##INFO=<ID=EVIDENCE,Number=1,Type=String,Description="Evidence layer">\n')
    out.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n')
    for i, r in enumerate(mixed, 1):
        svlen = int(float(r['SVLEN_Median_bp']))
        if r['SVTYPE'] == 'DEL': svlen = -abs(svlen)
        info = f"SVTYPE={r['SVTYPE']};END={r['End']};SVLEN={svlen};NS={r['Sample_Count']};FREQ={float(r['Frequency']):.6f};EVIDENCE={r['Evidence_Layer']}"
        out.write(f"{r['Chrom']}\t{r['Start']}\tOP39SV{i:07d}\tN\t<{r['SVTYPE']}>\t.\tPASS\t{info}\n")

counts = Counter(r['SVTYPE'] for r in mixed)
with open(outdir / 'oilpalm_hap38.final_highconfidence_sv.summary.tsv', 'w') as out:
    out.write('SVTYPE\tCount\n')
    for t in ('DEL','INS','INV','DUP','TRA'): out.write(f'{t}\t{counts[t]}\n')
    out.write(f'TOTAL\t{sum(counts.values())}\n')
print(f'[OK] strict={len(strict)} mixed_final={len(mixed)} counts={dict(counts)}')
