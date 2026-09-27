#!/usr/bin/env python3
"""Exact display summaries of F006: no calling, smoothing or ancestry inference."""
import csv
import gzip
import hashlib
import json
import sys
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
ORIGINAL = ROOT / 'results/10-ancestry/figures_v3_kmer_pedigree/F006_Kmer_Pedigree_Linear'
OUT = Path(sys.argv[1]) if len(sys.argv)>1 else ROOT / 'results/10-ancestry/F006_Kmer_Pedigree_Linear_v2'
ROWS = ['Dura_h1', 'Dura_h2', 'Pisifera_h1', 'Pisifera_h2', 'MZ4_h1', 'MZ4_h2',
        'EO12', 'EG11', 'Nigerian_h1', 'Nigerian_h2', 'TN_h1', 'TN_h2', 'FL_HapA', 'FL_HapB']
CLASSES = ['Dura-like', 'Pisifera-like', 'Oleifera-like', 'Mixed', 'Insufficient']
GROUPS = ['Dura', 'Dura', 'Pisifera', 'Pisifera', 'E. oleifera', 'E. oleifera',
          'EO12', 'EG11', 'Nigerian', 'Nigerian', 'TN', 'TN', 'FL', 'FL']


def write(name, header, records):
    with (OUT / 'source-data' / name).open('x') as f:
        w = csv.writer(f, delimiter='\t', lineterminator='\n')
        w.writerow(header)
        w.writerows(records)


def main():
    OUT.mkdir(exist_ok=False)
    (OUT / 'source-data').mkdir()
    (OUT / 'checks').mkdir()
    length_path = ROOT / 'config/10-ancestry/kmer_pedigree/Chrom_Lengths.tsv'
    with length_path.open() as f:
        lengths = {r['Chr']: int(r['Length']) for r in csv.DictReader(f, delimiter='\t')}
    genome = sum(lengths.values())
    calls = defaultdict(list)
    source = ORIGINAL / 'source-data/Consensus.tsv.gz'
    with gzip.open(source, 'rt') as f:
        for r in csv.DictReader(f, delimiter='\t'):
            s, c, a, b, state = r['Lookup_ID'], r['Chromosome'], int(r['Start0']), int(r['End0']), r['Consensus']
            assert s in ROWS and c in lengths and state in CLASSES
            assert 0 <= a < b <= lengths[c]
            calls[s, c].append((a, b, state))
    assert len(calls) == 14 * 16
    counts = {s: Counter() for s in ROWS}
    bins = defaultdict(Counter)
    runs, windows = [], 0
    for s in ROWS:
        for c, length in lengths.items():
            cursor = 0
            previous = None
            for a, b, state in sorted(calls[s, c]):
                assert a == cursor, (s, c, a, cursor)
                cursor = b
                windows += 1
                counts[s][state] += b-a
                if previous and previous[4] == state and previous[3] == a:
                    previous[3] = b
                else:
                    previous = [s, c, a, b, state]
                    runs.append(previous)
                pos = a
                while pos < b:
                    tile = pos // 2000000
                    stop = min(b, (tile+1)*2000000)
                    bins[c, tile][state] += stop-pos
                    pos = stop
            assert cursor == length
        assert sum(counts[s].values()) == genome
    write('Track_Runs_v2.tsv', ['Sample', 'Chromosome', 'Start0', 'End0', 'Class'], runs)
    write('Chromosome_Lengths_v2.tsv', ['Chromosome', 'Length_bp'], lengths.items())
    write('Track_Composition_v2.tsv', ['Sample', 'Group', 'Class', 'Length_bp', 'Denominator_bp', 'Fraction'],
          [[s, g, k, counts[s][k], genome, counts[s][k]/genome]
           for s, g in zip(ROWS, GROUPS) for k in CLASSES])
    binrows = []
    for c, length in lengths.items():
        for tile in range((length+1999999)//2000000):
            a, b = tile*2000000, min((tile+1)*2000000, length)
            den = (b-a)*len(ROWS)
            assert sum(bins[c, tile].values()) == den
            binrows.extend([c, a, b, k, bins[c, tile][k], den, bins[c, tile][k]/den] for k in CLASSES)
    write('Local_Composition_v2.tsv', ['Chromosome', 'Start0', 'End0', 'Class', 'Length_bp', 'Denominator_bp', 'Fraction'], binrows)
    group_sets = {'Dura': ROWS[:2], 'Pisifera': ROWS[2:4], 'E. oleifera': ROWS[4:6], 'Other tracks': ROWS[6:]}
    write('Group_Composition_v2.tsv', ['Group', 'Class', 'Length_bp', 'Denominator_bp', 'Fraction', 'N_tracks'],
          [[g, k, sum(counts[s][k] for s in samples), genome*len(samples),
            sum(counts[s][k] for s in samples)/(genome*len(samples)), len(samples)]
           for g, samples in group_sets.items() for k in CLASSES])
    check = {'Status': 'PASS', 'Windows': windows, 'Exact_contiguous_runs': len(runs),
             'Tracks': len(ROWS), 'Chromosomes': len(lengths), 'Per_track_denominator_bp': genome,
             'Original_consensus_sha256': hashlib.sha256(source.read_bytes()).hexdigest(),
             'Checks': ['coordinates bounded', 'full grid; no gaps or overlaps',
                        'all five categories retained', 'exact contiguous-state compression only',
                        'all compositions sum to recorded denominators'],
             'Local_summary': 'Disjoint 2-Mb bins; bp-weighted across 14 display tracks; no smoothing',
             'Group_summary': 'Pooled display-track composition, not population ancestry', 'Calibrated': False}
    (OUT / 'checks/Data_Validation_v2.json').write_text(json.dumps(check, indent=2)+'\n')
    print(json.dumps(check, indent=2))


if __name__ == '__main__':
    main()
