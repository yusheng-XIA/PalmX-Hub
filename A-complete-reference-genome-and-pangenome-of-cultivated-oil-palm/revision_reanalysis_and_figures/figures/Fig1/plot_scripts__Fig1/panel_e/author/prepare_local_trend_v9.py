#!/usr/bin/env python3
"""Exact bp-weighted local display trends; original ancestry calls are immutable."""
import bisect
import csv
import hashlib
import json
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
PACKAGE = ROOT / 'results/10-ancestry/F006_Kmer_Pedigree_Linear_v2'
STATES = ['Dura-like', 'Pisifera-like', 'Oleifera-like', 'Mixed', 'Insufficient']


def read(path):
    with path.open() as f:
        return list(csv.DictReader(f, delimiter='\t'))


def main():
    source = PACKAGE / 'source-data/Track_Runs_v2.tsv'
    lengths = {r['Chromosome']: int(r['Length_bp']) for r in
               read(PACKAGE / 'source-data/Chromosome_Lengths_v2.tsv')}
    events = defaultdict(lambda: defaultdict(lambda: [0]*5))
    samples = set()
    for row in read(source):
        samples.add(row['Sample'])
        ch, a, b, k = row['Chromosome'], int(row['Start0']), int(row['End0']), STATES.index(row['Class'])
        events[ch][a][k] += 1
        events[ch][b][k] -= 1
    assert len(samples) == 14
    result = []
    for ch, length in lengths.items():
        boundaries = sorted(events[ch])
        assert boundaries[0] == 0 and boundaries[-1] == length
        active = [0]*5
        area = [0]*5
        prefix, levels = [], []
        previous = 0
        for pos in boundaries:
            area = [area[k] + active[k]*(pos-previous) for k in range(5)]
            active = [active[k] + events[ch][pos][k] for k in range(5)]
            prefix.append(area.copy())
            levels.append(active.copy())
            assert min(active) >= 0 and sum(active) == (0 if pos == length else 14)
            previous = pos
        assert sum(area) == length*14

        def integral(pos):
            i = bisect.bisect_right(boundaries, pos)-1
            return [prefix[i][k] + (pos-boundaries[i])*levels[i][k] for k in range(5)]

        positions = list(range(0, length, 2000000)) + [length]
        for center in positions:
            a, b = max(0, center-15000000), min(length, center+15000000)
            left, right = integral(a), integral(b)
            amounts = [right[k]-left[k] for k in range(5)]
            denominator = 14*(b-a)
            assert sum(amounts) == denominator
            result.extend([ch, center, a, b, STATES[k], amounts[k], denominator,
                           amounts[k]/denominator] for k in range(5))
    output = PACKAGE / 'source-data/Local_Trend_30Mb_v9.tsv'
    with output.open('x') as f:
        writer = csv.writer(f, delimiter='\t', lineterminator='\n')
        writer.writerow(['Chromosome', 'Center_bp', 'Window_Start0', 'Window_End0',
                         'Class', 'Length_bp', 'Denominator_bp', 'Fraction'])
        writer.writerows(result)
    report = {'Status': 'PASS', 'Display_only': True, 'Main_calls_changed': False,
              'Window_bp': 30000000, 'Step_bp': 2000000, 'Tracks': 14,
              'Chromosomes': 16, 'Positions': len(result)//5,
              'Boundary_rule': 'Centered windows clipped at chromosome ends; never cross chromosomes',
              'Weighting': 'Exact overlap of original constant-state intervals with each window; all five classes retained',
              'Source_SHA256': hashlib.sha256(source.read_bytes()).hexdigest(),
              'Output_SHA256': hashlib.sha256(output.read_bytes()).hexdigest()}
    (PACKAGE / 'checks/Local_Trend_Validation_v9.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
