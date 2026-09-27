#!/usr/bin/env python3
"""Draw the new hap1 intervals in the original Figure 2a Oleifera cell style."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import shutil

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Rectangle
from matplotlib.text import Text

# Infer identities from the old Oleifera cell's single-colour chromosomes and
# saved old interval table, rather than the order of the five bars above the tree.
COLOURS = {'APK1': '#f8d093', 'APK2': '#dc6372', 'APK3': '#758bc4',
           'APK4': '#6cbee6', 'APK5': '#9e9f91'}


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run', type=Path, required=True)
    parser.add_argument('--reference-image', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    run, reference = args.run.resolve(strict=True), args.reference_image.resolve(strict=True)
    output = args.output.resolve()
    output.mkdir(exist_ok=False)
    source = run / 'karyotype_segments.tsv'
    with source.open() as handle:
        rows = list(csv.DictReader(handle, delimiter='\t'))
    with (run / 'model_counts.tsv').open() as handle:
        counts = {k: int(v) for k, v in next(csv.DictReader(handle, delimiter='\t')).items()}
    lengths = {r[0]: int(r[2]) for r in
               (line.split() for line in (run / 'Eoleifera.len').read_text().splitlines())}
    chromosomes = ['chr{:02d}A_RagTag'.format(i) for i in range(1, 17)]
    assert list(lengths) == chromosomes
    assert counts == {'A': 10, 'B': 38, 'C': 16, 'Model_fissions': 28, 'Model_fusions': 22}
    assert len(rows) == counts['B']
    for chrom in chromosomes:
        end = 0
        for row in [r for r in rows if r['Chromosome'] == chrom]:
            start = int(row['Start_gene_order'])
            assert start == end + 1
            end = int(row['End_gene_order'])
        assert end == lengths[chrom]

    # Coordinates follow the reference cell at a 1600-pixel full-figure display
    # width. Treat each design unit as one PDF point to preserve all proportions.
    width, height = 101.0, 218.0
    bar_width, pitch, first_x = 4.75, 10.15, 16.2
    baselines = [112.8, 198.4]
    scale = 83.0 / max(lengths.values())
    stem = 'Oleifera_hap1.reference_style'
    plt.rcParams.update({'font.family': 'Liberation Sans', 'font.size': 10,
                         'pdf.fonttype': 42, 'ps.fonttype': 42, 'svg.fonttype': 'none'})
    fig = plt.figure(figsize=(width / 72, height / 72))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set(xlim=(0, width), ylim=(height, 0))
    ax.set_axis_off()
    frame = FancyBboxPatch((6.8, 5.0), 91.2, 209.0,
                           boxstyle='round,pad=0,rounding_size=10', facecolor='none',
                           edgecolor='#737373', linewidth=.7, linestyle=(0, (4.8, 4.3)),
                           capstyle='round')
    frame.set_gid('species_frame')
    ax.add_patch(frame)
    ax.text(49.1, 18.0, 'Oleifera', ha='center', va='center',
            fontsize=11, fontweight='bold', color='#ce5b18')
    rectangles = []
    for chrom_index, chrom in enumerate(chromosomes):
        row_index, column = divmod(chrom_index, 8)
        x, baseline = first_x + pitch * column, baselines[row_index]
        for row in [r for r in rows if r['Chromosome'] == chrom]:
            start, end = int(row['Start_gene_order']), int(row['End_gene_order'])
            segment_height = (end - start + 1) * scale
            patch = Rectangle((x - bar_width / 2, baseline - end * scale),
                              bar_width, segment_height, facecolor=COLOURS[row['APK_component']],
                              edgecolor='none', linewidth=0)
            patch.set_gid('segment_{:02d}_{}_{}'.format(chrom_index + 1, start, row['APK_component']))
            ax.add_patch(patch)
            rectangles.append({'chromosome': chrom, 'start': start, 'end': end,
                               'component': row['APK_component'], 'height_points': segment_height})
    for y, start_label, end_label, arrow_start, arrow_end, label_start, label_x in [
            (120.3, '1', '8', 18.2, 84.4, 14.5, 87.8),
            (204.6, '9', '16', 19.5, 79.0, 15.6, 88.7)]:
        ax.add_patch(FancyArrowPatch((arrow_start, y), (arrow_end, y),
                     arrowstyle='-|>,head_length=6.5,head_width=3.2', mutation_scale=1,
                     linewidth=.8, color='#111111', shrinkA=0, shrinkB=0))
        ax.text(label_start, y + .7, start_label, ha='center', va='center', fontsize=10)
        ax.text(label_x, y + .7, end_label, ha='center', va='center', fontsize=10)

    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    labels = [(t.get_text(), t.get_window_extent(renderer)) for t in fig.findobj(Text)
              if t.get_visible() and t.get_text().strip()]
    clipped = [t for t, box in labels if box.x0 < 0 or box.y0 < 0 or
               box.x1 > fig.bbox.width or box.y1 > fig.bbox.height]
    overlaps = [(a, b) for i, (a, box) in enumerate(labels)
                for b, other in labels[i + 1:] if box.overlaps(other)]
    assert not clipped and not overlaps
    assert len(rectangles) == 38
    for rect in rectangles:
        assert abs(rect['height_points'] / (rect['end'] - rect['start'] + 1) - scale) < 1e-12
    for extension in ['pdf', 'svg', 'png']:
        fig.savefig(output / (stem + '.' + extension), dpi=900, facecolor='white')
    plt.close(fig)
    shutil.copy2(source, output / (stem + '_plotting_data.tsv'))
    shutil.copy2(run / 'model_counts.tsv', output / 'model_counts.tsv')
    shutil.copy2(__file__, output / 'plot_oleifera_reference_style.py')
    metadata = {
        'status': 'candidate_for_author_review', 'material': 'meizhou4 hap1',
        'display_label': 'Oleifera', 'reference_image': str(reference),
        'reference_image_sha256': sha256(reference),
        'reference_cell_crop_original_pixels': [807, 800, 1042, 1307],
        'reference_measurement_scale': '3720-pixel source width / 1600-pixel display width = 2.325',
        'source_data': str(source), 'source_data_sha256': sha256(source),
        'code_sha256': sha256(__file__), 'counts': counts, 'palette': COLOURS,
        'palette_basis': {'APK1': 'Old Chr8, yellow', 'APK2': 'Old Chr4, pink-red',
                          'APK3': 'Old Chr14, blue', 'APK4': 'Old Chr3, cyan',
                          'APK5': 'Old Chr13 dominant grey'},
        'palette_correction': 'V02 assigned APK3/4/5 grey/blue/cyan from the top schematic. '
                              'This version follows the actual Oleifera cell: blue/cyan/grey. '
                              'Component identities and intervals are unchanged.',
        'figure_spec': {'profile': 'user_reference_cell', 'width_mm': width * 25.4 / 72,
                        'height_mm': height * 25.4 / 72, 'font': 'Liberation Sans (Arial-compatible)',
                        'title_points': 11, 'range_label_points': 10, 'dpi': 900,
                        'bar_width_points': bar_width, 'bar_pitch_points': pitch,
                        'row_baselines_points': baselines, 'points_per_gene': scale,
                        'maximum_bar_height_points': 83, 'frame': 'grey rounded dashed border'},
        'metric_spec': {'unit': 'retained coding-gene order', 'source_coordinates': '1-based inclusive',
                        'display_coordinates': 'start-1/end boundaries',
                        'row_scales': 'identical points per gene in both rows; no per-chromosome stretching'},
        'pattern': 'Specialized chromosome interval tracks with user-requested reference-cell decoration',
        'design_comparison': 'Restore reference cell aspect ratio, thin columns, row spacing, orange label, '
                             'grey dashed frame, endpoint labels and black arrows. No Y axis or extra legend '
                             'inside the tile because the main figure supplies shared context.',
        'text_geometry': {'items': len(labels), 'clipped': clipped, 'overlaps': overlaps},
        'data_rectangles': len(rectangles),
        'limitations': ['Six segments have low matching-seed fractions; previous scientific review still applies',
                       'Fixed APK reference and assumed A=10; not historical event counts',
                       'Legacy palette retained by user request; monochrome discrimination is limited',
                       'Existing composite image has not been overwritten']}
    (output / (stem + '_metadata.json')).write_text(json.dumps(metadata, indent=2) + '\n')
    (output / (stem + '_notes.md')).write_text(
        '# Oleifera reference-style tile\n\n'
        'The user requested the style, proportions and colours of the Oleifera cell in Figure2-A4_1.5fold.jpg. '
        'The new drawing reproduces that cell geometry using the current meizhou4 hap1 data. '
        'Both rows have the same points-per-gene scale. Bars encode retained gene order, not base pairs. '
        'Only endpoint labels 1/8 and 9/16 are shown, matching the reference.\n\n'
        'The actual reference-cell palette maps APK1/2/3/4/5 to yellow/pink-red/blue/cyan/grey. '
        'This corrects the APK3/4/5 display permutation in V02 without changing data or identities. '
        'There are still 38 blocks and model counts of 28 fissions / 22 fusions at A=10, C=16. '
        'Counts and the shared APK legend belong outside this insertion tile in the composite. '
        'The original composite and previous figure versions remain unchanged.\n')
    with (output / 'file_manifest.tsv').open('w', newline='') as handle:
        writer = csv.writer(handle, delimiter='\t'); writer.writerow(['File', 'SHA256'])
        for path in sorted(output.iterdir()):
            if path.is_file() and path.name != 'file_manifest.tsv':
                writer.writerow([path.name, sha256(path)])
    print(json.dumps({'output': str(output), 'blocks': len(rectangles),
                      'text_checks': metadata['text_geometry'], 'palette': COLOURS}))


if __name__ == '__main__':
    main()
