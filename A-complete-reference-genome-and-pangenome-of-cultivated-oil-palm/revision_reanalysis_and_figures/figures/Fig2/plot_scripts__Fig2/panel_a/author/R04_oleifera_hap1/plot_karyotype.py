#!/usr/bin/env python3
"""Render the reviewed-format candidate from one run's interval and count tables."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import shutil

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch, Rectangle
from matplotlib.ticker import MaxNLocator

# Match component identities to the actual Oleifera cell in Figure2-A4_1.5fold.jpg.
PALETTE = {'APK1': '#f8d093', 'APK2': '#dc6372', 'APK3': '#758bc4',
           'APK4': '#6cbee6', 'APK5': '#9e9f91'}


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run', type=Path, required=True)
    parser.add_argument('--output', type=Path, help='New figure-version directory; must not exist')
    args = parser.parse_args()
    run = args.run.resolve(strict=True)
    output = args.output.resolve() if args.output else run / 'figures'
    output.mkdir(exist_ok=False)
    table = run / 'karyotype_segments.tsv'
    with table.open() as handle:
        rows = list(csv.DictReader(handle, delimiter='\t'))
    with (run / 'model_counts.tsv').open() as handle:
        counts = {k: int(v) for k, v in next(csv.DictReader(handle, delimiter='\t')).items()}
    chroms = [line.split()[0] for line in (run / 'Eoleifera.len').read_text().splitlines()]
    assert chroms == ['chr{:02d}A_RagTag'.format(i) for i in range(1, 17)]
    lengths = {line.split()[0]: int(line.split()[2])
               for line in (run / 'Eoleifera.len').read_text().splitlines()}
    assert counts['A'] == 10 and counts['B'] == len(rows) and counts['C'] == 16
    assert counts['A'] + counts['Model_fissions'] - counts['Model_fusions'] == counts['C']
    for chrom in chroms:
        end = 0
        for row in [r for r in rows if r['Chromosome'] == chrom]:
            assert int(row['Start_gene_order']) == end + 1
            end = int(row['End_gene_order'])
            assert row['APK_component'] in PALETTE
        assert end == lengths[chrom]
    stem = 'Eoleifera_hap1.karyotype'
    shutil.copy2(table, output / (stem + '_plotting_data.tsv'))
    shutil.copy2(run / 'model_counts.tsv', output / 'model_counts.tsv')
    shutil.copy2(__file__, output / 'plot_karyotype.py')
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 8,
                         'pdf.fonttype': 42, 'ps.fonttype': 42, 'svg.fonttype': 'none',
                         'axes.linewidth': .4, 'xtick.major.width': .4, 'ytick.major.width': .4})

    text_checks = []
    row_maxima = [max(lengths[c] for c in chroms[:8]), max(lengths[c] for c in chroms[8:])]
    for clean in (False, True):
        fig, axes = plt.subplots(2, 1, figsize=(3.5, 3.6) if clean else (7.0, 4.5),
                                 gridspec_kw={'height_ratios': row_maxima})
        fig.subplots_adjust(left=.04 if clean else .10, right=.98,
                            bottom=.08 if clean else .17, top=.91, hspace=.20 if clean else .26)
        for group_index, ax in enumerate(axes):
            group = chroms[group_index * 8:(group_index + 1) * 8]
            for x, chrom in enumerate(group):
                for row in [r for r in rows if r['Chromosome'] == chrom]:
                    start, end = int(row['Start_gene_order']), int(row['End_gene_order'])
                    ax.add_patch(Rectangle((x - .28, start - 1), .56, end - start + 1,
                                           facecolor=PALETTE[row['APK_component']], edgecolor='none'))
                ax.add_patch(Rectangle((x - .28, 0), .56, lengths[chrom],
                                       facecolor='none', edgecolor='#555555', linewidth=.3))
            ax.set(xlim=(-.7, 7.7), ylim=(0, row_maxima[group_index] * 1.05))
            ax.set_xticks(range(8), [str(i + 1 + group_index * 8) for i in range(8)])
            ax.tick_params(axis='x', length=0, pad=4)
            for side in ('top', 'right', 'bottom'):
                ax.spines[side].set_visible(False)
            if clean:
                ax.spines['left'].set_visible(False)
                ax.set_yticks([])
            else:
                ax.set_ylabel('Gene order', fontsize=8)
                ax.yaxis.set_major_locator(MaxNLocator(nbins=4, integer=True))
                ax.tick_params(axis='y', labelsize=7, length=2)
        fig.suptitle('Elaeis oleifera (hap1)', fontsize=9, fontstyle='italic', y=.98)
        if not clean:
            fig.legend(handles=[Patch(facecolor=colour, label=component)
                                for component, colour in PALETTE.items()],
                       loc='lower center', bbox_to_anchor=(.5, .067), ncol=5, frameon=False,
                       fontsize=8, handlelength=1.1, handleheight=.8, columnspacing=1.7)
            fig.text(.5, .025,
                     'ACPK = 10    Blocks = {B}    Model fissions = {Model_fissions}    '
                     'Model fusions = {Model_fusions}'.format(**counts), ha='center', fontsize=8)
        name = stem + ('.clean' if clean else '')
        # Different row heights remove empty coordinate space without changing pixels/gene.
        scales = [ax.get_position().height / ax.get_ylim()[1] for ax in axes]
        assert abs(scales[0] - scales[1]) < 1e-12
        fig.canvas.draw()
        renderer = fig.canvas.get_renderer()
        labels = [(text.get_text(), text.get_window_extent(renderer))
                  for text in fig.findobj(matplotlib.text.Text)
                  if text.get_visible() and text.get_text().strip()]
        clipped = [text for text, box in labels if box.x0 < 0 or box.y0 < 0
                   or box.x1 > fig.bbox.width or box.y1 > fig.bbox.height]
        overlaps = [(a, b) for i, (a, box_a) in enumerate(labels)
                    for b, box_b in labels[i + 1:] if box_a.overlaps(box_b)]
        assert not clipped, 'Text outside figure: ' + repr(clipped)
        assert not overlaps, 'Overlapping text: ' + repr(overlaps)
        text_checks.append({'figure': name, 'text_items': len(labels), 'clipped_text': clipped,
                            'overlapping_text': overlaps, 'normalized_height_per_gene': scales})
        for suffix in ('pdf', 'svg', 'png'):
            fig.savefig(output / (name + '.' + suffix), dpi=600, facecolor='white')
        plt.close(fig)

    metadata = {
        'material': 'meizhou4 hap1', 'species_label': 'Elaeis oleifera',
        'status': 'manuscript_candidate_requires_visual_and_scientific_review',
        'input': str(table), 'input_sha256': digest(table), 'code_sha256': digest(__file__),
        'coordinate_convention': '1-based inclusive gene-order intervals; drawn on start-1,end boundaries',
        'chromosome_order': chroms, 'chromosome_labels': list(range(1, 17)),
        'length_semantics': 'Number of retained coding representatives; identical scale in both rows',
        'counts': counts, 'palette': PALETTE,
        'palette_source': 'Oleifera cell in retained Figure2-A4_1.5fold.jpg; matched using old single-component chromosomes and JPEG RGB samples',
        'figure_spec': {'journal_profile': 'general_scientific', 'width_mm': 177.8,
                        'height_mm': 114.3, 'dpi': 600, 'font': 'DejaVu Sans',
                        'clean_inset_width_mm': 88.9, 'clean_inset_height_mm': 91.44},
        'pattern_selection': 'Specialized chromosome interval tracks; native rectangle geometry. '
                             'Manhattan pattern inspected but not used because there is no score axis.',
        'cognitive_load_exception': 'Five reference groups and 16 chromosomes are required; '
                                    'one shared legend, two equally scaled rows of eight chromosomes. '
                                    'Row heights are proportional to their maximum gene counts, preserving '
                                    'the same physical scale while removing empty coordinate space.',
        'text_geometry_checks': text_checks,
        'comparison_scope': 'Old material is different; preserve layout/colour identity, not old data values.',
        'limitations': ['RagTag reference-assisted scaffolding using old American_hap1',
                       'Liftoff annotation transfer and recorded protein QC exclusions',
                       'Existing five-component APK is a fixed reference, not newly independently inferred',
                       'WGDI boundary extension/recolouring; model counts are not historical events']}
    (output / (stem + '_metadata.json')).write_text(json.dumps(metadata, indent=2) + '\n')
    (output / (stem + '_notes.md')).write_text(
        '# meizhou4 hap1 chromosome composition\n\n'
        'The annotated figure uses two rows of eight chromosomes at the same gene-order scale. '
        'Row heights follow their maximum gene counts; pixels per gene are identical in both rows. '
        'The clean version is an insertion tile; use it with a shared APK legend and gene-order '
        'scale explanation in the composite caption. Chromosome labels 1–16 correspond to '
        'chr01A_RagTag–chr16A_RagTag. The plotting table is the exact interval source.\n\n'
        'Only translated representative proteins passing the declared QC contribute to gene order. '
        'Counts use the complete split-then-join model at A=10, not independently inferred events. '
        'Scaffold and annotation transfer dependence on the older American reference remains a '
        'limitation. No claim about a change from the old material follows from swapping this panel.\n')
    with (output / 'file_manifest.tsv').open('w', newline='') as handle:
        writer = csv.writer(handle, delimiter='\t')
        writer.writerow(['File', 'SHA256'])
        for p in sorted(output.iterdir()):
            if p.is_file() and p.name != 'file_manifest.tsv':
                writer.writerow([p.name, digest(p)])
    print(json.dumps({'figures': str(output), 'counts': counts}))


if __name__ == '__main__':
    main()
