#!/usr/bin/env python3
"""Export candidate Fig. 2a source data; standard library only, no overwrite.
Run from this project: python3 export_karyotype_source_data.py
This is a format/consistency export, not a WGDI rerun or scientific signoff.
"""
import csv
import hashlib
import json
import re
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SRC = ROOT / 'workingdir/self_blastp_out'
OUT = ROOT / 'results/01-karyotype/source-data'
CHECK = ROOT / 'results/01-karyotype/checks'
INPUTS = {}
TABLES = {}
DICTIONARY = []
WARNINGS = []


def read_rows(name):
    path = SRC / name
    data = path.read_bytes()
    INPUTS[name] = (len(data), hashlib.sha256(data).hexdigest())
    return [(i, line.split()) for i, line in enumerate(data.decode().splitlines(), 1) if line.strip()]


def table(name, fields, rows):
    # Each field has an explicit definition and unit; no multi-line pseudo-headers.
    headers = [x[0] for x in fields]
    assert len(headers) == len(set(headers))
    assert all(len(r) == len(headers) for r in rows)
    TABLES[name] = (headers, rows)
    DICTIONARY.extend((name, *f) for f in fields)


def f(name, description, unit='not_applicable'):
    return name, description, unit


FIG = f('Figure_Panel', 'Linked manuscript figure panel; this export covers Fig. 2a only')
STATUS = f('Review_Status', 'Candidate = consistency export only; Provisional_upstream = unresolved Oleifera provenance; neither is submission signoff')
SOURCE = f('Source_File', 'Original input basename; inputs are in workingdir/self_blastp_out in the owning project')
MATERIAL = f('Material', 'Source material label; Dura and Pisifera are fruit forms of Elaeis guineensis')
SPECIES = f('Species', 'Scientific species name')
APK = f('APK_Component', 'Five-component ancestral gene-order reference identifier; not a resolved duplicated-copy chromosome ID')
COLOR = f('Colour', 'Original WGDI colour token; see APK definition')


def natural(s):
    return [int(t) if t.isdigit() else t.lower() for t in re.split(r'(\d+)', s)]


def main():
    if OUT.exists() or CHECK.exists():
        raise SystemExit('Refusing to overwrite existing source-data/checks directories.')
    ancestor = read_rows('ancestor.txt')
    out_ancestor = read_rows('ancestor_out.txt')
    lens = {r[0]: (int(r[1]), int(r[2])) for _, r in read_rows('APK.len')}
    assert len(ancestor) == len(out_ancestor) == len(lens) == 5
    colour_origin = {r[3]: r[0] for _, r in ancestor}
    colour_apk = {r[3]: 'APK' + r[0] for _, r in out_ancestor}
    assert len(colour_origin) == len(colour_apk) == 5
    assert set(colour_origin) == set(colour_apk)
    apk_to_chr = {colour_apk[col]: chrom for col, chrom in colour_origin.items()}
    phoenix = {}
    for _, r in read_rows('Phoenix.gff'):
        assert len(r) == 7 and r[1] not in phoenix
        phoenix[r[1]] = r
    genes = read_rows('APK.gff')
    assert len(genes) == 7294
    assert len({r[1] for _, r in genes}) == 7294
    counts = Counter()
    orders = defaultdict(list)
    gene_rows = []
    for line, r in genes:
        assert len(r) == 7
        chrom, gene, start, end, strand, order, original = r
        apk = 'APK' + chrom
        p = phoenix[original]
        assert p[0] == apk_to_chr[apk]
        counts[chrom] += 1
        orders[chrom].append(int(order))
        same = (r[2:5] == p[2:5])
        if not same:
            WARNINGS.append(f'{gene}: APK coordinates/strand differ from Phoenix source record')
        gene_rows.append(['Fig. 2a', apk, gene, int(order), int(start), int(end), strand,
                          original, p[0], int(p[5]), int(p[2]), int(p[3]), p[4], 'APK.gff;Phoenix.gff', 'Candidate'])
    for chrom, (_, n) in lens.items():
        assert counts[chrom] == n
        assert sorted(orders[chrom]) == list(range(1, n + 1))
        assert max(int(r[3]) for _, r in genes if r[0] == chrom) == lens[chrom][0]
    definition = []
    for _, r in out_ancestor:
        chrom, start, end, colour, _ = r
        assert int(start) == 1 and int(end) == counts[chrom]
        definition.append(['Fig. 2a', 'APK' + chrom, 'Phoenix dactylifera', colour_origin[colour],
                           counts[chrom], int(start), int(end), colour,
                           'ancestor.txt;ancestor_out.txt;APK.len;APK.gff', 'Candidate'])
    table('Fig_2a_APK_Definition.csv', [FIG, APK, f('Reference_Species', 'Species providing APK gene-order coordinates'),
          f('Reference_Chromosome', 'Phoenix chromosome used to define this APK component'),
          f('Gene_Count', 'Number of retained APK genes', 'genes'),
          f('Start_Gene_Order', 'First APK gene order, 1-based inclusive', 'gene_order'),
          f('End_Gene_Order', 'Last APK gene order, 1-based inclusive', 'gene_order'), COLOR, SOURCE, STATUS], definition)
    table('Fig_2a_APK_Gene_Mapping.csv', [FIG, APK, f('APK_Gene_ID', 'Unique retained APK gene identifier'),
          f('APK_Gene_Order', 'Gene order within APK component, 1-based', 'gene_order'),
          f('APK_Start_bp', 'Start coordinate recorded in APK.gff; not an ancestral DNA sequence coordinate reconstruction', 'bp'),
          f('APK_End_bp', 'End coordinate recorded in APK.gff', 'bp'),
          f('APK_Strand', 'Strand token recorded in APK.gff'),
          f('Phoenix_Gene_ID', 'Original Phoenix gene ID in column 7 of APK.gff'),
          f('Phoenix_Chromosome', 'Chromosome from matched Phoenix.gff record'),
          f('Phoenix_Gene_Order', 'Original order in Phoenix.gff, before APK protein-availability exclusions', 'gene_order'),
          f('Phoenix_Start_bp', 'Start coordinate copied from Phoenix.gff without conversion', 'bp'),
          f('Phoenix_End_bp', 'End coordinate copied from Phoenix.gff without conversion', 'bp'),
          f('Phoenix_Strand', 'Strand token in Phoenix.gff'), SOURCE, STATUS], gene_rows)

    materials = [
        ('Phoenix', 'Phoenix dactylifera', 'Phoenix', 'Phoenix.len'),
        ('Areca', 'Areca catechu', 'Areca', 'Areca.len'),
        ('Cocos', 'Cocos nucifera', 'cocos', 'cocos_nucifera.len'),
        ('Oleifera', 'Elaeis oleifera', 'American', 'American.len'),
        ('Pisifera', 'Elaeis guineensis', 'pisifera', 'pisifera.len'),
        ('Dura', 'Elaeis guineensis', 'dura', 'dura.len'),
    ]
    segments, summaries = [], []
    expected = {'Phoenix': (26,18), 'Areca': (26,16), 'Cocos': (30,16),
                'Oleifera': (41,16), 'Pisifera': (42,16), 'Dura': (43,16)}
    for material, species, token, lenfile in materials:
        source = f'ancestor_{token}.txt'
        groups = defaultdict(list)
        lengths = {r[0]: int(r[2]) for _, r in read_rows(lenfile)}
        for line, r in read_rows(source):
            assert len(r) == 5
            c, s, e, col, extra = r
            s, e = int(s), int(e)
            assert col in colour_apk and 1 <= s <= e
            assert c in lengths
            assert e <= lengths[c], (source, line, e, lengths[c])
            groups[c].append((s, e, col, extra, line))
        status = 'Provisional_upstream' if material == 'Oleifera' else 'Candidate'
        b = 0
        for chrom in sorted(groups, key=natural):
            previous = None
            block_id = None
            for index, (s, e, col, extra, line) in enumerate(sorted(groups[chrom]), 1):
                if previous:
                    assert s > previous[1], (source, chrom, 'overlapping intervals')
                if previous is None or s != previous[1] + 1 or col != previous[2]:
                    b += 1
                    block_id = f'{material}_B{b:03d}'
                segments.append(['Fig. 2a', species, material, chrom, f'{material}_{chrom}_S{index:03d}',
                                 s, e, e-s+1, colour_apk[col], col, block_id, extra, line, source, status])
                previous = (s, e, col)
            if groups[chrom] and (min(r[0] for r in groups[chrom]) != 1 or max(r[1] for r in groups[chrom]) != lengths[chrom]):
                WARNINGS.append(f'{source}:{chrom}: mapping does not cover full chromosome gene-order span')
        c = len(groups)
        assert (b,c) == expected[material], (material,b,c)
        assert 10 + (b-10) - (b-c) == c
        note = 'Model-based counts; ten-copy identity/merging evidence requires review.'
        if material == 'Oleifera':
            note += ' Upstream BLAST overwritten according to HANDOFF.md; existing mapping only.'
        if material == 'Dura':
            note += ' Current files imply 33 fissions/27 fusions; HANDOFF reports figure label 31/25; reconcile figure.'
        summaries.append(['Fig. 2a', species, material, 10, b, c, b-10, b-c,
                          'B - A', 'B - C', 10+(b-10)-(b-c), source, status, note])
    table('Fig_2a_Chromosome_Segments.csv', [FIG, SPECIES, MATERIAL,
          f('Chromosome', 'Extant chromosome ID as recorded in source'),
          f('Segment_ID', 'Export identifier for an original colour interval; not an original WGDI collinear block ID'),
          f('Start_Gene_Order', 'Start of original mapped colour interval, 1-based inclusive', 'gene_order'),
          f('End_Gene_Order', 'End of original mapped colour interval, 1-based inclusive', 'gene_order'),
          f('Gene_Order_Span', 'End_Gene_Order - Start_Gene_Order + 1; not the number of supporting collinear gene pairs', 'gene_order_positions'),
          APK, COLOR, f('Continuous_Block_ID', 'Counting unit: consecutive intervals with same colour merge only if next start equals previous end + 1 on same chromosome'),
          f('WGDI_Field_5', 'Fifth source field preserved verbatim; not interpreted as duplicated ancestor identity'),
          f('Source_Line', '1-based physical line number in source file', 'line_number'), SOURCE, STATUS], segments)
    table('Fig_2a_Fission_Fusion_Counts.csv', [FIG, SPECIES, MATERIAL,
          f('Ancestral_Chromosomes_A', 'Assumed starting chromosome number in complete split-then-join scenario', 'chromosomes'),
          f('Continuous_Blocks_B', 'Number of unique Continuous_Block_ID values for this material', 'blocks'),
          f('Extant_Chromosomes_C', 'Number of distinct chromosomes represented in mapped interval file', 'chromosomes'),
          f('Model_Fissions', 'B - A; model-based operations, not independently reconstructed historical events', 'operations'),
          f('Model_Fusions', 'B - C; model-based operations, not independently reconstructed historical events', 'operations'),
          f('Fission_Formula', 'Formula used for Model_Fissions'), f('Fusion_Formula', 'Formula used for Model_Fusions'),
          f('Balance_Check_C', 'A + Model_Fissions - Model_Fusions; must equal Extant_Chromosomes_C', 'chromosomes'),
          SOURCE, STATUS, f('Notes', 'Material-specific limitations and figure reconciliation status')], summaries)

    readme = [
        ('Title', 'Source data for Fig. 2a: putative ancestral palm karyotype and chromosome composition'),
        ('Status', 'Candidate export for author review; not a submission-ready or independently audited analysis release.'),
        ('Scope', 'Fig. 2a only. Other panels of Fig. 2 are not covered. Four data CSV files plus data dictionary, README and file manifest.'),
        ('Encoding', 'UTF-8 CSV, comma delimiter, one header row, standard double-quote escaping; no merged cells or embedded spreadsheet formulas. NA means unavailable/not verified, not zero.'),
        ('Figure_link', 'Figure2-A4_1.5fold.jpg, panel a. Tables were exported from existing WGDI text inputs, not digitized from the composite image; exact visual agreement has not been verified.'),
        ('APK_definition', 'APK1-APK5 are five gene-order components anchored to Phoenix dactylifera chromosomes, not reconstructed ancestral DNA sequences.'),
        ('Ancestral_ten', 'A = 10 is the agreed starting scenario. The user reports mapping to ten ancestral chromosomes, but the ten-identity reference/mapping chain has not been independently verified. No duplicated-copy IDs were invented from the five colours.'),
        ('Literature_context', 'doi:10.1038/ng.3813 provides a published ten-chromosome ancestral-stage context; it does not validate these exact subtraction formulas or prove correspondence of ancestral nodes.'),
        ('Block_definition', 'Intervals are grouped by material/chromosome, sorted by start, and merged for counting only if directly adjacent and the same colour. All 208 source intervals remain separate blocks in these current files. Mapped colour intervals are not individual raw WGDI collinear blocks.'),
        ('Counting_model', 'Complete split-then-join model: Model_Fissions = B - A; Model_Fusions = B - C; A = 10. These are not demonstrated minimum historical event numbers or rearrangement rates. Counts are not independent measurements.'),
        ('Coordinate_units', 'Segment start/end are 1-based inclusive gene-order positions, NOT bp. APK/Phoenix bp fields are copied unchanged from WGDI-format GFF tables; original annotation coordinate convention must be checked against the release before reuse as BED.'),
        ('Material_scope', 'Six materials represent five species. Dura and Pisifera are fruit forms of Elaeis guineensis, not separate species. No replicate measurements, error bars or statistical tests are represented by these descriptive counts.'),
        ('Reference_versions', 'NA: assembly accessions and annotation releases for the six materials were not established in this bounded export. Confirm these in the manuscript/reference metadata before submission.'),
        ('Oleifera_limitation', 'HANDOFF.md reports APK_American.blastp.txt overwritten by Areca output. Oleifera rows reflect existing ancestor_American.txt only and carry Provisional_upstream status. No BLAST repair or rerun was performed here.'),
        ('Dura_limitation', 'Current source has 43 blocks, 33 model fissions and 27 model fusions. HANDOFF.md reports composite figure labels 31/25. Reconcile the chosen final mapping and figure; this export does not change the figure.'),
        ('Additional_evidence', 'Raw BLAST, protein sequences and filtered collinear-block/gene-pair evidence are not included in this compact figure-data export. These files do not constitute a complete WGDI reproducibility package.'),
        ('Source_location', 'Input basenames refer to workingdir/self_blastp_out in the owning project. No absolute server paths are included in this package.'),
        ('Journal_guidance', 'Nature Communications guidance requests numerical data underlying figures and accepts an Excel document or a zipped folder of suitable text files, with clearly linked figures/tables. CSV is the requested working format here; exact journal/editor upload requirements must be checked before submission.'),
        ('Guidance_URL', 'https://www.nature.com/ncomms/submit/how-to-submit'),
        ('Guidance_PDF', 'https://www.nature.com/documents/ncomms-submission-guide.pdf'),
        ('Guidance_scope', 'Nature Portfolio journals do not all have identical source-data upload specifications. No journal acceptance or universal Nature CSV template is claimed.'),
        ('Submission_action', 'Resolve flagged provenance/figure/ancestor/reference-version issues, review the CSVs and combine all relevant Fig. 2 source data in the journal-requested Excel/ZIP container. Add the source-data availability statement to the figure legend only when the final package exists.'),
        ('Reproduction', 'Run python3 export_karyotype_source_data.py from the owning project. The script refuses existing output directories and does not overwrite source data.'),
    ]
    table('Fig_2a_README.csv', [f('Topic', 'Documentation topic'), f('Description', 'Scope, interpretation or submission instruction')], readme)
    dictionary = list(DICTIONARY)
    table('Fig_2a_Data_Dictionary.csv', [f('File', 'CSV file described'), f('Column', 'Exact header'),
          f('Definition', 'Field definition'), f('Unit', 'Unit or not_applicable')], dictionary)

    # All biological parsing and integrity assertions above finish before any output writes.
    OUT.mkdir(parents=True, exist_ok=False)
    CHECK.mkdir(parents=True, exist_ok=False)
    manifest = []
    for name, (headers, rows) in TABLES.items():
        path = OUT / name
        with path.open('x', encoding='utf-8', newline='') as handle:
            w = csv.writer(handle)
            w.writerow(headers)
            w.writerows(rows)
        with path.open(encoding='utf-8', newline='') as handle:
            back = list(csv.reader(handle))
        assert back[0] == headers and len(back)-1 == len(rows)
        assert all(len(r) == len(headers) for r in back)
        assert '${DATA_ROOT}/' not in path.read_text()
        manifest.append([name, len(rows), path.stat().st_size, hashlib.sha256(path.read_bytes()).hexdigest(), 'Candidate'])
    with (OUT/'Fig_2a_File_Manifest.csv').open('x', encoding='utf-8', newline='') as handle:
        w = csv.writer(handle)
        w.writerow(['File','Data_Rows','Size_Bytes','SHA256','Review_Status'])
        w.writerows(manifest)
    with (CHECK/'Input_Checksums.csv').open('x', encoding='utf-8', newline='') as handle:
        w = csv.writer(handle)
        w.writerow(['Source_File','Size_Bytes','SHA256'])
        w.writerows([name,*values] for name,values in sorted(INPUTS.items()))
    report = {
        'export_consistency': 'passed', 'scientific_submission_status': 'candidate_pending_review',
        'apk_genes':len(genes), 'apk_components':len(definition), 'intervals':len(segments),
        'materials':len(summaries), 'warnings':WARNINGS,
        'checks':['unique APK gene IDs; 7294 genes', 'all APK IDs resolve to expected Phoenix chromosomes',
                  'APK orders contiguous and counts/max coordinates agree with APK.len',
                  'intervals positive, non-overlapping, within chromosome gene-order lengths',
                  'all colours resolve to APK definitions', 'counts match handoff; balance equations hold',
                  'all CSVs round-trip with stable column counts; no absolute server paths'],
        'limitations':['Oleifera upstream chain unresolved', 'Dura figure reconciliation pending',
                       'ten ancestral chromosome identity/merging chain unverified',
                       'assembly/annotation releases not established', 'full WGDI evidence not audited'],
        'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'tables':{name:len(rows) for name,(_,rows) in TABLES.items()},
    }
    with (CHECK/'Export_Validation.json').open('x', encoding='utf-8') as handle:
        json.dump(report,handle,indent=2)
        handle.write('\n')
    print(json.dumps(report,indent=2))


if __name__ == '__main__':
    main()
