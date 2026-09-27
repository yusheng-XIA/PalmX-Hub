#!/usr/bin/env python3
"""Map one chromosome-level E. oleifera assembly onto the retained APK reference."""
import argparse
import configparser
import csv
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
from collections import defaultdict
from datetime import datetime, timezone
from urllib.parse import unquote

PROJECT = Path(__file__).resolve().parents[2]
ENV = Path('${DATA_DIR2}/anaconda3/envs/WGD')
GFFREAD = Path('${DATA_DIR2}/tools/gffread/gffread')


def require(ok, message):
    if not ok:
        raise ValueError(message)


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b''):
            digest.update(chunk)
    return digest.hexdigest()


def read_fasta(path):
    records = {}
    name = None
    with Path(path).open() as handle:
        for line in handle:
            if line.startswith('>'):
                name = line[1:].split()[0]
                require(name not in records, 'Duplicate FASTA ID: ' + name)
                records[name] = []
            elif line.strip():
                require(name is not None, 'Sequence before FASTA header')
                records[name].append(line.strip())
    return {name: ''.join(seq).upper() for name, seq in records.items()}


def write_tsv(path, header, rows):
    with Path(path).open('w', newline='') as handle:
        writer = csv.writer(handle, delimiter='\t')
        if header:
            writer.writerow(header)
        writer.writerows(rows)


def annotation_models(path):
    """Parse gene -> mRNA/transcript -> CDS relationships, including shared CDS."""
    models = {}
    coding = defaultdict(list)
    with Path(path).open() as handle:
        for number, line in enumerate(handle, 1):
            if line.startswith('##FASTA'):
                break
            if line.startswith('#') or not line.strip():
                continue
            row = line.rstrip('\n').split('\t')
            require(len(row) == 9, 'Invalid GFF3 line {}'.format(number))
            chrom, _, feature, start, end, _, strand, phase, attr = row
            if feature not in ('mRNA', 'transcript', 'CDS'):
                continue
            attrs = dict(part.split('=', 1) for part in attr.split(';') if '=' in part)
            parents = [unquote(p) for p in attrs.get('Parent', '').split(',') if p]
            start, end = int(start), int(end)
            require(0 < start <= end and strand in ('+', '-'),
                    'Invalid coordinates/strand at GFF3 line {}'.format(number))
            if feature in ('mRNA', 'transcript'):
                identifier = unquote(attrs.get('ID', ''))
                require(identifier and len(parents) == 1,
                        'Transcript requires ID and one gene Parent at line {}'.format(number))
                require(identifier not in models, 'Duplicate transcript ID: ' + identifier)
                models[identifier] = (chrom, start, end, strand, parents[0])
            else:
                require(parents and phase in ('0', '1', '2'),
                        'CDS requires transcript Parent and phase at line {}'.format(number))
                for parent in parents:
                    coding[parent].append((chrom, start, end, strand))
    return models, coding


def prepare_inputs(annotation, proteins_path, fai, output, chromosomes=None, expected=16,
                   exclude_invalid_proteins=False):
    models, coding = annotation_models(annotation)
    proteins = read_fasta(proteins_path)
    lengths = {}
    for line in Path(fai).read_text().splitlines():
        row = line.split('\t')
        require(row[0] not in lengths, 'Duplicate genome sequence ID: ' + row[0])
        lengths[row[0]] = int(row[1])
    chosen_chroms = list(chromosomes) if chromosomes is not None else None
    if chosen_chroms is None:
        chosen_chroms = sorted({part[0] for parts in coding.values() for part in parts})
    require(len(set(chosen_chroms)) == len(chosen_chroms), 'Duplicate chromosome list entries')
    require(len(chosen_chroms) == expected,
            'Expected {} nuclear chromosomes, found {}. Supply --chromosomes with the '
            'exact primary chromosome IDs, one per line.'.format(expected, len(chosen_chroms)))
    require(set(chosen_chroms) <= set(lengths), 'Chromosome list does not match genome FASTA IDs')
    grouped = defaultdict(list)
    outside = []
    for transcript, parts in coding.items():
        if all(part[0] not in chosen_chroms for part in parts):
            outside.append((transcript, 'outside_selected_chromosomes'))
            continue
        require(transcript in models, 'CDS Parent has no mRNA/transcript record: ' + transcript)
        chrom, start, end, strand, gene = models[transcript]
        require(chrom in chosen_chroms and end <= lengths[chrom],
                'Transcript outside selected genome sequence: ' + transcript)
        parts.sort(key=lambda part: part[1])
        require(all(c == chrom and s == strand and start <= a <= b <= end
                    for c, a, b, s in parts), 'CDS/transcript coordinates disagree: ' + transcript)
        require(all(left[2] < right[1] for left, right in zip(parts, parts[1:])),
                'Duplicate or overlapping CDS features: ' + transcript)
        coding_length = sum(b - a + 1 for _, a, b, _ in parts)
        grouped[gene].append((transcript, coding_length))
    representatives = []
    exclusions = list(outside)
    for gene, candidates in grouped.items():
        require(len({(models[t][0], models[t][3]) for t, _ in candidates}) == 1,
                'One gene Parent spans multiple chromosomes/strands: ' + gene)
        # Longest annotated spliced CDS; transcript ID breaks ties deterministically.
        candidates.sort(key=lambda item: (-item[1], item[0]))
        transcript, cds_length = candidates[0]
        exclusions.extend((t, 'alternative_transcript') for t, _ in candidates[1:])
        require(transcript in proteins, 'Representative protein missing from extraction: ' + transcript)
        seq = proteins[transcript].replace('.', '*')
        if seq.endswith('*'):
            seq = seq[:-1]
        reason = 'empty_protein' if not seq else ('internal_stop' if '*' in seq else '')
        if reason and exclude_invalid_proteins:
            exclusions.append((transcript, reason))
            continue
        require(not reason, 'Representative protein empty/contains internal stops: ' + transcript)
        require(all('A' <= aa <= 'Z' for aa in seq), 'Unexpected protein symbol: ' + transcript)
        representatives.append((transcript, gene, cds_length, seq))
    order_chroms = {chrom: index for index, chrom in enumerate(chosen_chroms)}
    representatives.sort(key=lambda item: (order_chroms[models[item[0]][0]],
                                           models[item[0]][1], models[item[0]][2], item[0]))
    orders = defaultdict(int)
    rows, id_map, peptide_rows = [], [], []
    for index, (transcript, gene, cds_length, seq) in enumerate(representatives, 1):
        chrom, start, end, strand, _ = models[transcript]
        orders[chrom] += 1
        new_id = 'EOLg{:07d}'.format(index)
        rows.append((chrom, new_id, start, end, strand, orders[chrom], transcript))
        id_map.append((new_id, gene, transcript, chrom, start, end, strand, orders[chrom],
                       cds_length, len(seq)))
        peptide_rows.append('>{}\n{}\n'.format(new_id, seq))
    require(set(orders) == set(chosen_chroms), 'A selected chromosome has no coding representative')
    output = Path(output)
    write_tsv(output / 'Eoleifera.gff', None, rows)
    write_tsv(output / 'Eoleifera.len', None,
              [(c, lengths[c], orders[c]) for c in chosen_chroms])
    (output / 'Eoleifera.pep').write_text(''.join(peptide_rows))
    write_tsv(output / 'gene_transcript_map.tsv',
              ['WGDI_ID', 'Gene_ID', 'Transcript_ID', 'Chromosome', 'Start_bp', 'End_bp',
               'Strand', 'Gene_order', 'Annotated_CDS_bp', 'Protein_aa'], id_map)
    write_tsv(output / 'excluded_transcripts.tsv', ['Transcript_ID', 'Reason'], sorted(exclusions))
    return {'coding_representatives': len(rows), 'chromosomes': chosen_chroms,
            'alternative_transcripts': sum(reason == 'alternative_transcript' for _, reason in exclusions),
            'outside_selected_chromosomes': len(outside),
            'excluded_internal_stop': sum(reason == 'internal_stop' for _, reason in exclusions),
            'excluded_empty_protein': sum(reason == 'empty_protein' for _, reason in exclusions),
            'protein_qc_policy': 'exclude_internal_stop_and_empty' if exclude_invalid_proteins else 'strict_error',
            'partial_ORFs_without_internal_stops': 'retained; not excluded solely by Liftoff valid_ORF=False'}


def summarize(output):
    output = Path(output)
    lens = {r[0]: int(r[2]) for r in
            (line.split() for line in (output / 'Eoleifera.len').read_text().splitlines())}
    colour_apk = {r[3]: r[0] for r in
                  (line.split() for line in (output / 'ancestor_out.txt').read_text().splitlines())}
    grouped = defaultdict(list)
    for line in (output / 'ancestor_Eoleifera.txt').read_text().splitlines():
        if not line.strip():
            continue
        chrom, a, b, colour, classification = line.split()
        a, b = int(a), int(b)
        require(chrom in lens and 1 <= a <= b <= lens[chrom], 'Mapped interval out of bounds')
        require(colour in colour_apk and classification == '1', 'Unexpected APK colour/class')
        grouped[chrom].append((a, b, colour))
    require(set(grouped) == set(lens), 'Mapping omitted a selected chromosome')
    segments = []
    for chrom in lens:
        end = 0
        for a, b, colour in sorted(grouped[chrom]):
            require(a == end + 1, 'Overlapping or missing mapped gene orders: ' + chrom)
            if segments and segments[-1][0] == chrom and segments[-1][3] == colour:
                segments[-1][2] = b
            else:
                segments.append([chrom, a, b, colour])
            end = b
        require(end == lens[chrom], 'Mapping does not reach chromosome end: ' + chrom)
    seeds = defaultdict(set)
    with (output / 'APK_Eoleifera.correspondence.csv').open() as handle:
        blocks = list(csv.DictReader(handle))
    for row in blocks:
        require(float(row['pvalue']) <= .2 and int(row['length']) >= 5,
                'Correspondence filter differs from declared parameters')
        # WGDI 0.75 -km uses length > block_length, whereas -c uses >= block_length.
        if int(row['length']) > 5:
            seeds[(row['chr1'], row['chr2'])].update(int(x) for x in row['block1'].split('_'))
    table = []
    for chrom, a, b, colour in segments:
        component = colour_apk[colour]
        n = sum(a <= position <= b for position in seeds[(chrom, component)])
        table.append((chrom, a, b, b - a + 1, 'APK' + component, colour, n,
                      round(n / (b - a + 1), 6)))
    write_tsv(output / 'karyotype_segments.tsv',
              ['Chromosome', 'Start_gene_order', 'End_gene_order', 'Gene_order_span',
               'APK_component', 'Colour', 'Matching_seed_positions', 'Matching_seed_fraction'], table)
    # Native -k draws height=end-start. Convert inclusive orders to boundary coordinates
    # so adjacent blocks meet and a one-gene interval has nonzero height.
    write_tsv(output / 'ancestor_Eoleifera.plot.txt', None,
              [(chrom, a - 1, b, colour, 1) for chrom, a, b, colour in segments])
    B, C, A = len(segments), len(lens), 10
    require(B >= max(A, C), 'Complete split-then-join model is not applicable to these counts')
    counts = {'A': A, 'B': B, 'C': C, 'Model_fissions': B - A, 'Model_fusions': B - C}
    write_tsv(output / 'model_counts.tsv', list(counts), [list(counts.values())])
    return dict(counts, coordinate_unit='1-based inclusive gene order',
                interpretation='Model operations, not independently inferred historical events',
                support_note='Seed support precedes WGDI BLAST-based recolouring and gap extension; '
                             'not a posterior probability or the complete assignment evidence',
                scientific_status='Candidate requiring interval and rendered-figure review')


def wgdi_config(threads):
    shared = dict(gff1='Eoleifera.gff', gff2='APK.gff', lens1='Eoleifera.len',
                  lens2='APK.len', blast='APK_Eoleifera.blastp.tsv', score='100', evalue='1e-5')
    config = configparser.ConfigParser()
    config['collinearity'] = dict(shared, blast_reverse='false', comparison='genomes',
        multiple='1', process=str(threads), grading='50,30,25', mg='25,25', pvalue='1',
        repeat_number='10', position='order', savefile='APK_Eoleifera.collinearity')
    config['blockinfo'] = dict(shared, collinearity='APK_Eoleifera.collinearity',
        repeat_number='20', position='order', ks='none', ks_col='ks_NG86',
        savefile='APK_Eoleifera.block_info.csv')
    config['correspondence'] = dict(blockinfo='APK_Eoleifera.block_info.csv',
        lens1='Eoleifera.len', lens2='APK.len', tandem='true', tandem_length='200',
        pvalue='0.2', block_length='5', multiple='1', homo='0,1',
        savefile='APK_Eoleifera.correspondence.csv')
    config['karyotype_mapping'] = dict(blast='APK_Eoleifera.blastp.tsv', blast_reverse='false',
        gff1='Eoleifera.gff', gff2='APK.gff', score='100', evalue='1e-5', repeat_number='5',
        ancestor_top='ancestor_out.txt', the_other_lens='Eoleifera.len',
        blockinfo='APK_Eoleifera.correspondence.csv', blockinfo_reverse='false',
        block_length='5', limit_length='5', the_other_ancestor_file='ancestor_Eoleifera.txt')
    return config


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--genome', type=Path, required=True, help='Uncompressed genome FASTA')
    parser.add_argument('--annotation', type=Path, required=True, help='Matching gene/mRNA/CDS GFF3')
    parser.add_argument('--output', type=Path, required=True, help='New run directory; must not exist')
    parser.add_argument('--chromosomes', type=Path, help='Primary nuclear chromosome IDs, in display order')
    parser.add_argument('--exclude-invalid-proteins', action='store_true',
                        help='Exclude and record internal-stop/empty representative proteins; retain partial ORFs')
    parser.add_argument('--font-cache', type=Path,
                        help='Optional compatible Matplotlib fontlist-v330.json copied into the run')
    parser.add_argument('--reference', type=Path, default=PROJECT / 'workingdir/self_blastp_out')
    args = parser.parse_args()
    require(os.environ.get('SLURM_JOB_ID'), 'Run the full workflow in a SLURM allocation')
    threads = int(os.environ.get('SLURM_CPUS_PER_TASK', '1'))
    require(1 <= threads <= 8, 'This workflow is bounded to 1-8 CPUs')
    genome, annotation, reference = (p.resolve(strict=True) for p in
                                     (args.genome, args.annotation, args.reference))
    require(genome.suffix != '.gz', 'Use an uncompressed FASTA')
    paths = [genome, annotation] + [reference / name for name in
                                   ('APK.gff', 'APK.len', 'APK.pep', 'ancestor_out.txt')]
    require(all(p.is_file() for p in paths), 'Missing required input/reference file')
    require(all(p.is_file() for p in (GFFREAD, ENV / 'bin/python', ENV / 'bin/wgdi',
                                     ENV / 'bin/blastp', ENV / 'bin/makeblastdb')), 'Missing executable')
    chromosomes = None
    if args.chromosomes:
        chromosome_file = args.chromosomes.resolve(strict=True)
        paths.append(chromosome_file)
        chromosomes = [x.strip() for x in chromosome_file.read_text().splitlines()
                       if x.strip() and not x.startswith('#')]
    font_cache = None
    if args.font_cache:
        font_cache = args.font_cache.resolve(strict=True)
        require(json.loads(font_cache.read_text()).get('_version') == 330,
                'Expected Matplotlib font cache schema 330')
        paths.append(font_cache)
    output = args.output.absolute()
    require(not output.exists() and not output.is_symlink(), 'Output already exists: ' + str(output))
    require(output.parent == output.parent.resolve(), 'Output parent contains a symlink')
    output.parent.mkdir(parents=True, exist_ok=True)
    output.mkdir()
    state = {'status': 'running', 'started_utc': datetime.now(timezone.utc).isoformat(),
             'slurm_job_id': os.environ['SLURM_JOB_ID'], 'cpus': threads,
             'A': 10, 'reference_scope': 'Existing five APK components; assumed doubling to ACPK=10'}
    env = dict(os.environ, PATH=str(ENV / 'bin') + ':' + os.environ.get('PATH', ''),
               LD_LIBRARY_PATH=str(ENV / 'lib') + (':' + os.environ['LD_LIBRARY_PATH'] if os.environ.get('LD_LIBRARY_PATH') else ''),
               MPLBACKEND='Agg', MPLCONFIGDIR=str(output / 'mplconfig'),
               TMPDIR=str(output / 'tmp'), OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1',
               MKL_NUM_THREADS='1', PYTHONDONTWRITEBYTECODE='1', PYTHONNOUSERSITE='1',
               PYTHONPATH='')
    (output / 'tmp').mkdir()
    (output / 'mplconfig').mkdir()
    if font_cache:
        shutil.copy2(font_cache, output / 'mplconfig/fontlist-v330.json')
    commands = output / 'commands.jsonl'

    def run(command, log, timeout=None):
        command = [str(x) for x in command]
        with commands.open('a') as handle:
            handle.write(json.dumps(dict(argv=command, cwd=str(output))) + '\n')
        with (output / log).open('w') as handle:
            subprocess.run(command, cwd=output, env=env, stdout=handle,
                           stderr=subprocess.STDOUT, check=True, timeout=timeout)

    (output / 'run_status.json').write_text(json.dumps(state, indent=2) + '\n')
    try:
        state['inputs'] = [{'path': str(p), 'sha256': sha256(p)} for p in paths]
        shutil.copy2(__file__, output / 'run_pipeline.py')
        state['script_sha256'] = sha256(__file__)
        (output / 'genome.fa').symlink_to(genome)
        for p in paths[2:6]:
            shutil.copy2(p, output / p.name)
        # Full CLI health was checked separately; querying metadata avoids importing every
        # unused WGDI module solely to print its version at the start of each run.
        run([ENV / 'bin/python', '-c', 'from importlib.metadata import version; print(version("wgdi"))'],
            'wgdi.version.txt', timeout=120)
        require('0.75' in (output / 'wgdi.version.txt').read_text().splitlines(), 'Expected WGDI 0.75')
        for executable, flag, label in [(GFFREAD, '--version', 'gffread'),
                                       (ENV / 'bin/blastp', '-version', 'blastp')]:
            run([executable, flag], label + '.version.txt', timeout=120)
        run([ENV / 'bin/python', '-c',
             "import numpy,pandas,Bio,matplotlib,json; from importlib.metadata import version; "
             "print(json.dumps({'wgdi':version('wgdi'),'numpy':numpy.__version__,"
             "'pandas':pandas.__version__,'biopython':Bio.__version__,"
             "'matplotlib':matplotlib.__version__}))"], 'packages.txt', timeout=180)
        source = ENV / 'lib/python3.8/site-packages/wgdi'
        shutil.copytree(source, output / 'code/wgdi', ignore=shutil.ignore_patterns('__pycache__'))
        env['PYTHONPATH'] = str(output / 'code')
        run([GFFREAD, annotation, '-g', 'genome.fa', '-y', 'transcripts.pep.fa'], 'gffread.log')
        require((output / 'genome.fa.fai').is_file(), 'gffread did not create the run-local FASTA index')
        state['preparation'] = prepare_inputs(annotation, output / 'transcripts.pep.fa',
                                              output / 'genome.fa.fai', output, chromosomes,
                                              exclude_invalid_proteins=args.exclude_invalid_proteins)
        config = wgdi_config(threads)
        with (output / 'APK_Eoleifera.conf').open('w') as handle:
            config.write(handle)
        run([ENV / 'bin/makeblastdb', '-in', 'APK.pep', '-dbtype', 'prot',
             '-parse_seqids', '-out', 'APK.db'], 'makeblastdb.log')
        run([ENV / 'bin/blastp', '-query', 'Eoleifera.pep', '-db', 'APK.db',
             '-evalue', '1e-5', '-max_target_seqs', '20', '-outfmt', '6',
             '-num_threads', threads, '-out', 'APK_Eoleifera.blastp.tsv'], 'blastp.log')
        for step in ('icl', 'bi', 'c', 'km'):
            run([ENV / 'bin/wgdi', '-' + step, 'APK_Eoleifera.conf'], 'wgdi_' + step + '.log')
        state['summary'] = summarize(output)
        # Native WGDI exports are analysis figures; composite layout/visual QA is a later step.
        for extension in ('pdf', 'svg', 'png'):
            plot = configparser.ConfigParser()
            plot['karyotype'] = dict(ancestor='ancestor_Eoleifera.plot.txt', width='0.5',
                                    figsize='10,6.18', savefig='Eoleifera.karyotype.' + extension)
            filename = 'karyotype_' + extension + '.conf'
            with (output / filename).open('w') as handle:
                plot.write(handle)
            run([ENV / 'bin/wgdi', '-k', filename], 'plot_' + extension + '.log')
            require((output / ('Eoleifera.karyotype.' + extension)).stat().st_size > 0,
                    'Empty figure output')
        require(all(sha256(p) == record['sha256'] for p, record in zip(paths, state['inputs'])),
                'An input changed during execution')
        state['status'] = 'completed_candidate_not_scientific_signoff'
    except BaseException as error:
        state.update(status='failed', error=repr(error))
        raise
    finally:
        state['ended_utc'] = datetime.now(timezone.utc).isoformat()
        (output / 'run_status.json').write_text(json.dumps(state, indent=2) + '\n')
        write_tsv(output / 'output_checksums.tsv', ['File', 'SHA256'],
                  [(str(p.relative_to(output)), sha256(p)) for p in sorted(output.rglob('*'))
                   if p.is_file() and not p.is_symlink() and 'tmp' not in p.relative_to(output).parts
                   and 'mplconfig' not in p.relative_to(output).parts
                   and p.name != 'output_checksums.tsv'])


if __name__ == '__main__':
    main()
