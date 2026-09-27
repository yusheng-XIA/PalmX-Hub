#!/usr/bin/env python3
"""Prepare, validate, and package EG11-vs-FL inputs for local SyntenyViz use."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import math
import os
import re
import shutil
import statistics
import tarfile
from collections import Counter, defaultdict
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path


CHROMS = [f"chr{i:02d}" for i in range(1, 17)]
CHR_RANK = {chrom: index for index, chrom in enumerate(CHROMS)}
N_RUN = re.compile(b"[Nn]+")


@dataclass(frozen=True)
class Chromosome:
    source_id: str
    display_id: str
    length: int


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--parameters", required=True, type=Path)
    parser.add_argument("--paf", required=True, type=Path)
    parser.add_argument("--project-root", required=True, type=Path)
    parser.add_argument("--window-size", type=int, default=200_000)
    parser.add_argument("--display-gap-min", type=int, default=100)
    parser.add_argument("--min-mapq", type=int, default=20)
    parser.add_argument("--min-block", type=int, default=50_000)
    return parser.parse_args()


def attrs(text: str) -> dict[str, str]:
    result: dict[str, str] = {}
    for item in text.split(";"):
        if "=" in item:
            key, value = item.split("=", 1)
            result[key] = value
    return result


def write_tsv(path: Path, header: list[str] | None, rows) -> None:
    if path.exists():
        raise FileExistsError(f"Refusing to overwrite: {path}")
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        if header is not None:
            writer.writerow(header)
        writer.writerows(rows)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while True:
            chunk = handle.read(1024 * 1024)
            if not chunk:
                break
            digest.update(chunk)
    return digest.hexdigest()


def file_record(kind: str, sample: str, path: Path, digest: str) -> list[object]:
    stat = path.stat()
    return [
        kind,
        sample,
        str(path),
        str(path.resolve()),
        stat.st_size,
        datetime.fromtimestamp(stat.st_mtime, timezone.utc).isoformat(),
        digest,
    ]


def load_chromosome_map(path: Path) -> tuple[list[Chromosome], dict[str, str]]:
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    chromosomes = [Chromosome(row["Source_ID"], row["Standard_ID"], int(row["Length_bp"])) for row in rows]
    if [item.display_id for item in chromosomes] != CHROMS:
        raise ValueError(f"Chromosome map must contain chr01-chr16 in natural order: {path}")
    if len({item.source_id for item in chromosomes}) != 16:
        raise ValueError(f"Duplicate source chromosome in map: {path}")
    return chromosomes, {item.source_id: item.display_id for item in chromosomes}


def update_gaps(sequence: bytes, offset: int, open_start: int | None, gaps: list[tuple[int, int]]) -> int | None:
    if not sequence:
        return open_start
    start_search = 0
    leading = re.match(b"[Nn]+", sequence)
    if open_start is not None:
        if leading is None:
            gaps.append((open_start, offset))
            open_start = None
        elif leading.end() < len(sequence):
            gaps.append((open_start, offset + leading.end()))
            open_start = None
            start_search = leading.end()
        else:
            return open_start
    for match in N_RUN.finditer(sequence, start_search):
        start = offset + match.start()
        end = offset + match.end()
        if match.end() == len(sequence):
            open_start = start
        else:
            gaps.append((start, end))
    return open_start


def scan_fasta(path: Path, chromosomes: list[Chromosome]):
    expected = {item.source_id: item for item in chromosomes}
    observed: dict[str, int] = {}
    gaps_by_display: dict[str, list[tuple[int, int]]] = {chrom: [] for chrom in CHROMS}
    n_by_display: Counter[str] = Counter()
    invalid_by_display: Counter[str] = Counter()
    digest = hashlib.sha256()
    name: str | None = None
    display: str | None = None
    length = 0
    open_start: int | None = None

    def finish() -> None:
        nonlocal name, display, length, open_start
        if name is None or display is None:
            return
        if open_start is not None:
            gaps_by_display[display].append((open_start, length))
            open_start = None
        if name in observed:
            raise ValueError(f"Duplicate FASTA sequence name: {name}")
        observed[name] = length

    with path.open("rb") as handle:
        for raw in handle:
            digest.update(raw)
            if raw.startswith(b">"):
                finish()
                name = raw[1:].split(maxsplit=1)[0].decode("ascii")
                if name not in expected:
                    raise ValueError(f"FASTA contains sequence absent from chromosome map: {name}")
                display = expected[name].display_id
                length = 0
                open_start = None
                continue
            if name is None or display is None:
                raise ValueError(f"Sequence data precedes FASTA header: {path}")
            sequence = raw.strip()
            n_by_display[display] += sequence.count(b"N") + sequence.count(b"n")
            invalid_by_display[display] += len(re.findall(b"[^ACGTRYSWKMBDHVacgtryswkmbdhvNn]", sequence))
            open_start = update_gaps(sequence, length, open_start, gaps_by_display[display])
            length += len(sequence)
    finish()
    if set(observed) != set(expected):
        raise ValueError(f"FASTA/map sequence mismatch: missing={sorted(set(expected)-set(observed))}")
    for source_id, observed_length in observed.items():
        if observed_length != expected[source_id].length:
            raise ValueError(f"FASTA length mismatch {source_id}: {observed_length} != {expected[source_id].length}")
    if sum(invalid_by_display.values()):
        raise ValueError(f"FASTA contains {sum(invalid_by_display.values())} invalid sequence characters: {path}")
    return observed, gaps_by_display, n_by_display, invalid_by_display, digest.hexdigest()


def n50(lengths: list[int]) -> int:
    threshold = sum(lengths) / 2
    cumulative = 0
    for length in sorted(lengths, reverse=True):
        cumulative += length
        if cumulative >= threshold:
            return length
    raise RuntimeError("Cannot calculate N50")


def parse_gff(
    path: Path,
    sample: str,
    chromosomes: list[Chromosome],
    source_map: dict[str, str],
    map_mode: str,
    window_size: int,
):
    lengths = {item.display_id: item.length for item in chromosomes}
    windows = {chrom: [0] * math.ceil(lengths[chrom] / window_size) for chrom in CHROMS}
    feature_counts: Counter[str] = Counter()
    all_feature_counts: Counter[str] = Counter()
    excluded_features: Counter[str] = Counter()
    gene_ids: set[str] = set()
    gff_map: dict[str, str] = dict(source_map) if map_mode == "GFF_seqid_equals_FASTA_seqid" else {}
    gff_name_rows: list[list[object]] = []
    digest = hashlib.sha256()
    invalid_coordinates = 0
    mapped_gene_count = 0

    if map_mode == "GFF_region_chromosome_attribute":
        display_to_source: dict[str, str] = {}
        with path.open("rb") as handle:
            for line_number, raw in enumerate(handle, 1):
                if raw.startswith(b"#") or not raw.strip():
                    continue
                fields = raw.decode("utf-8").rstrip("\n\r").split("\t")
                if len(fields) != 9:
                    raise ValueError(f"Malformed GFF line {line_number}: {path}")
                if fields[2] != "region":
                    continue
                values = attrs(fields[8])
                chromosome = values.get("chromosome", "")
                if not (chromosome.isdigit() and 1 <= int(chromosome) <= 16):
                    continue
                seqid = fields[0]
                display = f"chr{int(chromosome):02d}"
                if int(fields[3]) != 1 or int(fields[4]) != lengths[display]:
                    raise ValueError(f"EG11 GFF region coordinate/length mismatch: {seqid} -> {display}")
                if seqid in gff_map or display in display_to_source:
                    raise ValueError(f"EG11 GFF chromosome mapping is not one-to-one: {seqid} -> {display}")
                gff_map[seqid] = display
                display_to_source[display] = seqid
        if set(gff_map.values()) != set(CHROMS) or len(gff_map) != 16:
            raise ValueError("EG11 GFF did not map exactly 16 source seqids one-to-one to chromosomes")
        for seqid, display in sorted(gff_map.items(), key=lambda item: CHR_RANK[item[1]]):
            gff_name_rows.append([
                sample, seqid, display, lengths[display], "Included",
                "GFF RefSeq seqid mapped by chromosome attribute and exact region length", "GFF",
            ])

    with path.open("rb") as handle:
        for line_number, raw in enumerate(handle, 1):
            digest.update(raw)
            if raw.startswith(b"#") or not raw.strip():
                continue
            line = raw.decode("utf-8").rstrip("\n\r")
            fields = line.split("\t")
            if len(fields) != 9:
                raise ValueError(f"Malformed GFF line {line_number}: {path}")
            seqid, feature = fields[0], fields[2]
            values = attrs(fields[8])
            all_feature_counts[feature] += 1
            if map_mode == "GFF_region_chromosome_attribute" and feature == "region":
                if seqid in gff_map:
                    feature_counts[feature] += 1
                else:
                    excluded_features[seqid] += 1
                continue
            display = gff_map.get(seqid)
            if display is None:
                excluded_features[seqid] += 1
                continue
            start = int(fields[3])
            end = int(fields[4])
            if start < 1 or end < start or end > lengths[display]:
                invalid_coordinates += 1
                continue
            feature_counts[feature] += 1
            if feature != "gene":
                continue
            gene_id = values.get("ID", "")
            if not gene_id:
                raise ValueError(f"Gene without ID at {path}:{line_number}")
            if gene_id in gene_ids:
                raise ValueError(f"Duplicate gene ID: {gene_id}")
            gene_ids.add(gene_id)
            midpoint0 = (start + end - 2) // 2
            windows[display][midpoint0 // window_size] += 1
            mapped_gene_count += 1
    if invalid_coordinates:
        raise ValueError(f"GFF has {invalid_coordinates} invalid mapped coordinates: {path}")
    if mapped_gene_count == 0:
        raise ValueError(f"No gene features on selected chromosomes: {path}")
    return windows, feature_counts, all_feature_counts, excluded_features, gff_name_rows, digest.hexdigest()


def write_sample_tracks(
    stage: Path,
    sample: str,
    chromosomes: list[Chromosome],
    gaps_by_display: dict[str, list[tuple[int, int]]],
    n_by_display: Counter[str],
    invalid_by_display: Counter[str],
    windows: dict[str, list[int]],
    feature_counts: Counter[str],
    window_size: int,
    display_gap_min: int,
):
    sample_dir = stage / sample
    sample_dir.mkdir()
    lengths = {item.display_id: item.length for item in chromosomes}
    write_tsv(sample_dir / f"{sample}.genome_sizes.tsv", None, [[chrom, lengths[chrom]] for chrom in CHROMS])

    density_rows = []
    density_summary = []
    for chrom in CHROMS:
        counts = windows[chrom]
        for index, count in enumerate(counts):
            start = index * window_size
            end = min(start + window_size, lengths[chrom])
            density_rows.append([chrom, start, end, count])
        density_summary.append([
            chrom,
            sum(counts),
            sum(value > 0 for value in counts),
            max(counts),
            f"{statistics.mean(counts):.6f}",
            f"{statistics.median(counts):.6f}",
            len(counts),
        ])
    if sum(int(row[3]) for row in density_rows) != feature_counts["gene"]:
        raise RuntimeError(f"Gene-density conservation failed for {sample}")
    write_tsv(sample_dir / f"{sample}.gene_density.200kb.bedgraph", None, density_rows)
    write_tsv(
        sample_dir / f"{sample}.gene_density.chromosome_summary.tsv",
        ["Chromosome", "Gene_Count", "Nonzero_Window_Count", "Maximum_Window_Gene_Count", "Mean_Window_Gene_Count", "Median_Window_Gene_Count", "Window_Count"],
        density_summary,
    )

    complete_rows = []
    display_rows = []
    gap_summary = []
    for chrom in CHROMS:
        complete = gaps_by_display[chrom]
        shown = [(start, end) for start, end in complete if end - start >= display_gap_min]
        complete_rows.extend([[chrom, start, end] for start, end in complete])
        display_rows.extend([[chrom, start, end] for start, end in shown])
        gap_summary.append([
            chrom,
            len(complete),
            sum(end - start for start, end in complete),
            len(shown),
            sum(end - start for start, end in shown),
            n_by_display[chrom],
            invalid_by_display[chrom],
        ])
    write_tsv(sample_dir / f"{sample}.N_gaps.complete.bed", None, complete_rows)
    write_tsv(sample_dir / f"{sample}.N_gaps.display_min100bp.bed", None, display_rows)
    write_tsv(
        sample_dir / f"{sample}.N_gaps.chromosome_summary.tsv",
        ["Chromosome", "Complete_Gap_Count", "Complete_Gap_bp", "Display_Gap_Count", "Display_Gap_bp", "FASTA_N_bp", "Invalid_FASTA_Character_Count"],
        gap_summary,
    )

    end_rows = []
    for chrom in CHROMS:
        end_rows.append([chrom, 0, 1, ".", "chromosome_end", "left"])
        end_rows.append([chrom, lengths[chrom] - 1, lengths[chrom], ".", "chromosome_end", "right"])
    write_tsv(sample_dir / f"{sample}.chromosome_ends.tsv", ["#chr", "start", "end", "score", "type", "label"], end_rows)
    return density_summary, gap_summary


def process_paf(
    source: Path,
    stage: Path,
    target_chromosomes: list[Chromosome],
    query_chromosomes: list[Chromosome],
    min_mapq: int,
    min_block: int,
):
    output_dir = stage / "alignment"
    output_dir.mkdir()
    target_map = {item.source_id: item.display_id for item in target_chromosomes}
    query_map = {item.source_id: item.display_id for item in query_chromosomes}
    target_lengths = {item.display_id: item.length for item in target_chromosomes}
    query_lengths = {item.display_id: item.length for item in query_chromosomes}
    digest = hashlib.sha256()
    records: list[tuple[tuple[int, int, int, int, str], list[str]]] = []
    source_type: Counter[str] = Counter()

    with source.open("rb") as handle:
        for line_number, raw in enumerate(handle, 1):
            digest.update(raw)
            if not raw.strip():
                raise ValueError(f"Unexpected blank PAF line: {line_number}")
            fields = raw.decode("utf-8").rstrip("\n\r").split("\t")
            if len(fields) < 12:
                raise ValueError(f"PAF line has fewer than 12 fields: {line_number}")
            query_source, target_source = fields[0], fields[5]
            if query_source not in query_map or target_source not in target_map:
                raise ValueError(f"PAF direction/name mismatch at line {line_number}: {query_source} -> {target_source}")
            query, target = query_map[query_source], target_map[target_source]
            qlen, qstart, qend = int(fields[1]), int(fields[2]), int(fields[3])
            tlen, tstart, tend = int(fields[6]), int(fields[7]), int(fields[8])
            if qlen != query_lengths[query] or tlen != target_lengths[target]:
                raise ValueError(f"PAF length mismatch at line {line_number}")
            if not (0 <= qstart < qend <= qlen and 0 <= tstart < tend <= tlen):
                raise ValueError(f"PAF coordinate error at line {line_number}")
            if fields[4] not in {"+", "-"}:
                raise ValueError(f"Invalid PAF strand at line {line_number}")
            fields[0] = query
            fields[5] = target
            tag = "NO_tp"
            for value in fields[12:]:
                if value.startswith("tp:A:"):
                    tag = value
                    break
            source_type[tag] += 1
            key = (CHR_RANK[target], tstart, CHR_RANK[query], qstart, fields[4])
            records.append((key, fields))
    records.sort(key=lambda item: item[0])

    raw_path = output_dir / "EG11_vs_FL.raw.paf.gz"
    filtered_path = output_dir / "EG11_vs_FL.filtered.mapq20.block50kb.paf"
    pair_stats: dict[tuple[str, str], Counter[str]] = defaultdict(Counter)
    total_stats: Counter[str] = Counter()
    with raw_path.open("wb") as raw_handle, gzip.GzipFile(filename="", mode="wb", fileobj=raw_handle, mtime=0) as zipped, filtered_path.open("w") as filtered:
        for _, fields in records:
            line = "\t".join(fields) + "\n"
            zipped.write(line.encode("utf-8"))
            total_stats["Raw_Alignment_Count"] += 1
            total_stats[f"Raw_{fields[4]}_Count"] += 1
            block = int(fields[10])
            mapq = int(fields[11])
            if mapq < min_mapq or block < min_block:
                continue
            filtered.write(line)
            target, query = fields[5], fields[0]
            pair_stats[(target, query)]["Alignment_Count"] += 1
            pair_stats[(target, query)]["Block_bp"] += block
            pair_stats[(target, query)][f"{fields[4]}_Count"] += 1
            pair_stats[(target, query)]["Matching_bp"] += int(fields[9])
            total_stats["Filtered_Alignment_Count"] += 1
            total_stats[f"Filtered_{fields[4]}_Count"] += 1
            total_stats["Filtered_Block_bp"] += block
            total_stats["Filtered_Matching_bp"] += int(fields[9])
    if total_stats["Filtered_Alignment_Count"] == 0:
        raise RuntimeError("PAF filtering removed all alignments")

    pair_rows = []
    for target in CHROMS:
        for query in CHROMS:
            values = pair_stats.get((target, query), Counter())
            if values["Alignment_Count"]:
                identity = values["Matching_bp"] / values["Block_bp"]
                pair_rows.append([target, query, values["Alignment_Count"], values["Block_bp"], values["+_Count"], values["-_Count"], f"{identity:.6f}"])
    write_tsv(
        output_dir / "EG11_vs_FL.chromosome_pair_summary.tsv",
        ["EG11_Chromosome", "FL_Chromosome", "Alignment_Count", "Alignment_Block_bp", "Forward_Block_Count", "Reverse_Block_Count", "Weighted_Match_Fraction"],
        pair_rows,
    )

    same_total = sum(values["Block_bp"] for (target, query), values in pair_stats.items() if target == query)
    all_total = sum(values["Block_bp"] for values in pair_stats.values())
    target_rankings: dict[str, list[tuple[int, str]]] = {}
    query_rankings: dict[str, list[tuple[int, str]]] = {}
    for target in CHROMS:
        target_rankings[target] = sorted(
            [(values["Block_bp"], query) for (candidate_target, query), values in pair_stats.items() if candidate_target == target],
            reverse=True,
        )
    for query in CHROMS:
        query_rankings[query] = sorted(
            [(values["Block_bp"], target) for (target, candidate_query), values in pair_stats.items() if candidate_query == query],
            reverse=True,
        )

    primary_query = {target: ranked[0][1] for target, ranked in target_rankings.items() if ranked}
    primary_target = {query: ranked[0][1] for query, ranked in query_rankings.items() if ranked}
    correspondence_rows = []
    correspondence_status: Counter[str] = Counter()

    def add_correspondence(direction: str, source: str, ranked: list[tuple[int, str]], reciprocal_map: dict[str, str]) -> None:
        if not ranked:
            correspondence_rows.append([direction, source, "NA", 0, 0, "0.000000", "NA", "0.000000", "No", "No", "FAIL", "No filtered alignment"])
            correspondence_status["FAIL"] += 1
            return
        main_bp, main_match = ranked[0]
        total_bp = sum(value for value, _ in ranked)
        second_bp, second_match = ranked[1] if len(ranked) > 1 else (0, "NA")
        reciprocal = reciprocal_map.get(main_match) == source
        same_name = main_match == source
        ambiguous = second_bp >= main_bp * 0.25
        status = "PASS" if reciprocal and same_name and not ambiguous else "WARN"
        notes = []
        if not same_name:
            notes.append("primary match is not same-name")
        if not reciprocal:
            notes.append("primary match is not reciprocal")
        if ambiguous:
            notes.append("secondary block length is >=25% of primary")
        if not notes:
            notes.append("same-name reciprocal primary match")
        correspondence_status[status] += 1
        correspondence_rows.append([
            direction, source, main_match, main_bp, total_bp, f"{main_bp/total_bp:.6f}",
            second_match, f"{second_bp/total_bp:.6f}", "Yes" if reciprocal else "No",
            "Yes" if same_name else "No", status, "; ".join(notes),
        ])

    for target in CHROMS:
        add_correspondence("EG11_to_FL", target, target_rankings[target], primary_target)
    for query in CHROMS:
        add_correspondence("FL_to_EG11", query, query_rankings[query], primary_query)
    write_tsv(
        output_dir / "EG11_vs_FL.primary_correspondence.tsv",
        ["Direction", "Source_Chromosome", "Primary_Match", "Primary_Block_bp", "All_Block_bp", "Primary_Block_Fraction", "Secondary_Match", "Secondary_Block_Fraction", "Reciprocal_Primary", "Same_Name_Primary", "Status", "Note"],
        correspondence_rows,
    )
    total_stats["Same_Name_Block_bp"] = same_total
    total_stats["All_Pair_Block_bp"] = all_total
    total_stats["Same_Name_Block_Fraction"] = f"{same_total / all_total:.6f}" if all_total else "0.000000"
    total_stats["Weighted_Match_Fraction"] = f"{total_stats['Filtered_Matching_bp'] / total_stats['Filtered_Block_bp']:.6f}"
    for status in ("PASS", "WARN", "FAIL"):
        total_stats[f"Correspondence_{status}_Count"] = correspondence_status[status]
    return digest.hexdigest(), total_stats, source_type


def validate_density(path: Path, lengths: dict[str, int], window_size: int) -> int:
    previous: dict[str, int] = {chrom: 0 for chrom in CHROMS}
    counts: Counter[str] = Counter()
    window_counts: Counter[str] = Counter()
    last_rank = -1
    rows = 0
    with path.open() as handle:
        for line in handle:
            if not line.strip():
                raise ValueError(f"Unexpected blank density line: {path}")
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 4:
                raise ValueError(f"Density file is not four-column TSV: {path}")
            chrom, start, end, count = fields[0], int(fields[1]), int(fields[2]), int(fields[3])
            if chrom not in lengths or count < 0 or start != previous[chrom] or not (0 <= start < end <= lengths[chrom]):
                raise ValueError(f"Density coordinate/continuity error: {path}: {line.rstrip()}")
            rank = CHR_RANK[chrom]
            if rank < last_rank:
                raise ValueError(f"Density chromosome sort error: {path}")
            if end - start != window_size and end != lengths[chrom]:
                raise ValueError(f"Non-terminal density window has wrong width: {path}")
            previous[chrom] = end
            counts[chrom] += count
            window_counts[chrom] += 1
            last_rank = rank
            rows += 1
    if any(previous[chrom] != lengths[chrom] for chrom in CHROMS):
        raise ValueError(f"Density does not cover complete chromosomes: {path}")
    if any(window_counts[chrom] != math.ceil(lengths[chrom] / window_size) for chrom in CHROMS):
        raise ValueError(f"Density window-count mismatch: {path}")
    return rows


def validate_bed(path: Path, lengths: dict[str, int], minimum: int = 1) -> int:
    previous_end: dict[str, int | None] = {chrom: None for chrom in CHROMS}
    last_rank = -1
    rows = 0
    with path.open() as handle:
        for line in handle:
            if not line.strip():
                raise ValueError(f"Unexpected blank BED line: {path}")
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 3:
                raise ValueError(f"BED is not three-column TSV: {path}")
            chrom, start, end = fields[0], int(fields[1]), int(fields[2])
            if chrom not in lengths or not (0 <= start < end <= lengths[chrom]) or end - start < minimum:
                raise ValueError(f"BED coordinate error: {path}: {line.rstrip()}")
            rank = CHR_RANK[chrom]
            if rank < last_rank or (previous_end[chrom] is not None and start <= previous_end[chrom]):
                raise ValueError(f"BED sort error: {path}")
            previous_end[chrom] = end
            last_rank = rank
            rows += 1
    return rows


def read_bed_set(path: Path) -> set[tuple[str, int, int]]:
    rows: set[tuple[str, int, int]] = set()
    with path.open() as handle:
        for line in handle:
            if line.strip():
                chrom, start, end = line.rstrip("\n").split("\t")
                rows.add((chrom, int(start), int(end)))
    return rows


def validate_genome_sizes(path: Path, chromosomes: list[Chromosome]) -> None:
    expected = [(item.display_id, item.length) for item in chromosomes]
    observed: list[tuple[str, int]] = []
    with path.open() as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 2 or not fields[1].isdigit():
                raise ValueError(f"Genome sizes must be two-column headerless TSV: {path}")
            observed.append((fields[0], int(fields[1])))
    if observed != expected:
        raise ValueError(f"Genome-size content/order mismatch: {path}")


def validate_chromosome_ends(path: Path, lengths: dict[str, int]) -> None:
    with path.open(newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader, None)
        if header != ["#chr", "start", "end", "score", "type", "label"]:
            raise ValueError(f"Chromosome-end header mismatch: {path}")
        observed = [tuple(row) for row in reader]
    expected = []
    for chrom in CHROMS:
        expected.extend([
            (chrom, "0", "1", ".", "chromosome_end", "left"),
            (chrom, str(lengths[chrom] - 1), str(lengths[chrom]), ".", "chromosome_end", "right"),
        ])
    if observed != expected:
        raise ValueError(f"Chromosome-end coordinates/order mismatch: {path}")


def validate_tsv_columns(path: Path, column_count: int, minimum_data_rows: int = 0) -> int:
    with path.open(newline="") as handle:
        rows = list(csv.reader(handle, delimiter="\t"))
    if not rows or len(rows[0]) != column_count:
        raise ValueError(f"TSV header/column-count error: {path}")
    for row_number, row in enumerate(rows, 1):
        if len(row) != column_count or any("\t" in field or "\r" in field or "\n" in field for field in row):
            raise ValueError(f"TSV shape/control-character error at {path}:{row_number}")
    data_rows = len(rows) - 1
    if data_rows < minimum_data_rows:
        raise ValueError(f"TSV has too few data rows: {path}")
    return data_rows


def validate_chromosome_summary(path: Path, column_count: int) -> None:
    validate_tsv_columns(path, column_count, 16)
    with path.open(newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        next(reader)
        observed = [row[0] for row in reader]
    if observed != CHROMS:
        raise ValueError(f"Chromosome summary is not complete/naturally ordered: {path}")


def validate_name_map(path: Path, sample_lengths: dict[str, dict[str, int]]) -> None:
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    covered: dict[str, set[str]] = defaultdict(set)
    seen: set[tuple[str, str, str]] = set()
    for row in rows:
        sample, display = row["Sample"], row["Display_Name"]
        key = (sample, row["Original_Name"], row["Source_Type"])
        if sample not in sample_lengths or display not in sample_lengths[sample]:
            raise ValueError(f"Chromosome name-map sample/display mismatch: {path}")
        if int(row["Length_bp"]) != sample_lengths[sample][display] or row["Status"] != "Included" or key in seen:
            raise ValueError(f"Chromosome name-map length/status/duplicate error: {path}")
        seen.add(key)
        covered[sample].add(display)
    if any(covered[sample] != set(CHROMS) for sample in ("EG11", "FL")):
        raise ValueError(f"Chromosome name map does not cover chr01-chr16: {path}")


def validate_paf(
    path: Path,
    target_lengths: dict[str, int],
    query_lengths: dict[str, int],
    min_mapq: int | None = None,
    min_block: int | None = None,
) -> tuple[int, Counter[str]]:
    opener = gzip.open if path.suffix == ".gz" else open
    rows = 0
    strands: Counter[str] = Counter()
    previous_key: tuple[int, int, int, int, str] | None = None
    with opener(path, "rt", encoding="utf-8", newline="") as handle:
        for line_number, line in enumerate(handle, 1):
            if not line.strip() or "\r" in line:
                raise ValueError(f"Blank/CRLF PAF line at {path}:{line_number}")
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 12:
                raise ValueError(f"PAF has fewer than 12 fields at {path}:{line_number}")
            query, target = fields[0], fields[5]
            if query not in query_lengths or target not in target_lengths:
                raise ValueError(f"PAF direction/name mismatch at {path}:{line_number}")
            qlen, qstart, qend = int(fields[1]), int(fields[2]), int(fields[3])
            tlen, tstart, tend = int(fields[6]), int(fields[7]), int(fields[8])
            matching, block, mapq = int(fields[9]), int(fields[10]), int(fields[11])
            if qlen != query_lengths[query] or tlen != target_lengths[target]:
                raise ValueError(f"PAF sequence-length mismatch at {path}:{line_number}")
            if not (0 <= qstart < qend <= qlen and 0 <= tstart < tend <= tlen):
                raise ValueError(f"PAF coordinate error at {path}:{line_number}")
            if fields[4] not in {"+", "-"} or not (0 <= matching <= block) or not (0 <= mapq <= 255):
                raise ValueError(f"PAF strand/count/MAPQ error at {path}:{line_number}")
            if min_mapq is not None and mapq < min_mapq:
                raise ValueError(f"Filtered PAF MAPQ violation at {path}:{line_number}")
            if min_block is not None and block < min_block:
                raise ValueError(f"Filtered PAF block-length violation at {path}:{line_number}")
            key = (CHR_RANK[target], tstart, CHR_RANK[query], qstart, fields[4])
            if previous_key is not None and key < previous_key:
                raise ValueError(f"PAF sort error at {path}:{line_number}")
            previous_key = key
            strands[fields[4]] += 1
            rows += 1
    if rows == 0:
        raise ValueError(f"Empty PAF: {path}")
    return rows, strands


def validate_text(path: Path, allow_empty: bool = False, forbid_blank_lines: bool = True) -> None:
    data = path.read_bytes()
    if not data and allow_empty:
        return
    if not data:
        raise ValueError(f"Unexpected empty text file: {path}")
    if b"\r" in data or b"\x00" in data:
        raise ValueError(f"CR/NUL detected: {path}")
    if forbid_blank_lines and b"\n\n" in data:
        raise ValueError(f"Unexpected blank line: {path}")
    data.decode("utf-8")


def write_readme(stage: Path, window_size: int, display_gap_min: int, min_mapq: int, min_block: int) -> None:
    text = f"""# EG11–FL SyntenyViz input bundle

This bundle contains derived tracks for EG11 (target) and FL Africa hap2 (query). It does not contain the source FASTA/GFF files, software environments, temporary files, or any rendered figure.

The initial FL haplotype ambiguity was resolved by the user: option 1 (`FL_Africa_hap2`) was selected. The alternative `FL_American_hap1` is documented only in `metadata/Candidate_Resolution.tsv` and was not processed.

## Coordinate and naming conventions

- All chromosome names are normalized to `chr01` through `chr16`; original identifiers are retained in `metadata/chromosome_name_map.tsv`.
- Genome-size files have two tab-separated columns and no header.
- Gene-density bedGraph files have four tab-separated columns (`chrom`, `start`, `end`, `gene_count`) and no header. Coordinates are 0-based, half-open. Windows are {window_size:,} bp, with a shorter terminal window when required. A gene is counted once using `midpoint0 = (GFF_start + GFF_end - 2) // 2`.
- N-gap BED files have three tab-separated columns and use 0-based, half-open coordinates. Complete files contain every maximal FASTA N/n run; display files retain runs at least {display_gap_min} bp long.
- Chromosome-end TSV files use 0-based, half-open 1-bp intervals. They are generic chromosome ends, not confirmed telomeres. No validated telomere detection result, motif, or threshold was found in the bounded project search.
- PAF uses FL as query (columns 1–4) and EG11 as target (columns 6–9), with 0-based coordinates. The raw PAF is unfiltered but chromosome names are normalized. The filtered PAF requires MAPQ >= {min_mapq} and alignment block length >= {min_block:,} bp and retains both strands.

## Evidence boundaries

- FASTA N-runs are described only as assembly gaps or N-gaps.
- No `filled-gaps` BED is included because no AGP, old/new assembly coordinate mapping, or validated gap-filling coordinate record was found.
- The original PAF was generated with minimap2 2.30-r1287 using `-x asm5 -c -t 16`; it contains primary and secondary records. Filtering follows only the declared MAPQ/block thresholds.

## Key files

- `EG11/` and `FL/`: genome sizes, 200-kb gene density, complete/display N-gaps, chromosome-level summaries, and chromosome ends.
- `alignment/`: raw gzipped PAF, filtered PAF, chromosome-pair statistics, and primary correspondence.
- `metadata/`: input/output checksums, source inventory, chromosome-name mapping, excluded annotation seqids, feature/assembly/PAF statistics, full command log, scripts, and automated QC.
"""
    (stage / "README.md").write_text(text, encoding="utf-8", newline="\n")


def reproducible_tar(source_dir: Path, archive: Path) -> None:
    temp_tar = archive.with_suffix("")
    if temp_tar.exists() or archive.exists():
        raise FileExistsError(f"Refusing to overwrite archive: {archive}")
    base = "EG11_FL_synteny_inputs"
    with tarfile.open(temp_tar, "w") as tar:
        for path in sorted(source_dir.rglob("*"), key=lambda item: item.relative_to(source_dir).as_posix()):
            arcname = f"{base}/{path.relative_to(source_dir).as_posix()}"
            info = tar.gettarinfo(str(path), arcname=arcname)
            info.uid = info.gid = 0
            info.uname = info.gname = ""
            info.mtime = 0
            if path.is_file():
                with path.open("rb") as handle:
                    tar.addfile(info, handle)
            else:
                tar.addfile(info)
    with temp_tar.open("rb") as source, archive.open("wb") as destination:
        with gzip.GzipFile(filename="", mode="wb", fileobj=destination, mtime=0) as zipped:
            shutil.copyfileobj(source, zipped)
    temp_tar.unlink()
    with tarfile.open(archive, "r:gz") as tar:
        names = tar.getnames()
        if not names or any(name.startswith("/") or ".." in Path(name).parts for name in names):
            raise RuntimeError("Unsafe or empty archive")


def main() -> int:
    args = parse_args()
    project = args.project_root.resolve()
    final_stage = project / "results" / "EG11_FL_synteny_inputs"
    build_stage = project / "tmp" / f"EG11_FL_synteny_inputs.build.{os.getpid()}"
    archive = project / "deliverables" / "EG11_FL_synteny_inputs.tar.gz"
    archive_sha = project / "deliverables" / "EG11_FL_synteny_inputs.tar.gz.sha256"
    for output in (final_stage, build_stage, archive, archive_sha):
        if output.exists():
            raise FileExistsError(f"Refusing to overwrite: {output}")
    for required in (args.manifest, args.parameters, args.paf):
        if not required.is_file() or required.stat().st_size == 0:
            raise FileNotFoundError(required)
    if (args.window_size, args.display_gap_min, args.min_mapq, args.min_block) != (200_000, 100, 20, 50_000):
        raise ValueError("This fixed-name delivery workflow requires 200kb/100bp/MAPQ20/block50kb parameters")

    with args.manifest.open(newline="") as handle:
        manifest_rows = list(csv.DictReader(handle, delimiter="\t"))
    manifest = {row["Sample_ID"]: row for row in manifest_rows}
    if len(manifest_rows) != 2 or len(manifest) != 2 or set(manifest) != {"EG11", "FL"} or manifest["EG11"]["Role"] != "Target" or manifest["FL"]["Role"] != "Query":
        raise ValueError("Manifest must define EG11 target and FL query")
    with args.parameters.open(newline="") as handle:
        parameter_rows = list(csv.DictReader(handle, delimiter="\t"))
    parameters = {row["Parameter"]: row["Value"] for row in parameter_rows}
    if len(parameter_rows) != len(parameters):
        raise ValueError("Parameters.tsv contains duplicate parameter names")
    expected_parameters = {
        "Window_Size_bp": str(args.window_size),
        "Display_Gap_Min_bp": str(args.display_gap_min),
        "PAF_Min_MAPQ": str(args.min_mapq),
        "PAF_Min_Block_bp": str(args.min_block),
        "PAF_Target": "EG11",
        "PAF_Query": "FL",
    }
    if any(parameters.get(key) != value for key, value in expected_parameters.items()):
        raise ValueError("CLI arguments and Parameters.tsv are inconsistent")
    for row in manifest_rows:
        for field in ("Genome_FASTA", "Annotation_GFF", "Chromosome_Map"):
            path = Path(row[field])
            if not path.is_file() or path.stat().st_size == 0:
                raise FileNotFoundError(path)

    build_stage.parent.mkdir(parents=True, exist_ok=True)
    build_stage.mkdir()
    stage = build_stage
    metadata = stage / "metadata"
    metadata.mkdir()

    input_rows = []
    assembly_rows = []
    feature_rows = []
    excluded_rows = []
    name_rows = []
    sample_objects = {}
    for sample in ("EG11", "FL"):
        row = manifest[sample]
        fasta = Path(row["Genome_FASTA"])
        gff = Path(row["Annotation_GFF"])
        map_path = Path(row["Chromosome_Map"])
        for path in (fasta, gff, map_path):
            if not path.is_file() or path.stat().st_size == 0:
                raise FileNotFoundError(path)
        chromosomes, source_map = load_chromosome_map(map_path)
        observed, gaps, n_counts, invalid_counts, fasta_sha = scan_fasta(fasta, chromosomes)
        windows, feature_counts, all_feature_counts, excluded, gff_name_rows, gff_sha = parse_gff(
            gff, sample, chromosomes, source_map, row["Annotation_Map_Mode"], args.window_size
        )
        input_rows.append(file_record("Genome_FASTA", sample, fasta, fasta_sha))
        input_rows.append(file_record("Annotation_GFF", sample, gff, gff_sha))
        input_rows.append(file_record("Chromosome_Map", sample, map_path, sha256_file(map_path)))
        lengths = [item.length for item in chromosomes]
        n_total = sum(n_counts.values())
        assembly_rows.append([
            sample, len(lengths), sum(lengths), n50(lengths), max(lengths), n_total,
            f"{n_total/sum(lengths):.10f}", 0, "All FASTA sequences are mapped chromosome-level sequences",
        ])
        preferred_features = ["region", "gene", "pseudogene", "mRNA", "transcript", "exon", "CDS"]
        features_to_report = preferred_features + sorted(set(all_feature_counts) - set(preferred_features))
        for feature in features_to_report:
            feature_rows.append([
                sample,
                feature,
                all_feature_counts[feature],
                feature_counts[feature],
                all_feature_counts[feature] - feature_counts[feature],
            ])
        for seqid, count in sorted(excluded.items()):
            excluded_rows.append([sample, "GFF", seqid, count, "Excluded because seqid is not one of the 16 chromosome-level sequences"])
        for item in chromosomes:
            exact_gff_names = row["Annotation_Map_Mode"] == "GFF_seqid_equals_FASTA_seqid"
            name_rows.append([
                sample, item.source_id, item.display_id, item.length, "Included",
                "FASTA chromosome mapped by audited chromosome map; GFF seqid is identical" if exact_gff_names else "FASTA chromosome mapped by audited chromosome map",
                "FASTA_and_GFF" if exact_gff_names else "FASTA",
            ])
        name_rows.extend(gff_name_rows)
        write_sample_tracks(
            stage, sample, chromosomes, gaps, n_counts, invalid_counts, windows, feature_counts,
            args.window_size, args.display_gap_min,
        )
        sample_objects[sample] = (chromosomes, source_map, feature_counts, gaps, n_counts)

    paf_digest, paf_stats, paf_types = process_paf(
        args.paf, stage, sample_objects["EG11"][0], sample_objects["FL"][0], args.min_mapq, args.min_block
    )
    input_rows.append(file_record("Whole_Genome_PAF", "EG11_vs_FL", args.paf, paf_digest))

    write_tsv(
        metadata / "Input_Files_And_Checksums.tsv",
        ["Input_Type", "Sample_ID", "Logical_Path", "Resolved_Path", "Size_Bytes", "Modified_UTC", "SHA256"],
        input_rows,
    )
    write_tsv(
        metadata / "Assembly_Statistics.tsv",
        ["Sample_ID", "Sequence_Count", "Total_Length_bp", "N50_bp", "Longest_Sequence_bp", "N_bp", "N_Fraction", "Unlocalized_Sequence_Count", "Inclusion_Note"],
        assembly_rows,
    )
    write_tsv(
        metadata / "Annotation_Feature_Counts.tsv",
        ["Sample_ID", "Feature_Type", "Total_GFF_Count", "Included_Chromosome_Count", "Excluded_Sequence_Count"],
        feature_rows,
    )
    write_tsv(metadata / "Excluded_Annotation_Seqids.tsv", ["Sample_ID", "Source_Type", "Original_SeqID", "Feature_Count", "Exclusion_Reason"], excluded_rows)
    write_tsv(
        metadata / "chromosome_name_map.tsv",
        ["Sample", "Original_Name", "Display_Name", "Length_bp", "Status", "Note", "Source_Type"],
        sorted(name_rows, key=lambda row: (row[0], CHR_RANK[row[2]], row[6], row[1])),
    )
    paf_stat_rows = [[key, value] for key, value in sorted(paf_stats.items())]
    paf_stat_rows.extend([[key, value] for key, value in sorted(paf_types.items())])
    write_tsv(metadata / "PAF_Statistics.tsv", ["Metric", "Value"], paf_stat_rows)
    write_tsv(
        metadata / "Telomere_And_Filled_Gap_Evidence.tsv",
        ["Evidence_Type", "Status", "Source", "Method_or_Threshold", "Interpretation"],
        [
            ["Telomere", "Not_detected_in_bounded_inventory", "Current project and 00_final_input_resources", "No software/motif/threshold available", "Chromosome-end markers are not confirmed telomeres"],
            ["Filled_Gap", "Not_detected_in_bounded_inventory", "Current project and 00_final_input_resources", "No validated historical coordinate mapping", "No filled-gap BED was generated"],
        ],
    )
    write_tsv(
        metadata / "Input_Compatibility_Evidence.tsv",
        ["Sample_ID", "Comparison", "Status", "Evidence"],
        [
            ["EG11", "FASTA_vs_GFF", "PASS", "GFF build GCF_000442705.2; NC_025993.2-NC_026008.2 chromosome attributes and lengths map one-to-one to FASTA CM002081.2-CM002096.2"],
            ["FL", "FASTA_vs_GFF", "PASS", "All 16 GFF seqids are identical to FASTA seqids chr01B-chr16B and chromosome lengths are in bounds"],
            ["EG11_vs_FL", "PAF_direction", "PASS", "Original PAF query names/lengths belong to selected FL Africa hap2; target names/lengths belong to EG11"],
        ],
    )

    commands = [
        "# Complete build command",
        " ".join([
            "${DATA_DIR2}/anaconda3/bin/python3",
            str((project / "scripts" / "build_synteny_inputs.py").resolve()),
            "--manifest", str(args.manifest.resolve()),
            "--parameters", str(args.parameters.resolve()),
            "--paf", str(args.paf.resolve()),
            "--project-root", str(project),
            "--window-size", str(args.window_size),
            "--display-gap-min", str(args.display_gap_min),
            "--min-mapq", str(args.min_mapq),
            "--min-block", str(args.min_block),
        ]),
        "# Source whole-genome alignment command recovered from project log",
        "minimap2 -x asm5 -c -t 16 EG11_chromosomes.fa Africa_hap2.fasta > Orientation.paf",
        "# Archive verification",
        f"tar -tzf {archive}",
        f"sha256sum {archive}",
        "# Read-only preflight and SLURM submission",
        f"bash {(project / 'scripts' / '00_preflight.sh').resolve()}",
        f"sbatch {(project / 'scripts' / '10_build_inputs.slurm').resolve()}",
    ]
    (metadata / "Commands.log").write_text("\n".join(commands) + "\n", encoding="utf-8", newline="\n")
    shutil.copy2(project / "scripts" / "build_synteny_inputs.py", metadata / "build_synteny_inputs.py")
    shutil.copy2(project / "scripts" / "00_preflight.sh", metadata / "00_preflight.sh")
    shutil.copy2(project / "scripts" / "10_build_inputs.slurm", metadata / "10_build_inputs.slurm")
    shutil.copy2(project / "config" / "Discovery_Commands.log", metadata / "Discovery_Commands.log")
    shutil.copy2(args.manifest, metadata / "Input_Manifest.tsv")
    shutil.copy2(args.parameters, metadata / "Parameters.tsv")
    shutil.copy2(project / "config" / "Candidate_Resolution.tsv", metadata / "Candidate_Resolution.tsv")
    write_readme(stage, args.window_size, args.display_gap_min, args.min_mapq, args.min_block)

    qc_rows = []
    for sample in ("EG11", "FL"):
        chromosomes = sample_objects[sample][0]
        lengths = {item.display_id: item.length for item in chromosomes}
        genome_sizes = stage / sample / f"{sample}.genome_sizes.tsv"
        density_path = stage / sample / f"{sample}.gene_density.200kb.bedgraph"
        complete_gap = stage / sample / f"{sample}.N_gaps.complete.bed"
        display_gap = stage / sample / f"{sample}.N_gaps.display_min100bp.bed"
        chromosome_ends = stage / sample / f"{sample}.chromosome_ends.tsv"
        validate_genome_sizes(genome_sizes, chromosomes)
        density_n = validate_density(density_path, lengths, args.window_size)
        validate_chromosome_summary(stage / sample / f"{sample}.gene_density.chromosome_summary.tsv", 7)
        complete_n = validate_bed(complete_gap, lengths, 1)
        display_n = validate_bed(display_gap, lengths, args.display_gap_min)
        validate_chromosome_summary(stage / sample / f"{sample}.N_gaps.chromosome_summary.tsv", 7)
        complete_set = read_bed_set(complete_gap)
        display_set = read_bed_set(display_gap)
        if not display_set.issubset(complete_set):
            raise ValueError(f"Display gap BED is not a subset of complete BED: {sample}")
        gap_bp = sum(end - start for _, start, end in complete_set)
        if gap_bp != sum(sample_objects[sample][4].values()):
            raise ValueError(f"N-gap/FASTA N-base conservation failed: {sample}")
        validate_chromosome_ends(chromosome_ends, lengths)
        qc_rows.extend([
            [f"{sample}_Genome_Size", "PASS", "16 unique naturally ordered chromosomes; lengths equal audited FASTA"],
            [f"{sample}_Gene_Density", "PASS", f"{density_n} contiguous windows cover complete chromosomes"],
            [f"{sample}_Complete_N_Gaps", "PASS", f"{complete_n} non-overlapping maximal N-runs within chromosome bounds; total length equals FASTA N count"],
            [f"{sample}_Display_N_Gaps", "PASS", f"{display_n} complete-BED subset records with length >= {args.display_gap_min} bp"],
            [f"{sample}_Chromosome_Ends", "PASS", "32 one-base chromosome-end markers; not claimed as telomeres"],
        ])
    target_lengths = {item.display_id: item.length for item in sample_objects["EG11"][0]}
    query_lengths = {item.display_id: item.length for item in sample_objects["FL"][0]}
    raw_paf_n, raw_strands = validate_paf(
        stage / "alignment" / "EG11_vs_FL.raw.paf.gz", target_lengths, query_lengths
    )
    filtered_paf_n, filtered_strands = validate_paf(
        stage / "alignment" / "EG11_vs_FL.filtered.mapq20.block50kb.paf",
        target_lengths,
        query_lengths,
        args.min_mapq,
        args.min_block,
    )
    if not ({"+", "-"} <= set(raw_strands) and {"+", "-"} <= set(filtered_strands)):
        raise ValueError("Expected both forward and reverse PAF alignments")
    correspondence_qc = "PASS" if paf_stats["Correspondence_WARN_Count"] == 0 and paf_stats["Correspondence_FAIL_Count"] == 0 else "WARN"
    qc_rows.extend([
        ["PAF_Direction", "PASS", f"{raw_paf_n} raw records: query names/lengths belong to FL; target names/lengths belong to EG11"],
        ["PAF_Filter", "PASS", f"{filtered_paf_n} records satisfy MAPQ >= {args.min_mapq} and block length >= {args.min_block}; both strands retained"],
        ["PAF_Sort", "PASS", "Raw and filtered derivatives use target rank/start then query rank/start natural ordering"],
        ["Chromosome_Correspondence", correspondence_qc, f"Bidirectional checks: {paf_stats['Correspondence_PASS_Count']} PASS, {paf_stats['Correspondence_WARN_Count']} WARN, {paf_stats['Correspondence_FAIL_Count']} FAIL; see primary_correspondence.tsv"],
        ["Filled_Gaps", "PASS", "No file generated because validated historical coordinates were unavailable"],
        ["Coordinate_System", "PASS", "BED/bedGraph/PAF coordinates are 0-based half-open"],
        ["Chromosome_Namespace", "PASS", "All derivative files use chr01-chr16"],
    ])
    write_tsv(metadata / "QC_Report.tsv", ["Check", "Status", "Evidence"], qc_rows)

    tsv_checks = [
        (stage / "alignment" / "EG11_vs_FL.chromosome_pair_summary.tsv", 7, 1),
        (stage / "alignment" / "EG11_vs_FL.primary_correspondence.tsv", 12, 32),
        (metadata / "Input_Files_And_Checksums.tsv", 7, 7),
        (metadata / "Assembly_Statistics.tsv", 9, 2),
        (metadata / "Annotation_Feature_Counts.tsv", 5, 2),
        (metadata / "Excluded_Annotation_Seqids.tsv", 5, 1),
        (metadata / "chromosome_name_map.tsv", 7, 32),
        (metadata / "PAF_Statistics.tsv", 2, 1),
        (metadata / "Telomere_And_Filled_Gap_Evidence.tsv", 5, 2),
        (metadata / "Input_Compatibility_Evidence.tsv", 4, 3),
        (metadata / "Input_Manifest.tsv", 6, 2),
        (metadata / "Parameters.tsv", 3, 1),
        (metadata / "Candidate_Resolution.tsv", 10, 2),
        (metadata / "QC_Report.tsv", 3, 1),
    ]
    for path, columns, minimum_rows in tsv_checks:
        validate_tsv_columns(path, columns, minimum_rows)
    validate_name_map(
        metadata / "chromosome_name_map.tsv",
        {"EG11": target_lengths, "FL": query_lengths},
    )

    text_files = [path for path in stage.rglob("*") if path.is_file() and path.suffix != ".gz"]
    for path in text_files:
        allow_empty = path.name.endswith(".bed")
        forbid_blank_lines = path.suffix in {".tsv", ".bed", ".bedgraph", ".paf", ".sha256"}
        validate_text(path, allow_empty=allow_empty, forbid_blank_lines=forbid_blank_lines)

    checksum_rows = []
    for path in sorted([item for item in stage.rglob("*") if item.is_file()]):
        if path.name == "Output_Checksums.sha256":
            continue
        checksum_rows.append(f"{sha256_file(path)}  {path.relative_to(stage).as_posix()}")
    (metadata / "Output_Checksums.sha256").write_text("\n".join(checksum_rows) + "\n", encoding="ascii", newline="\n")
    validate_text(metadata / "Output_Checksums.sha256")
    for record in checksum_rows:
        expected_digest, relative = record.split("  ", 1)
        if sha256_file(stage / relative) != expected_digest:
            raise RuntimeError(f"Output checksum verification failed: {relative}")

    final_stage.parent.mkdir(parents=True, exist_ok=True)
    stage.replace(final_stage)
    stage = final_stage
    reproducible_tar(stage, archive)
    archive_digest = sha256_file(archive)
    archive_sha.write_text(f"{archive_digest}  {archive.name}\n", encoding="ascii", newline="\n")
    print(f"PASS archive={archive} size={archive.stat().st_size} sha256={archive_digest}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
