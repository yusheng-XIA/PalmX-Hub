#!/usr/bin/env python3
"""Count REF/ALT support once per query-name fragment and unique gene."""

from __future__ import annotations

import argparse
import collections
import csv
import gzip
import json
from pathlib import Path

import pysam


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--analysis", required=True, choices=("TN", "FL"))
    p.add_argument("--sample", required=True)
    p.add_argument("--sites", required=True, type=Path)
    p.add_argument("--bam", required=True)
    p.add_argument("--output", required=True, type=Path)
    p.add_argument("--summary", required=True, type=Path)
    p.add_argument("--mapq", type=int, default=20)
    p.add_argument("--baseq", type=int, default=20)
    return p.parse_args()


def load_sites(path: Path, analysis: str):
    sites: dict[str, dict[int, tuple[str, str, str]]] = collections.defaultdict(dict)
    genes: set[str] = set()
    with gzip.open(path, "rt") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"analysis", "chrom", "pos", "ref", "alt", "gene_id"}
        if not reader.fieldnames or not required.issubset(reader.fieldnames):
            raise SystemExit(f"invalid site-table header: {path}")
        for row in reader:
            if row["analysis"] != analysis:
                raise SystemExit(f"site table contains analysis {row['analysis']}, expected {analysis}")
            pos0 = int(row["pos"]) - 1
            value = (row["ref"].upper(), row["alt"].upper(), row["gene_id"])
            old = sites[row["chrom"]].setdefault(pos0, value)
            if old != value:
                raise SystemExit(f"conflicting duplicate site: {row['chrom']}:{row['pos']}")
            genes.add(row["gene_id"])
    if not genes:
        raise SystemExit("no informative sites loaded")
    return sites, genes


def process_fragment(
    records,
    sites,
    counts,
    qc,
    min_mapq: int,
    min_baseq: int,
) -> None:
    qc["query_name_groups"] += 1
    primary = [
        record
        for record in records
        if not record.is_secondary and not record.is_supplementary
    ]
    read1 = [record for record in primary if record.is_read1]
    read2 = [record for record in primary if record.is_read2]
    if len(read1) != 1 or len(read2) != 1:
        qc["discard_incomplete_or_duplicate_primary_pair"] += 1
        return
    pair = (read1[0], read2[0])
    if any(
        record.is_unmapped
        or not record.is_paired
        or not record.is_proper_pair
        or record.is_duplicate
        or record.is_qcfail
        for record in pair
    ):
        qc["discard_pair_flag_filter"] += 1
        return
    if any(record.mapping_quality < min_mapq for record in pair):
        qc["discard_pair_mapq"] += 1
        return

    observations: dict[str, set[str]] = collections.defaultdict(set)
    for record in pair:
        chrom_sites = sites.get(record.reference_name)
        if not chrom_sites:
            continue
        sequence = record.query_sequence
        qualities = record.query_qualities
        if sequence is None or qualities is None:
            qc["alignments_without_sequence_or_quality"] += 1
            continue
        for query_pos, ref_pos in record.get_aligned_pairs(matches_only=True):
            site = chrom_sites.get(ref_pos)
            if site is None:
                continue
            if qualities[query_pos] < min_baseq:
                qc["site_observations_below_baseq"] += 1
                continue
            ref, alt, gene = site
            base = sequence[query_pos].upper()
            qc["site_observations_passing_quality"] += 1
            if base == ref:
                observations[gene].add("REF")
            elif base == alt:
                observations[gene].add("ALT")
            else:
                observations[gene].add("OTHER")

    if not observations:
        qc["fragments_without_informative_observation"] += 1
        return
    if len(observations) != 1:
        qc["discard_multigene_fragment"] += 1
        return

    gene, alleles = next(iter(observations.items()))
    if "OTHER" in alleles:
        qc["discard_ambiguous_other_allele"] += 1
    elif alleles == {"REF"}:
        counts[gene]["ref_fragments"] += 1
        qc["retained_ref_fragments"] += 1
    elif alleles == {"ALT"}:
        counts[gene]["alt_fragments"] += 1
        qc["retained_alt_fragments"] += 1
    else:
        qc["discard_ref_alt_conflict"] += 1


def main() -> None:
    args = parse_args()
    if not args.sites.is_file() or args.sites.stat().st_size == 0:
        raise SystemExit(f"missing site table: {args.sites}")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.summary.parent.mkdir(parents=True, exist_ok=True)
    sites, genes = load_sites(args.sites, args.analysis)
    qc = collections.Counter()
    qc["analysis"] = args.analysis
    qc["sample"] = args.sample
    qc["mapq_threshold"] = args.mapq
    qc["baseq_threshold"] = args.baseq
    qc["site_table_genes"] = len(genes)
    qc["site_table_sites"] = sum(len(values) for values in sites.values())
    counts: dict[str, collections.Counter] = collections.defaultdict(collections.Counter)

    with pysam.AlignmentFile(args.bam, "rb") as bam:
        previous_name = None
        group = []
        for record in bam.fetch(until_eof=True):
            if previous_name is None:
                previous_name = record.query_name
            if record.query_name != previous_name:
                process_fragment(group, sites, counts, qc, args.mapq, args.baseq)
                group = []
                previous_name = record.query_name
            group.append(record)
        if group:
            process_fragment(group, sites, counts, qc, args.mapq, args.baseq)

    with args.output.open("w", newline="") as out:
        writer = csv.writer(out, delimiter="\t", lineterminator="\n")
        writer.writerow(
            (
                "analysis",
                "sample",
                "gene_id",
                "ref_fragments",
                "alt_fragments",
                "informative_fragments",
            )
        )
        for gene in sorted(counts):
            ref_count = counts[gene]["ref_fragments"]
            alt_count = counts[gene]["alt_fragments"]
            writer.writerow(
                (
                    args.analysis,
                    args.sample,
                    gene,
                    ref_count,
                    alt_count,
                    ref_count + alt_count,
                )
            )
    qc["genes_with_retained_fragments"] = len(counts)
    qc["retained_informative_fragments"] = (
        qc["retained_ref_fragments"] + qc["retained_alt_fragments"]
    )
    with args.summary.open("w") as out:
        json.dump(dict(sorted(qc.items())), out, indent=2)
        out.write("\n")


if __name__ == "__main__":
    main()
