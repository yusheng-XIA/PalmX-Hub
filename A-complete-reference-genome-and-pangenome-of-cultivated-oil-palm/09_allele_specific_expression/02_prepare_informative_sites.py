#!/usr/bin/env python3
"""Assign biallelic diagnostic SNPs to exactly one exonic backbone gene."""

from __future__ import annotations

import argparse
import collections
import csv
import gzip
import heapq
import json
from pathlib import Path

import pysam


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--analysis", required=True, choices=("TN", "FL"))
    p.add_argument("--vcf", required=True, type=Path)
    p.add_argument("--gff", required=True, type=Path)
    p.add_argument("--output", required=True, type=Path)
    p.add_argument("--summary", required=True, type=Path)
    return p.parse_args()


def attributes(text: str) -> dict[str, str]:
    result: dict[str, str] = {}
    for item in text.rstrip().split(";"):
        if not item:
            continue
        key, sep, value = item.partition("=")
        if sep:
            result[key] = value
    return result


def load_exons(gff: Path):
    tx_to_genes: dict[str, set[str]] = {}
    exons: list[tuple[str, int, int, tuple[str, ...]]] = []
    stats = collections.Counter()

    with gff.open() as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                stats["malformed_gff_rows"] += 1
                continue
            chrom, _, feature, start, end, _, _, _, attr_text = fields
            attrs = attributes(attr_text)
            if feature in {"mRNA", "transcript"}:
                tx = attrs.get("ID")
                parents = {x for x in attrs.get("Parent", "").split(",") if x}
                if tx and parents:
                    tx_to_genes.setdefault(tx, set()).update(parents)
            elif feature == "exon":
                parents = tuple(x for x in attrs.get("Parent", "").split(",") if x)
                if parents:
                    exons.append((chrom, int(start) - 1, int(end), parents))
                else:
                    stats["exons_without_parent"] += 1

    intervals: dict[str, set[tuple[int, int, str]]] = collections.defaultdict(set)
    for chrom, start0, end0, parents in exons:
        genes: set[str] = set()
        for parent in parents:
            genes.update(tx_to_genes.get(parent, ()))
        if len(genes) != 1:
            stats["exons_unresolved_or_multigene"] += 1
            continue
        intervals[chrom].add((start0, end0, next(iter(genes))))
        stats["resolved_exon_rows"] += 1

    ordered = {chrom: sorted(values) for chrom, values in intervals.items()}
    stats["transcripts_with_gene_parent"] = len(tx_to_genes)
    stats["unique_exon_gene_intervals"] = sum(map(len, ordered.values()))
    stats["gff_contigs_with_resolved_exons"] = len(ordered)
    return ordered, stats


def main() -> None:
    args = parse_args()
    for path in (args.vcf, args.gff):
        if not path.is_file() or path.stat().st_size == 0:
            raise SystemExit(f"missing or empty input: {path}")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.summary.parent.mkdir(parents=True, exist_ok=True)

    intervals, stats = load_exons(args.gff)
    stats["analysis"] = args.analysis
    stats["vcf"] = str(args.vcf)
    stats["gff"] = str(args.gff)

    with pysam.VariantFile(str(args.vcf)) as vcf, gzip.open(
        args.output, "wt", newline=""
    ) as out:
        writer = csv.writer(out, delimiter="\t", lineterminator="\n")
        writer.writerow(("analysis", "chrom", "pos", "ref", "alt", "gene_id"))

        for chrom in vcf.header.contigs:
            chrom_intervals = intervals.get(chrom)
            if not chrom_intervals:
                continue
            next_interval = 0
            active_heap: list[tuple[int, int, str]] = []
            active_genes: collections.Counter[str] = collections.Counter()
            interval_id = 0

            for record in vcf.fetch(chrom):
                stats["vcf_records_examined"] += 1
                if (
                    len(record.ref) != 1
                    or not record.alts
                    or len(record.alts) != 1
                    or len(record.alts[0]) != 1
                ):
                    stats["non_biallelic_snv_records"] += 1
                    continue
                pos0 = record.pos - 1
                while (
                    next_interval < len(chrom_intervals)
                    and chrom_intervals[next_interval][0] <= pos0
                ):
                    start0, end0, gene = chrom_intervals[next_interval]
                    next_interval += 1
                    if end0 <= pos0:
                        continue
                    interval_id += 1
                    heapq.heappush(active_heap, (end0, interval_id, gene))
                    active_genes[gene] += 1
                while active_heap and active_heap[0][0] <= pos0:
                    _, _, gene = heapq.heappop(active_heap)
                    active_genes[gene] -= 1
                    if active_genes[gene] == 0:
                        del active_genes[gene]

                if len(active_genes) == 1:
                    gene = next(iter(active_genes))
                    writer.writerow(
                        (
                            args.analysis,
                            chrom,
                            record.pos,
                            record.ref.upper(),
                            record.alts[0].upper(),
                            gene,
                        )
                    )
                    stats["retained_unique_gene_snps"] += 1
                elif len(active_genes) > 1:
                    stats["excluded_overlapping_gene_snps"] += 1
                else:
                    stats["nonexonic_snps"] += 1

    if stats["retained_unique_gene_snps"] == 0:
        args.output.unlink(missing_ok=True)
        raise SystemExit("no unique-gene exonic SNPs retained")
    with args.summary.open("w") as handle:
        json.dump(dict(sorted(stats.items())), handle, indent=2)
        handle.write("\n")


if __name__ == "__main__":
    main()
