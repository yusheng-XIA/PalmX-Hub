#!/usr/bin/env python3
"""Build FL Africa hap2 BISER segmental-duplication BED and 200-kb density files."""

from __future__ import annotations

import argparse
import csv
import hashlib
import shlex
import statistics
import sys
from collections import defaultdict
from pathlib import Path


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def write_tsv(path: Path, fields: list[str], rows) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def read_sizes(path: Path) -> tuple[list[str], dict[str, int]]:
    order, sizes = [], {}
    with path.open(encoding="utf-8") as handle:
        for line_no, line in enumerate(handle, 1):
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 2:
                raise ValueError(f"Genome-size line {line_no} does not have two tab-separated fields")
            chrom, length_text = fields
            length = int(length_text)
            if chrom in sizes or length <= 0:
                raise ValueError(f"Invalid genome-size entry at line {line_no}")
            order.append(chrom)
            sizes[chrom] = length
    return order, sizes


def display_chrom(name: str, sizes: dict[str, int]) -> str:
    if name in sizes:
        return name
    prefix = "seedless_hap2_chr"
    if name.startswith(prefix) and name.endswith("B"):
        chrom = "chr" + name[len(prefix):-1]
        if chrom in sizes:
            return chrom
    raise ValueError(f"BISER chromosome is not an accepted FL Africa hap2 chromosome: {name}")


def merge_intervals(intervals: list[tuple[int, int]]) -> list[tuple[int, int]]:
    if not intervals:
        return []
    merged: list[list[int]] = []
    for start, end in sorted(intervals):
        if merged and start <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], end)
        else:
            merged.append([start, end])
    return [(start, end) for start, end in merged]


def add_coverage(intervals: list[tuple[int, int]], chrom_length: int, window_size: int) -> list[int]:
    values = [0] * ((chrom_length + window_size - 1) // window_size)
    for start, end in intervals:
        for index in range(start // window_size, (end - 1) // window_size + 1):
            window_start = index * window_size
            window_end = min(window_start + window_size, chrom_length)
            values[index] += max(0, min(end, window_end) - max(start, window_start))
    return values


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--biser", required=True, type=Path)
    parser.add_argument("--genome-sizes", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--min-length", required=True, type=int)
    parser.add_argument("--max-score", required=True, type=float)
    parser.add_argument("--window-size", required=True, type=int)
    args = parser.parse_args()

    for path in (args.biser, args.genome_sizes):
        if not path.is_file() or path.stat().st_size == 0:
            raise FileNotFoundError(path)
    if args.min_length <= 0 or args.max_score < 0 or args.window_size <= 0:
        raise ValueError("Thresholds and window size must be valid positive values")
    if args.output.exists():
        raise FileExistsError(f"Refusing to overwrite existing output: {args.output}")
    args.output.mkdir(parents=True)

    order, sizes = read_sizes(args.genome_sizes)
    rank = {chrom: index for index, chrom in enumerate(order)}
    sides: dict[str, list[tuple[int, int, str, str]]] = {chrom: [] for chrom in order}
    pair_rows = []
    raw_pairs = filtered_pairs = intra_pairs = inter_pairs = 0

    with args.biser.open(encoding="utf-8", errors="strict") as handle:
        for line_no, line in enumerate(handle, 1):
            if not line.strip() or line.startswith("#"):
                continue
            raw_pairs += 1
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 14:
                raise ValueError(f"BISER line {line_no} has {len(fields)} fields instead of 14")
            try:
                chrom_1 = display_chrom(fields[0], sizes)
                start_1, end_1 = int(fields[1]), int(fields[2])
                chrom_2 = display_chrom(fields[3], sizes)
                start_2, end_2 = int(fields[4]), int(fields[5])
                error_percent = float(fields[7])
            except ValueError as error:
                raise ValueError(f"Invalid BISER record at line {line_no}: {error}") from error
            length_1, length_2 = end_1 - start_1, end_2 - start_2
            if start_1 < 0 or end_1 <= start_1 or end_1 > sizes[chrom_1]:
                raise ValueError(f"Out-of-range side 1 at BISER line {line_no}")
            if start_2 < 0 or end_2 <= start_2 or end_2 > sizes[chrom_2]:
                raise ValueError(f"Out-of-range side 2 at BISER line {line_no}")
            if length_1 < args.min_length or length_2 < args.min_length or error_percent > args.max_score:
                continue

            filtered_pairs += 1
            pair_id = f"FLSD{filtered_pairs:06d}"
            pair_class = "Intra" if chrom_1 == chrom_2 else "Inter"
            intra_pairs += pair_class == "Intra"
            inter_pairs += pair_class == "Inter"
            sides[chrom_1].append((start_1, end_1, pair_id + "_A", fields[8]))
            sides[chrom_2].append((start_2, end_2, pair_id + "_B", fields[9]))
            pair_rows.append({
                "Pair_ID": pair_id,
                "Chrom_1": chrom_1,
                "Start_1": start_1,
                "End_1": end_1,
                "Length_1_Bp": length_1,
                "Chrom_2": chrom_2,
                "Start_2": start_2,
                "End_2": end_2,
                "Length_2_Bp": length_2,
                "Reference_Pair": fields[6],
                "Error_Percent": fields[7],
                "Strand_1": fields[8],
                "Strand_2": fields[9],
                "BISER_Field_11": fields[10],
                "BISER_Field_12": fields[11],
                "CIGAR": fields[12],
                "Tags": fields[13],
                "Pair_Class": pair_class,
                "Coordinate_System": "0-based_half-open",
            })

    if not pair_rows:
        raise ValueError("No BISER pairs passed the selected thresholds")
    pair_rows.sort(key=lambda row: (
        rank[row["Chrom_1"]], row["Start_1"], row["End_1"],
        rank[row["Chrom_2"]], row["Start_2"], row["End_2"],
    ))
    pair_fields = list(pair_rows[0])
    write_tsv(args.output / "FL_Africa_hap2.BISER_SD_pairs.filtered.tsv", pair_fields, pair_rows)

    sites_bed = args.output / "FL_Africa_hap2.BISER_SD_sites.bed"
    union_bed = args.output / "FL_Africa_hap2.BISER_SD_union.bed"
    unions: dict[str, list[tuple[int, int]]] = {}
    with sites_bed.open("w", encoding="utf-8") as sites, union_bed.open("w", encoding="utf-8") as union:
        for chrom in order:
            sides[chrom].sort(key=lambda item: (item[0], item[1], item[2]))
            for start, end, name, strand in sides[chrom]:
                sites.write(f"{chrom}\t{start}\t{end}\t{name}\t0\t{strand}\n")
            unions[chrom] = merge_intervals([(start, end) for start, end, _, _ in sides[chrom]])
            for start, end in unions[chrom]:
                union.write(f"{chrom}\t{start}\t{end}\n")

    density_path = args.output / "FL_Africa_hap2.BISER_SD_density.200kb.tsv"
    fraction_path = args.output / "FL_Africa_hap2.BISER_SD_fraction.200kb.bedgraph"
    percent_path = args.output / "FL_Africa_hap2.BISER_SD_percent.200kb.bedgraph"
    covered_path = args.output / "FL_Africa_hap2.BISER_SD_covered_bp.200kb.bedgraph"
    density_fields = [
        "Chrom", "Start", "End", "Window_Bp", "SD_Covered_Bp", "SD_Fraction", "SD_Percent",
        "SD_Union_Locus_Midpoint_Count", "SD_Pair_Side_Midpoint_Count",
    ]
    density_rows = []
    chromosome_rows = []
    total_union_bp = total_union_loci = 0
    with fraction_path.open("w", encoding="utf-8") as fraction, \
            percent_path.open("w", encoding="utf-8") as percent, \
            covered_path.open("w", encoding="utf-8") as covered:
        for chrom in order:
            coverage = add_coverage(unions[chrom], sizes[chrom], args.window_size)
            union_counts = [0] * len(coverage)
            side_counts = [0] * len(coverage)
            for start, end in unions[chrom]:
                union_counts[((start + end) // 2) // args.window_size] += 1
            for start, end, _, _ in sides[chrom]:
                side_counts[((start + end) // 2) // args.window_size] += 1
            chrom_union_bp = sum(coverage)
            total_union_bp += chrom_union_bp
            total_union_loci += len(unions[chrom])
            chromosome_rows.append({
                "Chrom": chrom,
                "Chromosome_Bp": sizes[chrom],
                "Filtered_SD_Pair_Side_Count": len(sides[chrom]),
                "Merged_SD_Union_Interval_Count": len(unions[chrom]),
                "SD_Union_Covered_Bp": chrom_union_bp,
                "SD_Union_Fraction": f"{chrom_union_bp / sizes[chrom]:.8f}",
                "SD_Union_Percent": f"{100 * chrom_union_bp / sizes[chrom]:.6f}",
            })
            for index, covered_bp in enumerate(coverage):
                start = index * args.window_size
                end = min(start + args.window_size, sizes[chrom])
                window_bp = end - start
                value_fraction = covered_bp / window_bp
                value_percent = 100 * value_fraction
                density_rows.append({
                    "Chrom": chrom,
                    "Start": start,
                    "End": end,
                    "Window_Bp": window_bp,
                    "SD_Covered_Bp": covered_bp,
                    "SD_Fraction": f"{value_fraction:.8f}",
                    "SD_Percent": f"{value_percent:.6f}",
                    "SD_Union_Locus_Midpoint_Count": union_counts[index],
                    "SD_Pair_Side_Midpoint_Count": side_counts[index],
                })
                fraction.write(f"{chrom}\t{start}\t{end}\t{value_fraction:.8f}\n")
                percent.write(f"{chrom}\t{start}\t{end}\t{value_percent:.6f}\n")
                covered.write(f"{chrom}\t{start}\t{end}\t{covered_bp}\n")

    write_tsv(density_path, density_fields, density_rows)
    summary_fields = [
        "Chrom", "Chromosome_Bp", "Filtered_SD_Pair_Side_Count", "Merged_SD_Union_Interval_Count",
        "SD_Union_Covered_Bp", "SD_Union_Fraction", "SD_Union_Percent",
    ]
    write_tsv(args.output / "FL_Africa_hap2.BISER_SD_summary_by_chromosome.tsv", summary_fields, chromosome_rows)

    genome_bp = sum(sizes.values())
    genome_summary = [{
        "Sample": "FL_Africa_hap2",
        "Genome_Bp": genome_bp,
        "Raw_BISER_Pair_Count": raw_pairs,
        "Filtered_BISER_Pair_Count": filtered_pairs,
        "Filtered_Intra_Chromosomal_Pair_Count": intra_pairs,
        "Filtered_Inter_Chromosomal_Pair_Count": inter_pairs,
        "Filtered_SD_Pair_Side_Count": sum(len(items) for items in sides.values()),
        "Merged_SD_Union_Interval_Count": total_union_loci,
        "SD_Union_Covered_Bp": total_union_bp,
        "SD_Union_Fraction": f"{total_union_bp / genome_bp:.8f}",
        "SD_Union_Percent": f"{100 * total_union_bp / genome_bp:.6f}",
    }]
    write_tsv(args.output / "FL_Africa_hap2.BISER_SD_summary_genome.tsv", list(genome_summary[0]), genome_summary)

    stats_rows = []
    for chrom in order:
        values = [float(row["SD_Fraction"]) for row in density_rows if row["Chrom"] == chrom]
        stats_rows.append({
            "Chrom": chrom,
            "Window_Count": len(values),
            "Maximum_SD_Fraction": f"{max(values):.8f}",
            "Mean_SD_Fraction": f"{statistics.mean(values):.8f}",
            "Median_SD_Fraction": f"{statistics.median(values):.8f}",
        })
    write_tsv(
        args.output / "FL_Africa_hap2.BISER_SD_density_statistics.tsv",
        ["Chrom", "Window_Count", "Maximum_SD_Fraction", "Mean_SD_Fraction", "Median_SD_Fraction"],
        stats_rows,
    )

    expected_windows = sum((length + args.window_size - 1) // args.window_size for length in sizes.values())
    if len(density_rows) != expected_windows:
        raise ValueError("Density window count mismatch")
    if sum(int(row["SD_Covered_Bp"]) for row in density_rows) != total_union_bp:
        raise ValueError("Window coverage does not sum to genome-wide SD union coverage")
    if sum(int(row["SD_Union_Locus_Midpoint_Count"]) for row in density_rows) != total_union_loci:
        raise ValueError("Merged-locus midpoint counts do not sum to the union interval count")
    if sum(int(row["SD_Pair_Side_Midpoint_Count"]) for row in density_rows) != 2 * filtered_pairs:
        raise ValueError("Pair-side midpoint counts do not sum to twice the filtered pair count")
    if any(int(row["SD_Covered_Bp"]) > int(row["Window_Bp"]) for row in density_rows):
        raise ValueError("SD coverage exceeds window length")

    parameters = [
        {"Parameter": "Sample", "Value": "FL_Africa_hap2", "Unit": "text"},
        {"Parameter": "BISER_Version", "Value": "1.4", "Unit": "text"},
        {"Parameter": "Minimum_Length_Each_Side", "Value": args.min_length, "Unit": "bp"},
        {"Parameter": "Maximum_BISER_Error_Score", "Value": args.max_score, "Unit": "percent"},
        {"Parameter": "Window_Size", "Value": args.window_size, "Unit": "bp"},
        {"Parameter": "Coordinate_System", "Value": "0-based_half-open", "Unit": "text"},
        {"Parameter": "Density_Definition", "Value": "merged_SD_union_bp_divided_by_window_bp", "Unit": "text"},
        {"Parameter": "Chromosome_Name_Conversion", "Value": "seedless_hap2_chrNNB_to_chrNN", "Unit": "text"},
    ]
    write_tsv(args.output / "parameters.tsv", ["Parameter", "Value", "Unit"], parameters)

    command = shlex.join([sys.executable, str(Path(__file__).resolve()), *sys.argv[1:]])
    (args.output / "commands.log").write_text(command + "\n", encoding="utf-8")
    (args.output / "README.md").write_text(f"""# FL Africa hap2 BISER segmental duplications

The BISER v1.4 source is `seedless_hap2_SD`. Its source assembly has the same 16 chromosome lengths as the accepted FL Africa hap2 genome. Original names `seedless_hap2_chr01B`–`seedless_hap2_chr16B` were converted to display names `chr01`–`chr16`.

Pairs were retained only when both sides were at least {args.min_length:,} bp and the BISER error score was at most {args.max_score:g}. All coordinates are 0-based half-open. `FL_Africa_hap2.BISER_SD_sites.bed` is BED6 with both sides of each retained pair. `FL_Africa_hap2.BISER_SD_union.bed` is the merged, nonredundant BED3 used to calculate density.

`FL_Africa_hap2.BISER_SD_density.200kb.tsv` is the comprehensive fixed-window table. The recommended plotting track is `FL_Africa_hap2.BISER_SD_fraction.200kb.bedgraph`; it reports nonredundant SD-covered bp divided by actual window length. The final window of each chromosome may be shorter than {args.window_size:,} bp, and all zero-density windows are explicit.

The analysis retained {filtered_pairs:,} of {raw_pairs:,} BISER pairs. Their two sides merged into {total_union_loci:,} nonredundant intervals covering {total_union_bp:,} bp ({100 * total_union_bp / genome_bp:.4f}% of the FL hap2 assembly).
""", encoding="utf-8")
    (args.output / "QC_Report.md").write_text(f"""# QC report

- Input identity: the 16 BISER chromosome names map uniquely to the accepted FL Africa hap2 chromosomes.
- Raw BISER pairs: {raw_pairs:,}; retained pairs: {filtered_pairs:,}.
- Retained pair sides: {2 * filtered_pairs:,}; merged union intervals: {total_union_loci:,}.
- Density windows: {len(density_rows):,}; expected windows: {expected_windows:,}.
- All coordinates are within chromosome bounds and use 0-based half-open convention.
- Windows are continuous, nonoverlapping, chromosome ordered, and include explicit zeros.
- Window coverage sums to the merged genome-wide SD coverage ({total_union_bp:,} bp).
- No window has SD coverage greater than its window length.
- Status: PASS.
""", encoding="utf-8")

    write_tsv(
        args.output / "input_checksums.tsv",
        ["Input_Type", "Size_Bytes", "SHA256", "Absolute_Path"],
        ({
            "Input_Type": label,
            "Size_Bytes": path.stat().st_size,
            "SHA256": sha256(path),
            "Absolute_Path": str(path.resolve()),
        } for label, path in (("BISER_Raw_Output", args.biser), ("FL_Genome_Sizes", args.genome_sizes))),
    )

    checksum_path = args.output / "output_checksums.sha256"
    output_files = sorted(path for path in args.output.iterdir() if path.is_file() and path != checksum_path)
    with checksum_path.open("w", encoding="utf-8") as handle:
        for path in output_files:
            handle.write(f"{sha256(path)}  {path.name}\n")
    for path in output_files:
        data = path.read_bytes()
        if b"\r" in data or b"\x00" in data:
            raise ValueError(f"Illegal CR or NUL character in {path}")

    print(f"OUTPUT\t{args.output.resolve()}")
    print(f"RAW_PAIRS\t{raw_pairs}")
    print(f"FILTERED_PAIRS\t{filtered_pairs}")
    print(f"UNION_INTERVALS\t{total_union_loci}")
    print(f"SD_COVERED_BP\t{total_union_bp}")
    print(f"SD_PERCENT\t{100 * total_union_bp / genome_bp:.6f}")
    print(f"WINDOWS\t{len(density_rows)}")


if __name__ == "__main__":
    main()
