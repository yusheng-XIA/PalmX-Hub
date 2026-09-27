#!/usr/bin/env python3
"""Build FL Africa hap2 TRF tandem-repeat BED and fixed-window density tracks."""

from __future__ import annotations

import argparse
import csv
import hashlib
import shlex
import statistics
import sys
from pathlib import Path


EXPECTED_PARAMETERS = "2 6 6 80 10 50 2000"


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
                raise ValueError(f"Genome-size line {line_no} must have two tab-separated fields")
            chrom, length_text = fields
            length = int(length_text)
            if chrom in sizes or length <= 0:
                raise ValueError(f"Invalid genome-size entry at line {line_no}")
            order.append(chrom)
            sizes[chrom] = length
    return order, sizes


def display_chrom(name: str, sizes: dict[str, int]) -> str:
    chrom = name[:-1] if name.startswith("chr") and name.endswith("B") else name
    if chrom not in sizes:
        raise ValueError(f"Unknown FL Africa hap2 chromosome in TRF output: {name}")
    return chrom


def merge_intervals(intervals: list[tuple[int, int]]) -> list[tuple[int, int]]:
    merged: list[list[int]] = []
    for start, end in sorted(intervals):
        if merged and start <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], end)
        else:
            merged.append([start, end])
    return [(start, end) for start, end in merged]


def coverage_by_window(
    intervals: list[tuple[int, int]], chrom_length: int, window_size: int
) -> list[int]:
    values = [0] * ((chrom_length + window_size - 1) // window_size)
    for start, end in intervals:
        for index in range(start // window_size, (end - 1) // window_size + 1):
            window_start = index * window_size
            window_end = min(window_start + window_size, chrom_length)
            values[index] += max(0, min(end, window_end) - max(start, window_start))
    return values


def contains_cr_or_nul(path: Path) -> bool:
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            if b"\r" in chunk or b"\x00" in chunk:
                return True
    return False


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--trf-dat", required=True, type=Path)
    parser.add_argument("--genome-sizes", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--window-size", required=True, type=int)
    args = parser.parse_args()

    assert merge_intervals([(10, 20), (15, 25), (30, 40)]) == [(10, 25), (30, 40)]
    for path in (args.trf_dat, args.genome_sizes):
        if not path.is_file() or path.stat().st_size == 0:
            raise FileNotFoundError(path)
    if args.window_size <= 0:
        raise ValueError("Window size must be positive")
    if args.output.exists():
        raise FileExistsError(f"Refusing to overwrite existing output: {args.output}")
    args.output.mkdir(parents=True)

    order, sizes = read_sizes(args.genome_sizes)
    intervals: dict[str, list[tuple[int, int]]] = {chrom: [] for chrom in order}
    current_chrom = ""
    sequence_headers = raw_calls = 0
    parameter_sets: set[str] = set()
    min_score = max_score = max_period = None

    with args.trf_dat.open(encoding="utf-8", errors="strict") as handle:
        for line_no, line in enumerate(handle, 1):
            if line.startswith("Sequence: "):
                current_chrom = display_chrom(line.split(None, 1)[1].strip(), sizes)
                sequence_headers += 1
                continue
            if line.startswith("Parameters: "):
                parameter_sets.add(line.split(":", 1)[1].strip())
                continue
            if not line or not line[0].isdigit():
                continue
            if not current_chrom:
                raise ValueError(f"TRF data record before a Sequence header at line {line_no}")
            fields = line.split(None, 13)
            if len(fields) < 14:
                raise ValueError(f"Malformed TRF data record at line {line_no}")
            try:
                start_1based, end_1based = int(fields[0]), int(fields[1])
                period, score = int(fields[2]), int(fields[7])
            except ValueError as error:
                raise ValueError(f"Invalid numeric TRF field at line {line_no}: {error}") from error
            start, end = start_1based - 1, end_1based
            if start < 0 or end <= start or end > sizes[current_chrom]:
                raise ValueError(f"Out-of-range TRF interval at line {line_no}")
            intervals[current_chrom].append((start, end))
            raw_calls += 1
            min_score = score if min_score is None else min(min_score, score)
            max_score = score if max_score is None else max(max_score, score)
            max_period = period if max_period is None else max(max_period, period)

    if sequence_headers != len(order) or set(chrom for chrom, values in intervals.items() if values) != set(order):
        raise ValueError("TRF output does not contain all 16 accepted FL chromosomes")
    if parameter_sets != {EXPECTED_PARAMETERS}:
        raise ValueError(f"Unexpected TRF parameters: {sorted(parameter_sets)}")
    if raw_calls == 0:
        raise ValueError("No TRF tandem-repeat calls were parsed")

    calls_bed = args.output / "FL_Africa_hap2.TRF_tandem_repeat_calls.bed"
    union_bed = args.output / "FL_Africa_hap2.TRF_tandem_repeat_union.bed"
    unions: dict[str, list[tuple[int, int]]] = {}
    call_id = 0
    with calls_bed.open("w", encoding="utf-8") as calls, union_bed.open("w", encoding="utf-8") as union:
        for chrom in order:
            intervals[chrom].sort()
            for start, end in intervals[chrom]:
                call_id += 1
                calls.write(f"{chrom}\t{start}\t{end}\tFLTRF{call_id:07d}\n")
            unions[chrom] = merge_intervals(intervals[chrom])
            for start, end in unions[chrom]:
                union.write(f"{chrom}\t{start}\t{end}\n")

    density_path = args.output / "FL_Africa_hap2.TRF_tandem_repeat_density.200kb.tsv"
    fraction_path = args.output / "FL_Africa_hap2.TRF_tandem_repeat_fraction.200kb.bedgraph"
    percent_path = args.output / "FL_Africa_hap2.TRF_tandem_repeat_percent.200kb.bedgraph"
    covered_path = args.output / "FL_Africa_hap2.TRF_tandem_repeat_covered_bp.200kb.bedgraph"
    density_fields = [
        "Chrom", "Start", "End", "Window_Bp", "Tandem_Repeat_Covered_Bp",
        "Tandem_Repeat_Fraction", "Tandem_Repeat_Percent",
        "Merged_Locus_Midpoint_Count", "Raw_TRF_Call_Midpoint_Count",
    ]
    density_rows = []
    chromosome_rows = []
    total_union_bp = total_union_loci = 0

    with fraction_path.open("w", encoding="utf-8") as fraction, \
            percent_path.open("w", encoding="utf-8") as percent, \
            covered_path.open("w", encoding="utf-8") as covered:
        for chrom in order:
            coverage = coverage_by_window(unions[chrom], sizes[chrom], args.window_size)
            union_counts = [0] * len(coverage)
            raw_counts = [0] * len(coverage)
            for start, end in unions[chrom]:
                union_counts[((start + end) // 2) // args.window_size] += 1
            for start, end in intervals[chrom]:
                raw_counts[((start + end) // 2) // args.window_size] += 1
            chrom_union_bp = sum(coverage)
            total_union_bp += chrom_union_bp
            total_union_loci += len(unions[chrom])
            chromosome_rows.append({
                "Chrom": chrom,
                "Chromosome_Bp": sizes[chrom],
                "Raw_TRF_Call_Count": len(intervals[chrom]),
                "Merged_Tandem_Repeat_Locus_Count": len(unions[chrom]),
                "Tandem_Repeat_Union_Covered_Bp": chrom_union_bp,
                "Tandem_Repeat_Union_Fraction": f"{chrom_union_bp / sizes[chrom]:.8f}",
                "Tandem_Repeat_Union_Percent": f"{100 * chrom_union_bp / sizes[chrom]:.6f}",
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
                    "Tandem_Repeat_Covered_Bp": covered_bp,
                    "Tandem_Repeat_Fraction": f"{value_fraction:.8f}",
                    "Tandem_Repeat_Percent": f"{value_percent:.6f}",
                    "Merged_Locus_Midpoint_Count": union_counts[index],
                    "Raw_TRF_Call_Midpoint_Count": raw_counts[index],
                })
                fraction.write(f"{chrom}\t{start}\t{end}\t{value_fraction:.8f}\n")
                percent.write(f"{chrom}\t{start}\t{end}\t{value_percent:.6f}\n")
                covered.write(f"{chrom}\t{start}\t{end}\t{covered_bp}\n")

    write_tsv(density_path, density_fields, density_rows)
    chromosome_fields = [
        "Chrom", "Chromosome_Bp", "Raw_TRF_Call_Count", "Merged_Tandem_Repeat_Locus_Count",
        "Tandem_Repeat_Union_Covered_Bp", "Tandem_Repeat_Union_Fraction",
        "Tandem_Repeat_Union_Percent",
    ]
    write_tsv(
        args.output / "FL_Africa_hap2.TRF_tandem_repeat_summary_by_chromosome.tsv",
        chromosome_fields,
        chromosome_rows,
    )

    genome_bp = sum(sizes.values())
    genome_summary = [{
        "Sample": "FL_Africa_hap2",
        "Genome_Bp": genome_bp,
        "TRF_Version": "4.09",
        "TRF_Parameters": EXPECTED_PARAMETERS,
        "Raw_TRF_Call_Count": raw_calls,
        "Merged_Tandem_Repeat_Locus_Count": total_union_loci,
        "Tandem_Repeat_Union_Covered_Bp": total_union_bp,
        "Tandem_Repeat_Union_Fraction": f"{total_union_bp / genome_bp:.8f}",
        "Tandem_Repeat_Union_Percent": f"{100 * total_union_bp / genome_bp:.6f}",
        "Minimum_TRF_Score": min_score,
        "Maximum_TRF_Score": max_score,
        "Maximum_Period": max_period,
    }]
    write_tsv(
        args.output / "FL_Africa_hap2.TRF_tandem_repeat_summary_genome.tsv",
        list(genome_summary[0]),
        genome_summary,
    )

    stats_rows = []
    for chrom in order:
        values = [float(row["Tandem_Repeat_Fraction"]) for row in density_rows if row["Chrom"] == chrom]
        stats_rows.append({
            "Chrom": chrom,
            "Window_Count": len(values),
            "Maximum_Tandem_Repeat_Fraction": f"{max(values):.8f}",
            "Mean_Tandem_Repeat_Fraction": f"{statistics.mean(values):.8f}",
            "Median_Tandem_Repeat_Fraction": f"{statistics.median(values):.8f}",
        })
    write_tsv(
        args.output / "FL_Africa_hap2.TRF_tandem_repeat_density_statistics.tsv",
        ["Chrom", "Window_Count", "Maximum_Tandem_Repeat_Fraction", "Mean_Tandem_Repeat_Fraction", "Median_Tandem_Repeat_Fraction"],
        stats_rows,
    )

    expected_windows = sum((length + args.window_size - 1) // args.window_size for length in sizes.values())
    if call_id != raw_calls:
        raise ValueError("BED call count does not match parsed TRF call count")
    if len(density_rows) != expected_windows:
        raise ValueError("Density window count mismatch")
    if sum(int(row["Tandem_Repeat_Covered_Bp"]) for row in density_rows) != total_union_bp:
        raise ValueError("Window coverage does not sum to genome-wide union coverage")
    if sum(int(row["Merged_Locus_Midpoint_Count"]) for row in density_rows) != total_union_loci:
        raise ValueError("Merged-locus midpoint counts do not sum to the union interval count")
    if sum(int(row["Raw_TRF_Call_Midpoint_Count"]) for row in density_rows) != raw_calls:
        raise ValueError("Raw-call midpoint counts do not sum to the parsed TRF call count")
    if any(int(row["Tandem_Repeat_Covered_Bp"]) > int(row["Window_Bp"]) for row in density_rows):
        raise ValueError("Tandem-repeat coverage exceeds window length")

    write_tsv(
        args.output / "parameters.tsv",
        ["Parameter", "Value", "Unit"],
        [
            {"Parameter": "Sample", "Value": "FL_Africa_hap2", "Unit": "text"},
            {"Parameter": "TRF_Version", "Value": "4.09", "Unit": "text"},
            {"Parameter": "TRF_Parameters", "Value": EXPECTED_PARAMETERS, "Unit": "text"},
            {"Parameter": "Window_Size", "Value": args.window_size, "Unit": "bp"},
            {"Parameter": "Input_Coordinate_System", "Value": "1-based_inclusive", "Unit": "text"},
            {"Parameter": "Output_Coordinate_System", "Value": "0-based_half-open", "Unit": "text"},
            {"Parameter": "Density_Definition", "Value": "merged_tandem_repeat_union_bp_divided_by_window_bp", "Unit": "text"},
            {"Parameter": "Overlap_Merge_Rule", "Value": "overlapping_or_bookended_intervals", "Unit": "text"},
            {"Parameter": "Chromosome_Name_Conversion", "Value": "chrNNB_to_chrNN", "Unit": "text"},
        ],
    )

    command = shlex.join([sys.executable, str(Path(__file__).resolve()), *sys.argv[1:]])
    (args.output / "commands.log").write_text(command + "\n", encoding="utf-8")
    (args.output / "README.md").write_text(f"""# FL Africa hap2 TRF tandem-repeat density

The source is the Tandem Repeats Finder v4.09 `.dat` output generated with parameters `{EXPECTED_PARAMETERS}`. This is tandem-repeat evidence, not transcriptome or whole-genome TE annotation. Source names `chr01B`–`chr16B` were converted to the accepted FL display names `chr01`–`chr16`.

`FL_Africa_hap2.TRF_tandem_repeat_calls.bed` contains all {raw_calls:,} TRF calls as BED4. `FL_Africa_hap2.TRF_tandem_repeat_union.bed` is the merged, nonredundant BED3 used for density. TRF 1-based inclusive coordinates were converted to 0-based half-open BED coordinates.

`FL_Africa_hap2.TRF_tandem_repeat_density.200kb.tsv` is the comprehensive 200-kb table. The recommended plotting track is `FL_Africa_hap2.TRF_tandem_repeat_fraction.200kb.bedgraph`, where density equals nonredundant tandem-repeat-covered bp divided by actual window length. All zero-density windows are explicit, and the final window of each chromosome may be shorter than {args.window_size:,} bp.

The raw calls merged into {total_union_loci:,} loci covering {total_union_bp:,} bp ({100 * total_union_bp / genome_bp:.4f}% of the FL Africa hap2 assembly).
""", encoding="utf-8")
    (args.output / "QC_Report.md").write_text(f"""# QC report

- TRF sequence headers: {sequence_headers}; all map uniquely to the 16 accepted FL chromosomes.
- TRF parameter set: `{EXPECTED_PARAMETERS}`; exact match confirmed.
- Raw calls parsed: {raw_calls:,}; BED calls written: {call_id:,}.
- Merged tandem-repeat loci: {total_union_loci:,}.
- Density windows: {len(density_rows):,}; expected windows: {expected_windows:,}.
- All intervals are within chromosome bounds and all output coordinates are 0-based half-open.
- Windows are chromosome ordered, continuous, nonoverlapping and include explicit zeros.
- Window coverage sums to the merged genome-wide coverage ({total_union_bp:,} bp).
- No window coverage exceeds its window length.
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
        } for label, path in (("TRF_DAT", args.trf_dat), ("FL_Genome_Sizes", args.genome_sizes))),
    )

    checksum_path = args.output / "output_checksums.sha256"
    output_files = sorted(path for path in args.output.iterdir() if path.is_file() and path != checksum_path)
    with checksum_path.open("w", encoding="utf-8") as handle:
        for path in output_files:
            handle.write(f"{sha256(path)}  {path.name}\n")
    for path in output_files:
        if contains_cr_or_nul(path):
            raise ValueError(f"Illegal CR or NUL character in {path}")

    print(f"OUTPUT\t{args.output.resolve()}")
    print(f"RAW_CALLS\t{raw_calls}")
    print(f"UNION_INTERVALS\t{total_union_loci}")
    print(f"TANDEM_REPEAT_COVERED_BP\t{total_union_bp}")
    print(f"TANDEM_REPEAT_PERCENT\t{100 * total_union_bp / genome_bp:.6f}")
    print(f"WINDOWS\t{len(density_rows)}")


if __name__ == "__main__":
    main()
