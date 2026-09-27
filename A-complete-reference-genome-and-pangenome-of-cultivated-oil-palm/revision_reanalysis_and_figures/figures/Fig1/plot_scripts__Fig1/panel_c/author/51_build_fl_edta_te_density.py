#!/usr/bin/env python3
"""Build 200-kb FL Africa hap2 TE density and BED files from EDTA GFF3."""

from __future__ import annotations

import argparse
import csv
import hashlib
import os
import re
import statistics
from collections import Counter, defaultdict
from pathlib import Path


WINDOW_SIZE = 200_000
CLASSES = ("LTR_Copia", "LTR_Gypsy", "LTR_Unknown", "DNA_TIR", "Helitron", "LINE", "Other")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def attributes(text: str) -> dict[str, str]:
    result = {}
    for item in text.split(";"):
        if "=" in item:
            key, value = item.split("=", 1)
            result[key] = value
    return result


def major_class(classification: str) -> str:
    if classification == "LTR/Copia":
        return "LTR_Copia"
    if classification == "LTR/Gypsy":
        return "LTR_Gypsy"
    if classification == "LTR/unknown":
        return "LTR_Unknown"
    if classification.startswith(("DNAauto/", "DNAnona/")):
        return "Helitron" if "Helitron" in classification else "DNA_TIR"
    if classification.startswith("LINE/"):
        return "LINE"
    return "Other"


def read_sizes(path: Path) -> tuple[list[str], dict[str, int]]:
    order, sizes = [], {}
    with path.open() as handle:
        for line_no, line in enumerate(handle, 1):
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 2:
                raise ValueError(f"Invalid genome-size line {line_no}")
            chrom, length_text = fields
            length = int(length_text)
            if chrom in sizes or length <= 0:
                raise ValueError(f"Invalid genome-size entry at line {line_no}")
            order.append(chrom)
            sizes[chrom] = length
    return order, sizes


def read_name_map(path: Path) -> dict[str, str]:
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        result = {
            row["Original_Name"]: row["Display_Name"]
            for row in reader
            if row["Sample"] == "FL" and row["Status"] == "Included"
        }
    if len(result) != 16:
        raise ValueError(f"Expected 16 FL chromosome mappings, found {len(result)}")
    return result


def merge_intervals(intervals: list[tuple[int, int]]) -> list[tuple[int, int]]:
    if not intervals:
        return []
    intervals.sort()
    merged = [list(intervals[0])]
    for start, end in intervals[1:]:
        if start <= merged[-1][1]:
            if end > merged[-1][1]:
                merged[-1][1] = end
        else:
            merged.append([start, end])
    return [(start, end) for start, end in merged]


def add_coverage(intervals: list[tuple[int, int]], chrom_length: int, values: list[int]) -> None:
    for start, end in intervals:
        first = start // WINDOW_SIZE
        last = (end - 1) // WINDOW_SIZE
        for index in range(first, last + 1):
            window_start = index * WINDOW_SIZE
            window_end = min(window_start + WINDOW_SIZE, chrom_length)
            values[index] += max(0, min(end, window_end) - max(start, window_start))


def write_tsv(path: Path, fields: list[str], rows) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--gff", required=True, type=Path)
    parser.add_argument("--sizes", required=True, type=Path)
    parser.add_argument("--name-map", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--script", required=True, type=Path)
    args = parser.parse_args()

    for path in (args.gff, args.sizes, args.name_map, args.script):
        if not path.is_file() or path.stat().st_size == 0:
            raise FileNotFoundError(path)
    if args.output.exists():
        raise FileExistsError(f"Refusing to overwrite: {args.output}")
    args.output.mkdir(parents=True)

    order, sizes = read_sizes(args.sizes)
    rank = {chrom: index for index, chrom in enumerate(order)}
    name_map = read_name_map(args.name_map)
    if set(name_map.values()) != set(order):
        raise ValueError("FL chromosome map does not match genome sizes")

    primary_bed = args.output / "FL_Africa_hap2.EDTA_TE.primary.bed"
    union_bed = args.output / "FL_Africa_hap2.EDTA_TE.union.bed"
    density_tsv = args.output / "FL_Africa_hap2.EDTA_TE_density.200kb.tsv"
    fraction_bg = args.output / "FL_Africa_hap2.EDTA_TE_fraction.200kb.bedgraph"
    percent_bg = args.output / "FL_Africa_hap2.EDTA_TE_percent.200kb.bedgraph"

    intervals: dict[str, list[tuple[int, int]]] = {chrom: [] for chrom in order}
    class_intervals: dict[str, dict[str, list[tuple[int, int]]]] = {
        chrom: {group: [] for group in CLASSES} for chrom in order
    }
    feature_counts: dict[str, list[int]] = {
        chrom: [0] * ((sizes[chrom] + WINDOW_SIZE - 1) // WINDOW_SIZE) for chrom in order
    }
    class_feature_counts: dict[str, dict[str, list[int]]] = {
        chrom: {group: [0] * len(feature_counts[chrom]) for group in CLASSES} for chrom in order
    }
    exact_class_counts = Counter()
    method_counts = Counter()
    primary_count = 0
    skipped_child_count = 0
    previous_key = (-1, -1)

    with args.gff.open(encoding="utf-8") as source, primary_bed.open("w", encoding="utf-8") as bed:
        for line_no, line in enumerate(source, 1):
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                raise ValueError(f"GFF field count != 9 at line {line_no}")
            seqid, _, feature, start_text, end_text, _, strand, _, attr_text = fields
            attr = attributes(attr_text)
            method = attr.get("method", "unknown")
            is_primary = method == "homology" or (method == "structural" and feature == "repeat_region")
            if not is_primary:
                skipped_child_count += 1
                continue
            if seqid not in name_map:
                raise ValueError(f"Unknown FL EDTA seqid at line {line_no}: {seqid}")
            chrom = name_map[seqid]
            start, end = int(start_text) - 1, int(end_text)
            if start < 0 or end <= start or end > sizes[chrom]:
                raise ValueError(f"Out-of-range GFF interval at line {line_no}")
            # EDTA may emit nested records with the same start in arbitrary end order.
            key = rank[chrom], start
            if key < previous_key:
                raise ValueError(f"EDTA GFF is not naturally coordinate-sorted at line {line_no}")
            previous_key = key
            classification = attr.get("classification", "Unclassified")
            group = major_class(classification)
            te_id = attr.get("ID", f"EDTA_line_{line_no}")
            midpoint = (start + end) // 2
            index = midpoint // WINDOW_SIZE

            intervals[chrom].append((start, end))
            class_intervals[chrom][group].append((start, end))
            feature_counts[chrom][index] += 1
            class_feature_counts[chrom][group][index] += 1
            exact_class_counts[classification] += 1
            method_counts[method] += 1
            primary_count += 1
            bed.write(f"{chrom}\t{start}\t{end}\t{te_id}\t{group}\t{classification}\t{method}\t{strand}\n")

    total_coverage = {}
    class_coverage = {}
    union_count = 0
    with union_bed.open("w", encoding="utf-8") as handle:
        for chrom in order:
            merged = merge_intervals(intervals[chrom])
            union_count += len(merged)
            total_coverage[chrom] = [0] * len(feature_counts[chrom])
            add_coverage(merged, sizes[chrom], total_coverage[chrom])
            for start, end in merged:
                handle.write(f"{chrom}\t{start}\t{end}\n")
            class_coverage[chrom] = {}
            for group in CLASSES:
                class_coverage[chrom][group] = [0] * len(feature_counts[chrom])
                add_coverage(
                    merge_intervals(class_intervals[chrom][group]),
                    sizes[chrom],
                    class_coverage[chrom][group],
                )

    fields = [
        "Chrom", "Start", "End", "Window_Bp", "TE_Covered_Bp", "TE_Fraction",
        "TE_Percent", "TE_Primary_Feature_Midpoint_Count",
    ]
    for group in CLASSES:
        fields.extend((f"{group}_Covered_Bp", f"{group}_Fraction", f"{group}_Percent", f"{group}_Feature_Midpoint_Count"))

    density_rows = []
    chromosome_rows = []
    with fraction_bg.open("w", encoding="utf-8") as frac, percent_bg.open("w", encoding="utf-8") as pct:
        for chrom in order:
            for index, covered_bp in enumerate(total_coverage[chrom]):
                start = index * WINDOW_SIZE
                end = min(start + WINDOW_SIZE, sizes[chrom])
                window_bp = end - start
                row = {
                    "Chrom": chrom, "Start": start, "End": end, "Window_Bp": window_bp,
                    "TE_Covered_Bp": covered_bp, "TE_Fraction": f"{covered_bp / window_bp:.8f}",
                    "TE_Percent": f"{100 * covered_bp / window_bp:.6f}",
                    "TE_Primary_Feature_Midpoint_Count": feature_counts[chrom][index],
                }
                for group in CLASSES:
                    class_bp = class_coverage[chrom][group][index]
                    row.update({
                        f"{group}_Covered_Bp": class_bp,
                        f"{group}_Fraction": f"{class_bp / window_bp:.8f}",
                        f"{group}_Percent": f"{100 * class_bp / window_bp:.6f}",
                        f"{group}_Feature_Midpoint_Count": class_feature_counts[chrom][group][index],
                    })
                density_rows.append(row)
                frac.write(f"{chrom}\t{start}\t{end}\t{covered_bp / window_bp:.8f}\n")
                pct.write(f"{chrom}\t{start}\t{end}\t{100 * covered_bp / window_bp:.6f}\n")

            chrom_covered = sum(total_coverage[chrom])
            chromosome_rows.append({
                "Chrom": chrom,
                "Chromosome_Bp": sizes[chrom],
                "TE_Covered_Bp": chrom_covered,
                "TE_Fraction": f"{chrom_covered / sizes[chrom]:.8f}",
                "TE_Percent": f"{100 * chrom_covered / sizes[chrom]:.6f}",
                "TE_Primary_Feature_Count": sum(feature_counts[chrom]),
            })

    write_tsv(density_tsv, fields, density_rows)
    write_tsv(
        args.output / "FL_Africa_hap2.EDTA_TE_summary_by_chromosome.tsv",
        ["Chrom", "Chromosome_Bp", "TE_Covered_Bp", "TE_Fraction", "TE_Percent", "TE_Primary_Feature_Count"],
        chromosome_rows,
    )

    genome_bp = sum(sizes.values())
    genome_covered = sum(sum(values) for values in total_coverage.values())
    summary = [{
        "Sample": "FL_Africa_hap2", "Genome_Bp": genome_bp, "TE_Covered_Bp": genome_covered,
        "TE_Fraction": f"{genome_covered / genome_bp:.8f}",
        "TE_Percent": f"{100 * genome_covered / genome_bp:.6f}",
        "TE_Primary_Feature_Count": primary_count, "TE_Union_Interval_Count": union_count,
        "Skipped_Structural_Child_Feature_Count": skipped_child_count,
    }]
    write_tsv(
        args.output / "FL_Africa_hap2.EDTA_TE_summary_genome.tsv",
        list(summary[0]), summary,
    )
    write_tsv(
        args.output / "FL_Africa_hap2.EDTA_TE_class_counts.tsv",
        ["Classification", "Major_Class", "Primary_Feature_Count"],
        ({"Classification": key, "Major_Class": major_class(key), "Primary_Feature_Count": value}
         for key, value in sorted(exact_class_counts.items())),
    )

    if len(density_rows) != sum((length + WINDOW_SIZE - 1) // WINDOW_SIZE for length in sizes.values()):
        raise ValueError("Density row count mismatch")
    if sum(row["TE_Primary_Feature_Midpoint_Count"] for row in density_rows) != primary_count:
        raise ValueError("Primary feature midpoint counts do not sum to primary features")
    if any(int(row["TE_Covered_Bp"]) > int(row["Window_Bp"]) for row in density_rows):
        raise ValueError("Total TE coverage exceeds window length")

    stats_rows = []
    for chrom in order:
        values = [total_coverage[chrom][i] / min(WINDOW_SIZE, sizes[chrom] - i * WINDOW_SIZE)
                  for i in range(len(total_coverage[chrom]))]
        stats_rows.append({
            "Chrom": chrom, "Window_Count": len(values), "Maximum_TE_Fraction": f"{max(values):.8f}",
            "Mean_TE_Fraction": f"{statistics.mean(values):.8f}",
            "Median_TE_Fraction": f"{statistics.median(values):.8f}",
        })
    write_tsv(args.output / "FL_Africa_hap2.EDTA_TE_density_statistics.tsv",
              ["Chrom", "Window_Count", "Maximum_TE_Fraction", "Mean_TE_Fraction", "Median_TE_Fraction"], stats_rows)

    parameters = [
        ("Sample", "FL_Africa_hap2", "text"),
        ("EDTA_Version", "2.2.2", "text"),
        ("Window_Size", str(WINDOW_SIZE), "bp"),
        ("Coordinate_System", "0-based_half-open", "text"),
        ("Total_Density", "union_bp_per_window_divided_by_window_bp", "text"),
        ("Structural_Record_Rule", "method=structural AND feature=repeat_region", "text"),
        ("Homology_Record_Rule", "method=homology", "text"),
        ("Feature_Count_Assignment", "primary_feature_midpoint", "text"),
    ]
    write_tsv(args.output / "parameters.tsv", ["Parameter", "Value", "Unit"],
              ({"Parameter": key, "Value": value, "Unit": unit} for key, value, unit in parameters))

    command = (
        f"python3 {args.script.resolve()} --gff {args.gff.resolve()} --sizes {args.sizes.resolve()} "
        f"--name-map {args.name_map.resolve()} --output {args.output.resolve()} --script {args.script.resolve()}"
    )
    (args.output / "commands.log").write_text(command + "\n", encoding="utf-8")
    (args.output / "README.md").write_text(f"""# FL Africa hap2 EDTA TE density

The source is the final EDTA v2.2.2 `TEanno.gff3`. Chromosome names `chr01B`–`chr16B` were converted with the accepted SyntenyViz chromosome map to `chr01`–`chr16`.

`FL_Africa_hap2.EDTA_TE_density.200kb.tsv` is the comprehensive 200-kb density table. Coordinates are 0-based half-open; the final window may be shorter. Total and class coverage are independently merged interval unions. Because different TE classes can be nested, class coverage sums may exceed total union coverage.

`FL_Africa_hap2.EDTA_TE.primary.bed` contains one accepted primary EDTA record per line: `Chrom Start End TE_ID Major_Class Classification Method Strand`. Structural LTR child features such as LTRs and target-site duplications are excluded; their parent `repeat_region` is retained. `FL_Africa_hap2.EDTA_TE.union.bed` is the three-column merged union used for total density.

Genome-wide total TE union coverage is {genome_covered:,} bp ({100 * genome_covered / genome_bp:.4f}%). The analysis retained {primary_count:,} primary EDTA records and skipped {skipped_child_count:,} structural child records.
""", encoding="utf-8")
    (args.output / "QC_Report.md").write_text(f"""# QC report

- GFF seqids: all 16 mapped to accepted FL display chromosomes.
- Primary EDTA records: {primary_count:,}.
- Structural child records excluded from element counting: {skipped_child_count:,}.
- Merged total TE intervals: {union_count:,}.
- Density windows: {len(density_rows):,}; all authoritative 200-kb windows present, including zero values.
- Primary feature midpoint sum: {sum(row['TE_Primary_Feature_Midpoint_Count'] for row in density_rows):,}.
- Coordinates: all within chromosome bounds; BED and bedGraph are 0-based half-open.
- Total coverage: union-based and never exceeds window length.
- Status: PASS.
""", encoding="utf-8")

    input_rows = []
    for label, path in (("EDTA_GFF3", args.gff), ("FL_Genome_Sizes", args.sizes), ("Chromosome_Name_Map", args.name_map)):
        input_rows.append({"Input_Type": label, "Size_Bytes": path.stat().st_size, "SHA256": sha256(path), "Absolute_Path": str(path.resolve())})
    write_tsv(args.output / "input_checksums.tsv", ["Input_Type", "Size_Bytes", "SHA256", "Absolute_Path"], input_rows)

    output_checksum = args.output / "output_checksums.sha256"
    files = sorted(path for path in args.output.iterdir() if path.is_file() and path != output_checksum)
    with output_checksum.open("w", encoding="utf-8") as handle:
        for path in files:
            handle.write(f"{sha256(path)}  {path.name}\n")

    for path in files:
        data = path.read_bytes()
        if b"\r" in data or b"\x00" in data:
            raise ValueError(f"Illegal CR or NUL character: {path}")
    print(f"OUTPUT\t{args.output.resolve()}")
    print(f"PRIMARY_FEATURES\t{primary_count}")
    print(f"UNION_INTERVALS\t{union_count}")
    print(f"TE_COVERED_BP\t{genome_covered}")
    print(f"TE_PERCENT\t{100 * genome_covered / genome_bp:.6f}")
    print(f"WINDOWS\t{len(density_rows)}")


if __name__ == "__main__":
    main()
