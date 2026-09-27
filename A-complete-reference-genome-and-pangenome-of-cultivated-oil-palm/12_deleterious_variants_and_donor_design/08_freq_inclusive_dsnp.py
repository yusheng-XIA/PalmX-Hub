#!/usr/bin/env python3
"""Prepare frequency-inclusive Phoenix ALT-derived conserved SNPs for Stage 09.

The source population catalog is streamed. Only Filter=PASS and
Conserved_Region_Proxy=Yes sites are retained. Existing same-reference Phoenix
PAF cs tags are reused; no alignment or GWAS input is used.
"""

import argparse
import bisect
import csv
import os
import re
from collections import Counter, defaultdict
from pathlib import Path


CS_TOKEN_RE = re.compile(
    r"(:\d+|=[A-Za-z]+|\*[A-Za-z][A-Za-z]|\+[A-Za-z]+|-[A-Za-z]+|~[A-Za-z]{2}\d+[A-Za-z]{2})"
)
VALID_BASES = {"A", "C", "G", "T"}
YES_VALUES = {"yes", "true", "1"}
NA_VALUES = {"", "na", "nan", "none"}


def parse_args():
    parser = argparse.ArgumentParser(
        description="Extract all-frequency conserved SNPs and polarize them with Phoenix."
    )
    parser.add_argument("--catalog", required=True)
    parser.add_argument("--paf", required=True)
    parser.add_argument("--sample-manifest", required=True)
    parser.add_argument("--reference-fai", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--min-mapq", type=int, default=20)
    parser.add_argument("--expected-catalog-rows", type=int, default=14_409_445)
    parser.add_argument("--expected-all38-samples", type=int, default=38)
    parser.add_argument("--expected-african35-samples", type=int, default=35)
    return parser.parse_args()


def parse_samples(value):
    if value is None:
        return []
    samples = []
    for token in re.split(r"[;,|]", str(value)):
        sample = token.strip()
        if sample and sample.lower() not in NA_VALUES:
            samples.append(sample)
    return sorted(set(samples))


def as_int(value, default=None):
    try:
        if value is None or str(value).strip().lower() in NA_VALUES:
            return default
        return int(float(str(value)))
    except (TypeError, ValueError):
        return default


def read_manifest(path):
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = reader.fieldnames or []
        required = ["Sample_ID", "Include_All38", "Include_African35", "Status"]
        missing = [field for field in required if field not in fields]
        if missing:
            raise ValueError(f"Manifest missing columns: {','.join(missing)}")
        rows = list(reader)
    if len(rows) != len({row["Sample_ID"] for row in rows}):
        raise ValueError("Duplicate Sample_ID in manifest")
    all38 = []
    african35 = []
    for row in rows:
        sample = row["Sample_ID"].strip()
        if row.get("Status", "").strip().upper() != "PASS":
            raise ValueError(f"Manifest status is not PASS: {sample}")
        if row["Include_All38"].strip().lower() in YES_VALUES:
            all38.append(sample)
        if row["Include_African35"].strip().lower() in YES_VALUES:
            african35.append(sample)
    return sorted(all38), sorted(african35)


def read_fai(path):
    chrom_lengths = {}
    with open(path) as handle:
        for line in handle:
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            length = as_int(fields[1] if len(fields) > 1 else None)
            if length is None or length <= 0:
                raise ValueError(f"Invalid FAI line: {line.rstrip()}")
            chrom_lengths[fields[0]] = length
    expected = {f"chr{i:02d}B" for i in range(1, 17)}
    if set(chrom_lengths) != expected:
        raise ValueError(
            "Reference chromosome set mismatch: " + ",".join(sorted(chrom_lengths))
        )
    return chrom_lengths


def write_tsv(path, rows, fields):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        raise FileExistsError(f"Refusing to overwrite existing output: {path}")
    tmp = path.with_name(path.name + f".tmp.{os.getpid()}")
    with open(tmp, "w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fields, delimiter="\t", lineterminator="\n", extrasaction="ignore"
        )
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "NA") for field in fields})
    os.replace(tmp, path)


def open_temp_writer(path, fields):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        raise FileExistsError(f"Refusing to overwrite existing output: {path}")
    tmp = path.with_name(path.name + f".tmp.{os.getpid()}")
    handle = open(tmp, "w", newline="")
    writer = csv.DictWriter(
        handle, fieldnames=fields, delimiter="\t", lineterminator="\n", extrasaction="ignore"
    )
    writer.writeheader()
    return path, tmp, handle, writer


def close_temp_writer(path, tmp, handle):
    handle.flush()
    os.fsync(handle.fileno())
    handle.close()
    os.replace(tmp, path)


def stream_catalog(catalog, output_path, all_samples, chrom_lengths):
    compact_fields = [
        "SNP_ID", "Chrom", "Pos", "Ref", "Alt", "Sample_Count", "Frequency", "Samples",
        "Source_Filter", "Conserved_Region_Proxy",
    ]
    path, tmp, output, writer = open_temp_writer(output_path, compact_fields)
    rows = []
    stats = Counter()
    examples = defaultdict(list)
    with open(catalog, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = reader.fieldnames or []
        required = [
            "SNP_ID", "Chrom", "Pos", "Ref", "Alt", "Sample_Count", "Frequency", "Samples",
            "Filter", "Conserved_Region_Proxy",
        ]
        missing = [field for field in required if field not in fields]
        if missing:
            raise ValueError(f"Catalog missing columns: {','.join(missing)}")
        for row_number, row in enumerate(reader, 2):
            stats["Catalog_Rows"] += 1
            conserved = row.get("Conserved_Region_Proxy", "").strip().lower() in YES_VALUES
            if not conserved:
                continue
            stats["Conserved_Rows_All_Filters"] += 1
            if row.get("Filter", "").strip() != "PASS":
                stats["Excluded_Conserved_NonPASS"] += 1
                continue
            chrom = row.get("Chrom", "")
            pos = as_int(row.get("Pos"))
            ref = row.get("Ref", "").upper()
            alt = row.get("Alt", "").upper()
            carriers = parse_samples(row.get("Samples"))
            sample_count = as_int(row.get("Sample_Count"))
            valid = True
            if chrom not in chrom_lengths or pos is None or not (1 <= pos <= chrom_lengths.get(chrom, 0)):
                stats["Invalid_Coordinate"] += 1
                if len(examples["Invalid_Coordinate"]) < 10:
                    examples["Invalid_Coordinate"].append(f"row{row_number}:{chrom}:{pos}")
                valid = False
            if ref not in VALID_BASES or alt not in VALID_BASES or ref == alt:
                stats["Invalid_SNV_Alleles"] += 1
                if len(examples["Invalid_SNV_Alleles"]) < 10:
                    examples["Invalid_SNV_Alleles"].append(f"row{row_number}:{ref}>{alt}")
                valid = False
            unknown = sorted(set(carriers) - all_samples)
            if unknown:
                stats["Unknown_Carriers"] += 1
                if len(examples["Unknown_Carriers"]) < 10:
                    examples["Unknown_Carriers"].append(f"row{row_number}:{','.join(unknown)}")
                valid = False
            if sample_count is None or sample_count != len(carriers) or sample_count <= 0:
                stats["Carrier_Count_Mismatch"] += 1
                if len(examples["Carrier_Count_Mismatch"]) < 10:
                    examples["Carrier_Count_Mismatch"].append(
                        f"row{row_number}:declared={sample_count}:parsed={len(carriers)}"
                    )
                valid = False
            if not valid:
                continue
            compact = {
                "SNP_ID": row["SNP_ID"],
                "Chrom": chrom,
                "Pos": pos,
                "Ref": ref,
                "Alt": alt,
                "Sample_Count": sample_count,
                "Frequency": row["Frequency"],
                "Samples": ";".join(carriers),
                "Source_Filter": row["Filter"],
                "Conserved_Region_Proxy": "Yes",
            }
            writer.writerow(compact)
            rows.append(
                (
                    compact["SNP_ID"], chrom, pos, ref, alt, sample_count,
                    compact["Frequency"], compact["Samples"], tuple(carriers),
                )
            )
            stats["Retained_PASS_Conserved"] += 1
    errors = (
        stats["Invalid_Coordinate"]
        + stats["Invalid_SNV_Alleles"]
        + stats["Unknown_Carriers"]
        + stats["Carrier_Count_Mismatch"]
    )
    if errors:
        output.close()
        detail = " | ".join(f"{key}={';'.join(value)}" for key, value in examples.items())
        raise ValueError(f"Retained conserved-site validation failed ({errors} rows): {detail}")
    close_temp_writer(path, tmp, output)
    return rows, stats


def get_cs(fields):
    for field in fields[12:]:
        if field.startswith("cs:Z:"):
            return field[5:]
    return None


def build_site_index(rows):
    positions_by_chrom = defaultdict(list)
    index_by_site = defaultdict(list)
    refs = []
    for index, row in enumerate(rows):
        chrom, pos, ref = row[1], row[2], row[3]
        positions_by_chrom[chrom].append(pos)
        index_by_site[(chrom, pos)].append(index)
        refs.append(ref)
    for chrom in positions_by_chrom:
        positions_by_chrom[chrom] = sorted(set(positions_by_chrom[chrom]))
    return positions_by_chrom, index_by_site, refs


def add_call(calls, index_by_site, chrom, pos, base, mapq):
    for index in index_by_site.get((chrom, pos), []):
        calls[index].append((base, mapq))


def add_match_range(
    calls, positions_by_chrom, index_by_site, refs, chrom, start0, end0, mapq
):
    positions = positions_by_chrom.get(chrom)
    if not positions or end0 <= start0:
        return
    left = bisect.bisect_left(positions, start0 + 1)
    right = bisect.bisect_right(positions, end0)
    for pos in positions[left:right]:
        for index in index_by_site.get((chrom, pos), []):
            calls[index].append((refs[index], mapq))


def parse_paf(paf_path, rows, min_mapq):
    positions_by_chrom, index_by_site, refs = build_site_index(rows)
    calls = defaultdict(list)
    stats = Counter()
    with open(paf_path) as handle:
        for line_number, line in enumerate(handle, 1):
            if not line.strip():
                continue
            stats["PAF_Rows"] += 1
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 12:
                stats["Malformed_PAF"] += 1
                continue
            mapq = as_int(fields[11], -1)
            if mapq < min_mapq:
                stats["Skipped_Low_MAPQ"] += 1
                continue
            chrom = fields[5]
            positions = positions_by_chrom.get(chrom)
            if not positions:
                stats["Skipped_Target_Not_In_Conserved_SNPs"] += 1
                continue
            target_start = as_int(fields[7])
            target_end = as_int(fields[8])
            if target_start is None or target_end is None or target_end <= target_start:
                stats["Malformed_Target_Coordinate"] += 1
                continue
            left = bisect.bisect_left(positions, target_start + 1)
            right = bisect.bisect_right(positions, target_end)
            if left >= right:
                stats["Skipped_No_Target_SNPs"] += 1
                continue
            cs = get_cs(fields)
            if not cs:
                stats["Missing_CS_At_Candidate_Overlap"] += 1
                continue
            stats["Used_Alignments"] += 1
            target0 = target_start
            for match in CS_TOKEN_RE.finditer(cs):
                token = match.group(0)
                operation = token[0]
                if operation == ":":
                    length = int(token[1:])
                    add_match_range(
                        calls, positions_by_chrom, index_by_site, refs, chrom,
                        target0, target0 + length, mapq,
                    )
                    target0 += length
                elif operation == "=":
                    length = len(token) - 1
                    add_match_range(
                        calls, positions_by_chrom, index_by_site, refs, chrom,
                        target0, target0 + length, mapq,
                    )
                    target0 += length
                elif operation == "*":
                    query_base = token[2].upper()
                    add_call(calls, index_by_site, chrom, target0 + 1, query_base, mapq)
                    target0 += 1
                elif operation == "-":
                    length = len(token) - 1
                    gap_left = bisect.bisect_left(positions, target0 + 1)
                    gap_right = bisect.bisect_right(positions, target0 + length)
                    for pos in positions[gap_left:gap_right]:
                        add_call(calls, index_by_site, chrom, pos, "-", mapq)
                    target0 += length
                elif operation == "+":
                    continue
                elif operation == "~":
                    intron = re.match(r"~[A-Za-z]{2}(\d+)[A-Za-z]{2}", token)
                    if intron:
                        length = int(intron.group(1))
                        gap_left = bisect.bisect_left(positions, target0 + 1)
                        gap_right = bisect.bisect_right(positions, target0 + length)
                        for pos in positions[gap_left:gap_right]:
                            add_call(calls, index_by_site, chrom, pos, "N", mapq)
                        target0 += length
    return calls, stats


def summarize_call(row, site_calls):
    if not site_calls:
        return {
            "Phoenix_Base": "NA",
            "Phoenix_Base_Status": "No_Phoenix_Coverage",
            "Phoenix_Alignment_Count": 0,
            "Phoenix_MapQ_Max": "NA",
            "ALT_Polarity_Phoenix": "Unknown_Polarity",
            "ALT_Polarity_Phoenix_Evidence": "No_Phoenix_Alignment_Covers_SNP",
        }
    base_counts = Counter(base.upper() for base, _ in site_calls)
    mapq_max = max(mapq for _, mapq in site_calls)
    if len(base_counts) > 1:
        return {
            "Phoenix_Base": ";".join(f"{base}:{count}" for base, count in sorted(base_counts.items())),
            "Phoenix_Base_Status": "Ambiguous_Multiple_Bases",
            "Phoenix_Alignment_Count": len(site_calls),
            "Phoenix_MapQ_Max": mapq_max,
            "ALT_Polarity_Phoenix": "Unknown_Polarity",
            "ALT_Polarity_Phoenix_Evidence": "Multiple_Phoenix_Bases_At_Site",
        }
    base = next(iter(base_counts))
    status = "Single_Alignment" if len(site_calls) == 1 else "Multiple_Alignments_Same_Base"
    ref = row[3]
    alt = row[4]
    if base == ref:
        polarity = "ALT_Derived"
        evidence = "Phoenix_Base_Equals_Reference"
    elif base == alt:
        polarity = "Ref_Derived_or_Ancestral_ALT"
        evidence = "Phoenix_Base_Equals_ALT"
    elif base in VALID_BASES:
        polarity = "Other_Outgroup_Base"
        evidence = "Phoenix_Base_Is_Third_Allele"
    elif base == "-":
        polarity = "Unknown_Polarity"
        evidence = "Phoenix_Gap_At_SNP"
    else:
        polarity = "Unknown_Polarity"
        evidence = "Phoenix_Base_Not_ACGT"
    return {
        "Phoenix_Base": base,
        "Phoenix_Base_Status": status,
        "Phoenix_Alignment_Count": len(site_calls),
        "Phoenix_MapQ_Max": mapq_max,
        "ALT_Polarity_Phoenix": polarity,
        "ALT_Polarity_Phoenix_Evidence": evidence,
    }


def emit_polarity(rows, calls, polarity_dir, all38, african35):
    fields = [
        "SNP_ID", "Chrom", "Pos", "Ref", "Alt", "Sample_Count", "Frequency", "Samples",
        "Phoenix_Base", "Phoenix_Base_Status", "Phoenix_Alignment_Count", "Phoenix_MapQ_Max",
        "ALT_Polarity_Phoenix", "ALT_Polarity_Phoenix_Evidence",
    ]
    all_path, all_tmp, all_handle, all_writer = open_temp_writer(
        polarity_dir / "conserved_snp.phoenix_polarized.tsv", fields
    )
    der_path, der_tmp, der_handle, der_writer = open_temp_writer(
        polarity_dir / "alt_derived_conserved_snp.tsv", fields
    )
    polarity_counts = Counter()
    status_counts = Counter()
    sample_burden = Counter()
    panel_counts = {
        "All38": Counter(),
        "African35": Counter(),
    }
    panel_sets = {"All38": set(all38), "African35": set(african35)}
    derived_count = 0
    for index, row in enumerate(rows):
        summary = summarize_call(row, calls.get(index, []))
        output = {
            "SNP_ID": row[0], "Chrom": row[1], "Pos": row[2], "Ref": row[3], "Alt": row[4],
            "Sample_Count": row[5], "Frequency": row[6], "Samples": row[7], **summary,
        }
        all_writer.writerow(output)
        polarity_counts[summary["ALT_Polarity_Phoenix"]] += 1
        status_counts[summary["Phoenix_Base_Status"]] += 1
        if summary["ALT_Polarity_Phoenix"] != "ALT_Derived":
            continue
        derived_count += 1
        der_writer.writerow(output)
        carriers = set(row[8])
        for sample in carriers:
            sample_burden[sample] += 1
        for panel, panel_set in panel_sets.items():
            selected = carriers & panel_set
            if not selected:
                panel_counts[panel]["No_Panel_Carrier"] += 1
            else:
                panel_counts[panel]["Derived_Site_Count"] += 1
                if selected == panel_set:
                    panel_counts[panel]["Fixed_Unavoidable_Site_Count"] += 1
                else:
                    panel_counts[panel]["Variable_Site_Count"] += 1
    close_temp_writer(all_path, all_tmp, all_handle)
    close_temp_writer(der_path, der_tmp, der_handle)
    return (
        derived_count, polarity_counts, status_counts, sample_burden, panel_counts
    )


def add_check(checks, check, passed, observed, expected, notes):
    checks.append(
        {
            "Check": check,
            "Status": "PASS" if passed else "FAIL",
            "Observed": observed,
            "Expected": expected,
            "Notes": notes,
        }
    )


def main():
    args = parse_args()
    outroot = Path(args.outdir)
    if outroot.exists() and any(outroot.iterdir()):
        raise FileExistsError(f"Output directory is non-empty; refusing to overwrite: {outroot}")
    outroot.mkdir(parents=True, exist_ok=True)
    data_dir = outroot / "data"
    polarity_dir = data_dir / "phoenix_polarity_conserved"
    qa_dir = polarity_dir / "qa"
    qa_dir.mkdir(parents=True, exist_ok=False)

    all38, african35 = read_manifest(args.sample_manifest)
    if len(all38) != args.expected_all38_samples:
        raise ValueError(f"All38 sample count {len(all38)} != {args.expected_all38_samples}")
    if len(african35) != args.expected_african35_samples:
        raise ValueError(f"African35 sample count {len(african35)} != {args.expected_african35_samples}")
    chrom_lengths = read_fai(args.reference_fai)

    conserved_path = data_dir / "conserved_snp_allfreq.tsv"
    rows, catalog_stats = stream_catalog(
        args.catalog, conserved_path, set(all38), chrom_lengths
    )
    if not rows:
        raise RuntimeError("No PASS conserved SNPs retained")
    calls, paf_stats = parse_paf(args.paf, rows, args.min_mapq)
    (
        derived_count, polarity_counts, status_counts, sample_burden, panel_counts
    ) = emit_polarity(rows, calls, polarity_dir, all38, african35)

    summary_rows = []
    for key, value in sorted(catalog_stats.items()):
        summary_rows.append(
            {"Summary_Type": "Catalog_Filter", "Group": key, "Count": value, "Notes": "Stage05 population catalog"}
        )
    for key, value in sorted(paf_stats.items()):
        summary_rows.append(
            {"Summary_Type": "PAF_Parse", "Group": key, "Count": value, "Notes": f"Min_MAPQ={args.min_mapq}"}
        )
    for key, value in sorted(polarity_counts.items()):
        summary_rows.append(
            {"Summary_Type": "Polarity", "Group": key, "Count": value, "Notes": "PASS conserved all-frequency SNPs"}
        )
    for key, value in sorted(status_counts.items()):
        summary_rows.append(
            {"Summary_Type": "Phoenix_Base_Status", "Group": key, "Count": value, "Notes": "PASS conserved all-frequency SNPs"}
        )
    write_tsv(
        polarity_dir / "polarity_summary.tsv",
        summary_rows,
        ["Summary_Type", "Group", "Count", "Notes"],
    )
    burden_rows = [
        {"Sample_ID": sample, "ALT_Derived_Conserved_SNP_Count": sample_burden[sample]}
        for sample in all38
    ]
    write_tsv(
        polarity_dir / "derived_sample_burden.tsv",
        burden_rows,
        ["Sample_ID", "ALT_Derived_Conserved_SNP_Count"],
    )
    panel_rows = []
    for panel, sample_count in (("All38", len(all38)), ("African35", len(african35))):
        counts = panel_counts[panel]
        panel_rows.append(
            {
                "Panel": panel,
                "Sample_Count": sample_count,
                "Derived_Site_Count": counts["Derived_Site_Count"],
                "Variable_Site_Count": counts["Variable_Site_Count"],
                "Fixed_Unavoidable_Site_Count": counts["Fixed_Unavoidable_Site_Count"],
                "No_Panel_Carrier_Site_Count": counts["No_Panel_Carrier"],
            }
        )
    write_tsv(
        polarity_dir / "panel_derived_summary.tsv",
        panel_rows,
        [
            "Panel", "Sample_Count", "Derived_Site_Count", "Variable_Site_Count",
            "Fixed_Unavoidable_Site_Count", "No_Panel_Carrier_Site_Count",
        ],
    )

    checks = []
    add_check(
        checks, "Catalog_row_count", catalog_stats["Catalog_Rows"] == args.expected_catalog_rows,
        catalog_stats["Catalog_Rows"], args.expected_catalog_rows, "Stage05 global catalog count",
    )
    add_check(
        checks, "PASS_conserved_nonempty", len(rows) > 0, len(rows), ">0",
        "Filter=PASS and Conserved_Region_Proxy=Yes",
    )
    add_check(
        checks, "Retained_site_validation", all(
            catalog_stats[key] == 0
            for key in ("Invalid_Coordinate", "Invalid_SNV_Alleles", "Unknown_Carriers", "Carrier_Count_Mismatch")
        ),
        sum(catalog_stats[key] for key in ("Invalid_Coordinate", "Invalid_SNV_Alleles", "Unknown_Carriers", "Carrier_Count_Mismatch")),
        0, "coordinates, alleles and carriers",
    )
    add_check(
        checks, "PAF_missing_cs_at_candidate_overlap",
        paf_stats["Missing_CS_At_Candidate_Overlap"] == 0,
        paf_stats["Missing_CS_At_Candidate_Overlap"], 0, "used candidate-overlapping alignments",
    )
    add_check(
        checks, "PAF_used_alignments_nonempty", paf_stats["Used_Alignments"] > 0,
        paf_stats["Used_Alignments"], ">0", f"Min_MAPQ={args.min_mapq}",
    )
    add_check(
        checks, "Polarity_category_sum", sum(polarity_counts.values()) == len(rows),
        sum(polarity_counts.values()), len(rows), "all retained conserved SNPs classified",
    )
    add_check(
        checks, "ALT_derived_nonempty", derived_count > 0, derived_count, ">0",
        "strict Phoenix base equals Africa_hap2 reference",
    )
    add_check(
        checks, "All38_sample_burden_complete", len(burden_rows) == args.expected_all38_samples,
        len(burden_rows), args.expected_all38_samples, "zero counts retained if present",
    )
    add_check(
        checks, "All38_panel_derived_count", panel_counts["All38"]["Derived_Site_Count"] == derived_count,
        panel_counts["All38"]["Derived_Site_Count"], derived_count, "all source carriers belong to All38",
    )
    add_check(
        checks, "African35_panel_derived_nonempty", panel_counts["African35"]["Derived_Site_Count"] > 0,
        panel_counts["African35"]["Derived_Site_Count"], ">0", "at least one African35 carrier",
    )
    write_tsv(
        qa_dir / "final_integrity_summary.tsv",
        checks,
        ["Check", "Status", "Observed", "Expected", "Notes"],
    )
    if any(row["Status"] == "FAIL" for row in checks):
        raise RuntimeError("Frequency-inclusive dSNP preparation integrity failure")

    print(
        "[OK] Frequency-inclusive Phoenix polarization complete | "
        f"Catalog={catalog_stats['Catalog_Rows']} | PASS_conserved={len(rows)} | "
        f"ALT_Derived={derived_count} | Fixed_All38={panel_counts['All38']['Fixed_Unavoidable_Site_Count']} | "
        f"Fixed_African35={panel_counts['African35']['Fixed_Unavoidable_Site_Count']}"
    )


if __name__ == "__main__":
    main()
