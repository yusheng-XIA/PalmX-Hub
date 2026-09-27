#!/usr/bin/env python3
import argparse
import bisect
import csv
import re
from collections import Counter, defaultdict
from pathlib import Path


BASE = Path("dsv_analysis")
DEFAULT_V1_DIR = BASE / "results_hap38/05_dSNP_minimap_hap38"
CS_TOKEN_RE = re.compile(r"(:\d+|=[A-Za-z]+|\*[A-Za-z][A-Za-z]|\+[A-Za-z]+|-[A-Za-z]+|~[A-Za-z]{2}\d+[A-Za-z]{2})")
VALID_BASES = {"A", "C", "G", "T"}


def parse_args():
    parser = argparse.ArgumentParser(description="Reproduce Phoenix polarization for the hap38 legacy-compatible rare functional dSNP set.")
    parser.add_argument("--rare-functional", default=str(DEFAULT_V1_DIR / "rare_functional_snp.tsv"))
    parser.add_argument("--paf", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--min-mapq", type=int, default=20)
    return parser.parse_args()


def write_tsv(path, rows, fieldnames):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "NA") for field in fieldnames})


def load_snps(path):
    rows = []
    positions_by_chrom = defaultdict(list)
    index_by_site = defaultdict(list)
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fieldnames = list(reader.fieldnames or [])
        for idx, row in enumerate(reader):
            chrom = row["Chrom"]
            pos = int(row["Pos"])
            rows.append(row)
            positions_by_chrom[chrom].append(pos)
            index_by_site[(chrom, pos)].append(idx)
    for chrom in positions_by_chrom:
        positions_by_chrom[chrom] = sorted(set(positions_by_chrom[chrom]))
    return rows, fieldnames, positions_by_chrom, index_by_site


def get_cs(fields):
    for field in fields[12:]:
        if field.startswith("cs:Z:"):
            return field[5:]
    return None


def add_call(calls, index_by_site, chrom, pos, base, state, mapq, aln_id):
    for idx in index_by_site.get((chrom, pos), []):
        calls[idx].append((base, state, mapq, aln_id))


def add_match_range(calls, positions_by_chrom, index_by_site, chrom, start0, end0, mapq, aln_id):
    positions = positions_by_chrom.get(chrom)
    if not positions or end0 <= start0:
        return
    left = bisect.bisect_left(positions, start0 + 1)
    right = bisect.bisect_right(positions, end0)
    for pos in positions[left:right]:
        for idx in index_by_site.get((chrom, pos), []):
            ref_base = rows_ref_cache[idx]
            calls[idx].append((ref_base, "Reference_Match", mapq, aln_id))


rows_ref_cache = []


def parse_paf_calls(paf_path, rows, positions_by_chrom, index_by_site, min_mapq):
    calls = defaultdict(list)
    stats = Counter()
    global rows_ref_cache
    rows_ref_cache = [row["Ref"].upper() for row in rows]
    with open(paf_path) as handle:
        for aln_id, line in enumerate(handle, start=1):
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 12:
                stats["Malformed_PAF"] += 1
                continue
            mapq = int(fields[11])
            if mapq < min_mapq:
                stats["Skipped_Low_MAPQ"] += 1
                continue
            chrom = fields[5]
            if chrom not in positions_by_chrom:
                stats["Skipped_Target_Not_In_SNPs"] += 1
                continue
            tstart = int(fields[7])
            tend = int(fields[8])
            positions = positions_by_chrom[chrom]
            if not positions:
                continue
            if bisect.bisect_left(positions, tstart + 1) == bisect.bisect_left(positions, tend + 1):
                stats["Skipped_No_Target_SNPs"] += 1
                continue
            cs = get_cs(fields)
            if not cs:
                stats["Missing_CS"] += 1
                continue
            stats["Used_Alignments"] += 1
            target0 = tstart
            for token in CS_TOKEN_RE.findall(cs):
                op = token[0]
                if op == ":":
                    length = int(token[1:])
                    add_match_range(calls, positions_by_chrom, index_by_site, chrom, target0, target0 + length, mapq, aln_id)
                    target0 += length
                elif op == "=":
                    seq = token[1:].upper()
                    add_match_range(calls, positions_by_chrom, index_by_site, chrom, target0, target0 + len(seq), mapq, aln_id)
                    target0 += len(seq)
                elif op == "*":
                    ref_base = token[1].upper()
                    query_base = token[2].upper()
                    add_call(calls, index_by_site, chrom, target0 + 1, query_base, f"Substitution_ref_{ref_base}", mapq, aln_id)
                    target0 += 1
                elif op == "-":
                    seq = token[1:].upper()
                    positions = positions_by_chrom.get(chrom, [])
                    left = bisect.bisect_left(positions, target0 + 1)
                    right = bisect.bisect_right(positions, target0 + len(seq))
                    for pos in positions[left:right]:
                        add_call(calls, index_by_site, chrom, pos, "-", "Phoenix_Gap", mapq, aln_id)
                    target0 += len(seq)
                elif op == "+":
                    continue
                elif op == "~":
                    match = re.match(r"~[A-Za-z]{2}(\d+)[A-Za-z]{2}", token)
                    if match:
                        length = int(match.group(1))
                        positions = positions_by_chrom.get(chrom, [])
                        left = bisect.bisect_left(positions, target0 + 1)
                        right = bisect.bisect_right(positions, target0 + length)
                        for pos in positions[left:right]:
                            add_call(calls, index_by_site, chrom, pos, "N", "Skipped_Target_Gap", mapq, aln_id)
                        target0 += length
    return calls, stats


def summarize_call(row, calls):
    if not calls:
        return {
            "Phoenix_Base": "NA",
            "Phoenix_Base_Status": "No_Phoenix_Coverage",
            "Phoenix_Alignment_Count": 0,
            "Phoenix_MapQ_Max": "NA",
            "ALT_Polarity_Phoenix": "Unknown_Polarity",
            "ALT_Polarity_Phoenix_Evidence": "No_Phoenix_Alignment_Covers_SNP",
        }
    bases = [call[0].upper() for call in calls]
    base_counts = Counter(bases)
    mapq_max = max(call[2] for call in calls)
    if len(base_counts) > 1:
        return {
            "Phoenix_Base": ";".join(f"{base}:{count}" for base, count in sorted(base_counts.items())),
            "Phoenix_Base_Status": "Ambiguous_Multiple_Bases",
            "Phoenix_Alignment_Count": len(calls),
            "Phoenix_MapQ_Max": mapq_max,
            "ALT_Polarity_Phoenix": "Unknown_Polarity",
            "ALT_Polarity_Phoenix_Evidence": "Multiple_Phoenix_Bases_At_Site",
        }
    base = bases[0]
    status = "Single_Alignment" if len(calls) == 1 else "Multiple_Alignments_Same_Base"
    ref = row["Ref"].upper()
    alt = row["Alt"].upper()
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
        "Phoenix_Alignment_Count": len(calls),
        "Phoenix_MapQ_Max": mapq_max,
        "ALT_Polarity_Phoenix": polarity,
        "ALT_Polarity_Phoenix_Evidence": evidence,
    }


def write_vcf(path, rows):
    with open(path, "w") as out:
        out.write("##fileformat=VCFv4.2\n")
        out.write("##source=41_polarize_dsnp_hap38_phoenix.py\n")
        out.write("##INFO=<ID=PHOENIX_BASE,Number=1,Type=String,Description=\"Phoenix base projected onto Africa_hap2 coordinate\">\n")
        out.write("##INFO=<ID=POLARITY,Number=1,Type=String,Description=\"ALT polarity inferred from Phoenix base\">\n")
        out.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for row in rows:
            info = f"PHOENIX_BASE={row['Phoenix_Base']};POLARITY={row['ALT_Polarity_Phoenix']};EVIDENCE={row['ALT_Polarity_Phoenix_Evidence']}"
            out.write(f"{row['Chrom']}\t{row['Pos']}\t{row['SNP_ID']}\t{row['Ref']}\t{row['Alt']}\t.\tPASS\t{info}\n")


def main():
    args = parse_args()
    outdir = Path(args.outdir)
    if outdir.exists() and any(outdir.iterdir()):
        raise SystemExit(f"[FATAL] Output directory exists and is non-empty: {outdir}")
    outdir.mkdir(parents=True, exist_ok=True)
    qa_dir = outdir / "qa"
    qa_dir.mkdir(parents=True, exist_ok=True)

    rows, input_fields, positions_by_chrom, index_by_site = load_snps(args.rare_functional)
    calls, paf_stats = parse_paf_calls(args.paf, rows, positions_by_chrom, index_by_site, args.min_mapq)

    extra_fields = [
        "Phoenix_Base",
        "Phoenix_Base_Status",
        "Phoenix_Alignment_Count",
        "Phoenix_MapQ_Max",
        "ALT_Polarity_Phoenix",
        "ALT_Polarity_Phoenix_Evidence",
        "dSNP_v2_Phoenix_Flag",
        "dSNP_v2_Phoenix_Evidence",
    ]
    fieldnames = input_fields + extra_fields
    base_fields = [
        "SNP_ID", "Chrom", "Pos", "Ref", "Alt", "Samples", "Functional_Class",
        "Phoenix_Base", "Phoenix_Base_Status", "Phoenix_Alignment_Count", "Phoenix_MapQ_Max",
        "ALT_Polarity_Phoenix", "ALT_Polarity_Phoenix_Evidence",
    ]

    polarized = []
    base_rows = []
    dsnp_v2 = []
    polarity_counts = Counter()
    status_counts = Counter()
    sample_burden = Counter()
    for idx, row in enumerate(rows):
        call_summary = summarize_call(row, calls.get(idx, []))
        out_row = dict(row)
        out_row.update(call_summary)
        if out_row["ALT_Polarity_Phoenix"] == "ALT_Derived":
            out_row["dSNP_v2_Phoenix_Flag"] = "Yes"
            out_row["dSNP_v2_Phoenix_Evidence"] = "Rare_functional_SNP_with_Phoenix_REF_like_base"
            dsnp_v2.append(out_row)
            for sample in row.get("Samples", "").split(";"):
                if sample:
                    sample_burden[sample] += 1
        else:
            out_row["dSNP_v2_Phoenix_Flag"] = "No"
            out_row["dSNP_v2_Phoenix_Evidence"] = "ALT_not_supported_as_derived_by_Phoenix"
        polarized.append(out_row)
        base_rows.append({field: out_row.get(field, "NA") for field in base_fields})
        polarity_counts[out_row["ALT_Polarity_Phoenix"]] += 1
        status_counts[out_row["Phoenix_Base_Status"]] += 1

    write_tsv(outdir / "rare_functional_snp.phoenix_polarized.tsv", polarized, fieldnames)
    write_tsv(outdir / "phoenix_outgroup_base_table.tsv", base_rows, base_fields)
    write_tsv(outdir / "dsnp_v2_phoenix_alt_derived_candidates.tsv", dsnp_v2, fieldnames)
    write_vcf(outdir / "phoenix_outgroup_base_sites.vcf", polarized)

    summary_rows = []
    for key, count in sorted(polarity_counts.items()):
        summary_rows.append({"Summary_Type": "Polarity", "Group": key, "Count": count, "Notes": "rare_functional_snp input"})
    for key, count in sorted(status_counts.items()):
        summary_rows.append({"Summary_Type": "Phoenix_Base_Status", "Group": key, "Count": count, "Notes": "rare_functional_snp input"})
    for sample, count in sorted(sample_burden.items()):
        summary_rows.append({"Summary_Type": "Sample_Burden_ALT_Derived", "Group": sample, "Count": count, "Notes": "dSNP_v2_Phoenix"})
    for key, count in sorted(paf_stats.items()):
        summary_rows.append({"Summary_Type": "PAF_Parse", "Group": key, "Count": count, "Notes": f"min_mapq={args.min_mapq}"})
    summary_rows.append({"Summary_Type": "Global", "Group": "Input_Rare_Functional_SNP", "Count": len(rows), "Notes": str(args.rare_functional)})
    summary_rows.append({"Summary_Type": "Global", "Group": "dSNP_v2_Phoenix_ALT_Derived", "Count": len(dsnp_v2), "Notes": "Strict Phoenix-polarized v2 set"})
    write_tsv(outdir / "phoenix_polarity_summary.tsv", summary_rows, ["Summary_Type", "Group", "Count", "Notes"])

    integrity = [
        {"Check": "input_nonempty", "Status": "PASS" if rows else "FAIL", "Observed": len(rows), "Expected": ">0", "Notes": "rare_functional_snp rows"},
        {"Check": "paf_has_cs", "Status": "PASS" if paf_stats.get("Missing_CS", 0) == 0 else "FAIL", "Observed": paf_stats.get("Missing_CS", 0), "Expected": 0, "Notes": "all used candidate-overlapping alignments must have cs"},
        {"Check": "polarity_output_nonempty", "Status": "PASS" if polarized else "FAIL", "Observed": len(polarized), "Expected": ">0", "Notes": "polarized table rows"},
        {"Check": "alt_derived_nonempty", "Status": "PASS" if dsnp_v2 else "WARN", "Observed": len(dsnp_v2), "Expected": ">0", "Notes": "strict Phoenix ALT-derived candidates"},
    ]
    write_tsv(qa_dir / "final_integrity_summary.tsv", integrity, ["Check", "Status", "Observed", "Expected", "Notes"])
    if any(row["Status"] == "FAIL" for row in integrity):
        raise SystemExit("[FATAL] Integrity check failed")
    (outdir / "_DONE").write_text("Phoenix SNP polarity completed\n")
    print(f"[INFO] Phoenix polarity complete | input={len(rows)} | ALT_Derived={len(dsnp_v2)}")


if __name__ == "__main__":
    main()
