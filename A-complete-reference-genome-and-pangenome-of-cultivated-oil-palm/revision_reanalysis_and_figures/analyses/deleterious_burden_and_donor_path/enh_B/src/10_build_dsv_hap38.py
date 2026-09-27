#!/usr/bin/env python3
"""Reproduce the previous dSV-v5 logic on the user-selected hap38 mixed SV catalog."""

import argparse
import bisect
import csv
import json
import re
from collections import Counter, defaultdict
from pathlib import Path


CONSIDERED_SVTYPES = {"DEL", "INS", "DUP", "INV"}
PRIMARY_OUTGROUP = "Phoenix"
CG_RE = re.compile(r"(\d+)([MIDNSHP=X])")


class PolarityCall:
    def __init__(self, outgroup_state, polarity, evidence):
        self.outgroup_state = outgroup_state
        self.polarity = polarity
        self.evidence = evidence


class InsertionEvidence:
    def __init__(self, state, max_query_insertion, spanning_alignment_count, evidence):
        self.state = state
        self.max_query_insertion = max_query_insertion
        self.spanning_alignment_count = spanning_alignment_count
        self.evidence = evidence


class PafAlignment:
    def __init__(self, chrom, start, end, mapq, cg, insertion_events):
        self.chrom = chrom
        self.start = start
        self.end = end
        self.mapq = mapq
        self.cg = cg
        self.insertion_events = insertion_events


def parse_args():
    parser = argparse.ArgumentParser(
        description="Apply the previous dSV-v5 polarity/function filters to the hap38 mixed SV catalog."
    )
    parser.add_argument("--sv-catalog", required=True, help="Carrier-aware companion TSV for the declared VCF")
    parser.add_argument("--source-vcf", required=True, help="User-declared site-only VCF, retained for provenance")
    parser.add_argument("--phoenix-bed", required=True)
    parser.add_argument("--conserved-bed", required=True)
    parser.add_argument("--gene-bed", required=True)
    parser.add_argument("--exon-bed", required=True)
    parser.add_argument("--cds-bed", required=True)
    parser.add_argument("--promoter-bed", required=True)
    parser.add_argument("--phoenix-paf", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--sample-count", type=int, required=True)
    parser.add_argument("--expected-catalog-rows", type=int, required=True)
    parser.add_argument("--reference-like-threshold", type=float, default=0.80)
    parser.add_argument("--poor-coverage-threshold", type=float, default=0.20)
    parser.add_argument("--rare-frequency-threshold", type=float, default=0.05)
    parser.add_argument("--ins-flank-bp", type=int, default=100)
    parser.add_argument("--ins-site-window-bp", type=int, default=50)
    parser.add_argument("--ins-alt-like-min-fraction", type=float, default=0.50)
    parser.add_argument("--ins-alt-like-min-bp", type=int, default=20)
    parser.add_argument("--bin-size", type=int, default=100000)
    return parser.parse_args()


def as_int(value, default=0):
    try:
        return int(float(str(value)))
    except (TypeError, ValueError):
        return default


def as_float(value, default=0.0):
    try:
        return float(str(value))
    except (TypeError, ValueError):
        return default


def fmt_float(value, digits=7):
    return f"{float(value):.{digits}f}"


def yes_no(value):
    return "Yes" if value else "No"


def split_list(value):
    text = str(value or "").strip()
    if not text or text.upper() == "NA":
        return []
    return [item for item in re.split(r"[;,|]", text) if item]


def coverage_state(fraction, high, low):
    if fraction >= high:
        return "Reference_Like"
    if fraction <= low:
        return "Poorly_Covered"
    return "Partial_Or_Ambiguous"


def read_tsv(path):
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return list(reader), reader.fieldnames or []


def read_vcf_identity_rows(path):
    rows = []
    with open(path) as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 8:
                raise SystemExit(f"[ERROR] Expected site-only 8-column VCF record, observed {len(fields)} columns")
            info = {}
            for item in fields[7].split(";"):
                if "=" in item:
                    key, value = item.split("=", 1)
                    info[key] = value
            rows.append(
                {
                    "SV_Key": fields[2],
                    "Chrom": fields[0],
                    "Start": as_int(fields[1], -1),
                    "End": as_int(info.get("END"), -1),
                    "SVTYPE": info.get("SVTYPE", "NA"),
                    "Sample_Count": as_int(info.get("NS"), -1),
                    "Frequency": as_float(info.get("FREQ"), -1.0),
                    "Evidence_Layer": info.get("EVIDENCE", "NA"),
                }
            )
    return rows


def write_tsv(path, rows, fieldnames):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "NA") for field in fieldnames})


def write_simple_tsv(path, rows, fieldnames):
    with Path(path).open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def read_bed_intervals(path, name_col=None, chrom_whitelist=None):
    intervals = defaultdict(list)
    with open(path) as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                continue
            chrom = fields[0]
            if chrom_whitelist and chrom not in chrom_whitelist:
                continue
            start = as_int(fields[1], None)
            end = as_int(fields[2], None)
            if start is None or end is None or end <= start:
                continue
            name = fields[name_col] if name_col is not None and len(fields) > name_col else "NA"
            intervals[chrom].append((start, end, name))
    for chrom in intervals:
        intervals[chrom].sort(key=lambda item: (item[0], item[1], item[2]))
    return intervals


def build_interval_bins(intervals, bin_size):
    bins = defaultdict(lambda: defaultdict(list))
    for chrom, chrom_intervals in intervals.items():
        for idx, (start, end, _name) in enumerate(chrom_intervals):
            first_bin = start // bin_size
            last_bin = max(end - 1, start) // bin_size
            for bin_id in range(first_bin, last_bin + 1):
                bins[chrom][bin_id].append(idx)
    return bins


def query_interval_hits(intervals, bins, chrom, start, end, bin_size):
    if chrom not in intervals or end <= start:
        return []
    candidate_indexes = set()
    first_bin = start // bin_size
    last_bin = max(end - 1, start) // bin_size
    for bin_id in range(first_bin, last_bin + 1):
        candidate_indexes.update(bins.get(chrom, {}).get(bin_id, []))
    hits = []
    for idx in candidate_indexes:
        iv_start, iv_end, name = intervals[chrom][idx]
        ov_start = max(start, iv_start)
        ov_end = min(end, iv_end)
        if ov_end > ov_start:
            hits.append((iv_start, iv_end, name, ov_end - ov_start))
    return hits


def covered_bases_from_hits(hits, query_start, query_end):
    clipped = []
    for iv_start, iv_end, _name, _overlap in hits:
        start = max(query_start, iv_start)
        end = min(query_end, iv_end)
        if end > start:
            clipped.append((start, end))
    if not clipped:
        return 0
    clipped.sort()
    merged = []
    for start, end in clipped:
        if not merged or start > merged[-1][1]:
            merged.append([start, end])
        else:
            merged[-1][1] = max(merged[-1][1], end)
    return sum(end - start for start, end in merged)


def parse_cg_insertions(tstart, cg):
    target_pos = tstart
    insertions = []
    for length_text, op in CG_RE.findall(cg):
        length = int(length_text)
        if op in {"M", "=", "X"}:
            target_pos += length
        elif op in {"D", "N"}:
            target_pos += length
        elif op == "I":
            insertions.append((target_pos, length))
    return insertions


def read_paf_alignments(path, min_mapq=20):
    by_chrom = defaultdict(list)
    with open(path) as handle:
        for line in handle:
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 12:
                continue
            mapq = as_int(fields[11], -1)
            if mapq < min_mapq:
                continue
            cg = ""
            for tag in fields[12:]:
                if tag.startswith("cg:Z:"):
                    cg = tag[5:]
                    break
            if not cg:
                continue
            chrom = fields[5]
            start = as_int(fields[7])
            end = as_int(fields[8])
            if end <= start:
                continue
            by_chrom[chrom].append(PafAlignment(chrom, start, end, mapq, cg, parse_cg_insertions(start, cg)))
    return by_chrom


def build_alignment_bins(by_chrom, bin_size):
    indexes = {}
    for chrom, alignments in by_chrom.items():
        bins = defaultdict(list)
        for idx, aln in enumerate(alignments):
            first_bin = aln.start // bin_size
            last_bin = max(aln.end - 1, aln.start) // bin_size
            for bin_id in range(first_bin, last_bin + 1):
                bins[bin_id].append(idx)
        starts = sorted((aln.start, idx) for idx, aln in enumerate(alignments))
        indexes[chrom] = {
            "alignments": alignments,
            "bins": bins,
            "starts": [item[0] for item in starts],
            "start_indexes": [item[1] for item in starts],
        }
    return indexes


def collect_candidate_alignments(index, left, right, bin_size):
    first_bin = max(0, left // bin_size)
    last_bin = max(0, right // bin_size)
    candidate_indexes = set()
    for bin_id in range(first_bin, last_bin + 1):
        candidate_indexes.update(index["bins"].get(bin_id, []))
    if candidate_indexes:
        return [index["alignments"][idx] for idx in candidate_indexes]
    upto = bisect.bisect_right(index["starts"], left)
    hits = []
    for idx in reversed(index["start_indexes"][:upto]):
        aln = index["alignments"][idx]
        if aln.end < right:
            continue
        hits.append(aln)
        if len(hits) >= 20:
            break
    return hits


def insertion_alt_like_min_size(abs_svlen, min_fraction, min_bp):
    if abs_svlen <= 0:
        return min_bp
    return max(min_bp, int(round(abs_svlen * min_fraction)))


def evaluate_insertion_evidence(row, paf_index, args):
    chrom = row.get("Chrom", "")
    if chrom not in paf_index:
        return InsertionEvidence("Unknown", 0, 0, "No_Phoenix_PAF_Alignment_On_Chrom")
    abs_svlen = as_int(row.get("Abs_SVLEN"))
    insert_boundary = as_int(row.get("Pos"))
    left = max(0, insert_boundary - args.ins_flank_bp)
    right = insert_boundary + args.ins_flank_bp
    min_alt_size = insertion_alt_like_min_size(abs_svlen, args.ins_alt_like_min_fraction, args.ins_alt_like_min_bp)

    spanning = 0
    max_query_insertion = 0
    for aln in collect_candidate_alignments(paf_index[chrom], left, right, args.bin_size):
        if aln.start > left or aln.end < right:
            continue
        spanning += 1
        for event_pos, event_len in aln.insertion_events:
            if abs(event_pos - insert_boundary) <= args.ins_site_window_bp:
                max_query_insertion = max(max_query_insertion, event_len)

    if spanning == 0:
        return InsertionEvidence("Unknown", max_query_insertion, 0, "No_Phoenix_Alignment_Spans_Insertion_Flanks")
    if max_query_insertion >= min_alt_size:
        return InsertionEvidence("ALT_Like_Insertion_Present", max_query_insertion, spanning, "Phoenix_Has_Query_Insertion_At_Insert_Site")
    return InsertionEvidence("Reference_Like_No_Insertion", max_query_insertion, spanning, "Phoenix_Spans_Insert_Site_Without_ALT_Size_Query_Insertion")


def classify_alt_polarity_v5(row, phoenix_ref_state, insertion_evidence):
    svtype = row.get("SVTYPE", "")
    if svtype == "TRA":
        return PolarityCall("Not_Considered", "Not_Considered", "TRA_Not_Considered_In_NatureStyle_dSV_v5")
    if svtype not in CONSIDERED_SVTYPES:
        return PolarityCall("Not_Considered", "Not_Considered", "SVTYPE_Not_Considered_In_dSV_v5")
    if row.get("Has_Valid_Interval") != "Yes":
        return PolarityCall("Unknown", "Unknown", "No_Valid_Reference_Interval")

    if svtype == "INS":
        if insertion_evidence is None:
            return PolarityCall("Unknown", "Unknown", "Insertion_Local_Evidence_Not_Available")
        if insertion_evidence.state == "Reference_Like_No_Insertion":
            return PolarityCall("REF_Like", "ALT_Derived", insertion_evidence.evidence)
        if insertion_evidence.state == "ALT_Like_Insertion_Present":
            return PolarityCall("ALT_Like", "ALT_Ancestral", insertion_evidence.evidence)
        return PolarityCall("Unknown", "Unknown", insertion_evidence.evidence)

    if svtype == "DEL":
        if phoenix_ref_state == "Reference_Like":
            return PolarityCall("REF_Like", "ALT_Derived", "Phoenix_Reference_Like_Presence_Across_DEL_Interval")
        if phoenix_ref_state == "Poorly_Covered":
            return PolarityCall("ALT_Like", "ALT_Ancestral", "Phoenix_Lacks_Reference_DEL_Interval")
        return PolarityCall("Unknown", "Unknown", "Partial_Phoenix_DEL_Interval_Coverage")

    if svtype == "DUP":
        if phoenix_ref_state == "Reference_Like":
            return PolarityCall("REF_Like", "ALT_Derived", "Phoenix_Reference_Like_Single_Copy_Proxy_Across_DUP_Interval")
        return PolarityCall("Unknown", "Unknown", "DUP_Copy_State_Not_Resolved_By_Phoenix_Reference_Coverage")

    if svtype == "INV":
        if phoenix_ref_state == "Reference_Like":
            return PolarityCall("REF_Like", "ALT_Derived", "Phoenix_Reference_Like_Orientation_Proxy_Across_INV_Interval")
        return PolarityCall("Unknown", "Unknown", "INV_Breakpoint_Orientation_Not_Resolved_By_Phoenix_Coverage")

    return PolarityCall("Unknown", "Unknown", "Unhandled_SVTYPE")


def evidence_class(row, conserved_overlap):
    has_cds = as_int(row.get("CDS_Overlap_Count")) > 0
    if has_cds and conserved_overlap:
        return "CDS_And_Conserved_Region_Proxy"
    if has_cds:
        return "CDS_Overlap"
    if conserved_overlap:
        return "Conserved_Region_Proxy"
    return "Other"


def is_formal_dsv_v5(row, alt_polarity, conserved_overlap, rare_threshold):
    has_cds = as_int(row.get("CDS_Overlap_Count")) > 0
    has_functional_proxy = has_cds or conserved_overlap
    return (
        row.get("SVTYPE") in CONSIDERED_SVTYPES
        and as_float(row.get("Frequency"), 1.0) <= rare_threshold
        and alt_polarity == "ALT_Derived"
        and has_functional_proxy
    )


def feature_overlap_summary(row, feature_sets):
    chrom = row["Chrom"]
    start = as_int(row["Interval_Start0"])
    end = as_int(row["Interval_End0"])
    result = {}
    for label, (intervals, bins) in feature_sets.items():
        hits = query_interval_hits(intervals, bins, chrom, start, end, 100000)
        genes = sorted({hit[2] for hit in hits if hit[2] and hit[2] != "NA"})
        result[f"{label}_Overlap_Count"] = len(genes)
        result[f"{label}_Gene_IDs"] = ";".join(genes) if genes else "NA"
    return result


def impact_class(gene_ids, exon_ids, cds_ids, promoter_ids):
    if cds_ids:
        return "CDS_Overlap"
    if exon_ids:
        return "Exon_Overlap"
    if gene_ids:
        return "Gene_Body_Overlap"
    if promoter_ids:
        return "Promoter_2kb_Overlap"
    return "No_Gene_Context"


def main():
    args = parse_args()
    outdir = Path(args.outdir)
    if outdir.exists() and any(outdir.iterdir()):
        raise SystemExit(f"[ERROR] Refusing to write into non-empty output directory: {outdir}")
    outdir.mkdir(parents=True, exist_ok=True)
    if args.sample_count <= 0:
        raise SystemExit("[ERROR] --sample-count must be positive")

    catalog_rows, catalog_fields = read_tsv(args.sv_catalog)
    required_fields = {
        "Cluster_ID", "Chrom", "Start", "End", "SVTYPE", "SVTYPE_Group",
        "SVLEN_Median_bp", "SVLEN_Min_bp", "SVLEN_Max_bp", "Call_Count",
        "Sample_Count", "Frequency", "Frequency_Bin", "Samples",
        "Dominant_Caller_Combo", "Interval_Start0", "Interval_End0",
        "Repeat_Overlap_bp", "Repeat_Overlap_Fraction", "Dominant_Repeat_Group",
        "Dominant_Repeat_Superfamily", "Evidence_Layer",
    }
    missing_fields = sorted(required_fields - set(catalog_fields))
    if missing_fields:
        raise SystemExit(f"[ERROR] SV catalog is missing required fields: {','.join(missing_fields)}")
    vcf_identity_rows = read_vcf_identity_rows(args.source_vcf)
    if len(vcf_identity_rows) != len(catalog_rows):
        raise SystemExit(
            f"[ERROR] VCF/TSV row-count mismatch: VCF={len(vcf_identity_rows)} TSV={len(catalog_rows)}"
        )
    vcf_tsv_mismatch_count = 0
    for row, identity in zip(catalog_rows, vcf_identity_rows):
        if not (
            identity["Chrom"] == row["Chrom"]
            and identity["Start"] == as_int(row["Start"], -1)
            and identity["End"] == as_int(row["End"], -1)
            and identity["SVTYPE"] == row["SVTYPE"]
            and identity["Sample_Count"] == as_int(row["Sample_Count"], -1)
            and abs(identity["Frequency"] - as_float(row["Frequency"], -1.0)) <= 5.1e-7
            and identity["Evidence_Layer"] == row["Evidence_Layer"]
        ):
            vcf_tsv_mismatch_count += 1
    chroms = {row["Chrom"] for row in catalog_rows}
    unique_source_clusters = {row["Cluster_ID"] for row in catalog_rows}
    unique_sv_keys = {row["SV_Key"] for row in vcf_identity_rows}
    svtype_counts = Counter(row["SVTYPE"] for row in catalog_rows)

    phoenix_intervals = read_bed_intervals(args.phoenix_bed, chrom_whitelist=chroms)
    conserved_intervals = read_bed_intervals(args.conserved_bed, chrom_whitelist=chroms)
    phoenix_bins = build_interval_bins(phoenix_intervals, args.bin_size)
    conserved_bins = build_interval_bins(conserved_intervals, args.bin_size)
    feature_sets = {
        "Gene": (read_bed_intervals(args.gene_bed, name_col=3, chrom_whitelist=chroms), None),
        "Exon": (read_bed_intervals(args.exon_bed, name_col=3, chrom_whitelist=chroms), None),
        "CDS": (read_bed_intervals(args.cds_bed, name_col=3, chrom_whitelist=chroms), None),
        "Promoter_2kb": (read_bed_intervals(args.promoter_bed, name_col=3, chrom_whitelist=chroms), None),
    }
    feature_sets = {label: (intervals, build_interval_bins(intervals, args.bin_size)) for label, (intervals, _bins) in feature_sets.items()}
    paf_index = build_alignment_bins(read_paf_alignments(args.phoenix_paf), args.bin_size)

    all_rows = []
    candidate_rows = []
    summary = Counter()
    by_svtype = defaultdict(Counter)
    by_evidence = defaultdict(Counter)
    polarity_audit = Counter()
    frequency_mismatch_count = 0
    carrier_mismatch_count = 0
    sample_count_range_bad = 0

    for row, identity in zip(catalog_rows, vcf_identity_rows):
        sv_key = identity["SV_Key"]
        source_cluster_id = row["Cluster_ID"]
        svtype = row["SVTYPE"]
        interval_start = as_int(row.get("Interval_Start0"))
        interval_end = as_int(row.get("Interval_End0"))
        interval_len = max(0, interval_end - interval_start)
        has_valid_interval = interval_len > 0
        abs_svlen = abs(as_int(row.get("SVLEN_Median_bp")))
        pos = as_int(row.get("Start"))

        base = {
            "SV_Key": sv_key,
            "SV_ID": sv_key,
            "Source_Cluster_ID": source_cluster_id,
            "Original_Chrom": row.get("Chrom", "NA"),
            "Chrom": row.get("Chrom", "NA"),
            "Pos": pos,
            "Start": row.get("Start", "NA"),
            "End": row.get("End", "NA"),
            "Original_Chr2": "NA",
            "Chr2": "NA",
            "SVTYPE": svtype,
            "SVTYPE_Group": row.get("SVTYPE_Group", svtype),
            "SVLEN": row.get("SVLEN_Median_bp", "NA"),
            "Abs_SVLEN": abs_svlen,
            "SVLEN_Min_bp": row.get("SVLEN_Min_bp", "NA"),
            "SVLEN_Max_bp": row.get("SVLEN_Max_bp", "NA"),
            "SUPP": row.get("Sample_Count", "NA"),
            "Call_Count": row.get("Call_Count", "NA"),
            "Sample_Count": row.get("Sample_Count", "NA"),
            "Frequency": row.get("Frequency", "NA"),
            "Frequency_Bin": row.get("Frequency_Bin", "NA"),
            "Samples": row.get("Samples", "NA"),
            "Dominant_Caller_Combo": row.get("Dominant_Caller_Combo", "NA"),
            "Evidence_Layer": row.get("Evidence_Layer", "NA"),
            "Is_Private": yes_no(as_int(row.get("Sample_Count")) == 1),
            "Is_Rare_05": yes_no(as_float(row.get("Frequency"), 1.0) <= args.rare_frequency_threshold),
            "Is_LowFreq_10": yes_no(as_float(row.get("Frequency"), 1.0) <= 0.10),
            "Is_Intra_Chrom": "Yes",
            "Has_Valid_Interval": yes_no(has_valid_interval),
            "Interval_Start0": interval_start,
            "Interval_End0": interval_end,
            "Interval_Length_bp": interval_len,
            "Repeat_Overlap_bp": row.get("Repeat_Overlap_bp", "NA"),
            "Repeat_Overlap_Fraction": row.get("Repeat_Overlap_Fraction", "NA"),
            "Dominant_Repeat_Group": row.get("Dominant_Repeat_Group", "NA"),
            "Dominant_Repeat_Superfamily": row.get("Dominant_Repeat_Superfamily", "NA"),
            "Note": "Hap38_User_Selected_Mixed_HighConfidence_Catalog",
        }

        sample_count = as_int(row.get("Sample_Count"), -1)
        carriers = split_list(row.get("Samples"))
        if sample_count < 1 or sample_count > args.sample_count:
            sample_count_range_bad += 1
        if len(carriers) != sample_count or len(set(carriers)) != len(carriers):
            carrier_mismatch_count += 1
        if sample_count >= 0 and abs(as_float(row.get("Frequency"), -1.0) - sample_count / args.sample_count) > 1e-9:
            frequency_mismatch_count += 1

        phoenix_hits = query_interval_hits(phoenix_intervals, phoenix_bins, base["Chrom"], interval_start, interval_end, args.bin_size)
        phoenix_bases = covered_bases_from_hits(phoenix_hits, interval_start, interval_end)
        phoenix_frac = phoenix_bases / interval_len if interval_len else 0.0
        phoenix_state = coverage_state(phoenix_frac, args.reference_like_threshold, args.poor_coverage_threshold)

        conserved_hits = query_interval_hits(conserved_intervals, conserved_bins, base["Chrom"], interval_start, interval_end, args.bin_size)
        conserved_bases = covered_bases_from_hits(conserved_hits, interval_start, interval_end)
        conserved_frac = conserved_bases / interval_len if interval_len else 0.0
        conserved_overlap = conserved_bases > 0

        feature_summary = feature_overlap_summary(base, feature_sets)
        gene_ids = split_list(feature_summary["Gene_Gene_IDs"])
        exon_ids = split_list(feature_summary["Exon_Gene_IDs"])
        cds_ids = split_list(feature_summary["CDS_Gene_IDs"])
        promoter_ids = split_list(feature_summary["Promoter_2kb_Gene_IDs"])

        base.update(
            {
                "Phoenix_Overlap_Count": len(phoenix_hits),
                "Phoenix_Covered_Bases": phoenix_bases,
                "SV_Bed_Length": interval_len,
                "Phoenix_Ref_Coverage_Fraction": fmt_float(phoenix_frac),
                "Phoenix_RefCov_v5": fmt_float(phoenix_frac),
                "Phoenix_Reference_State_v5": phoenix_state,
                "Conserved_Region_Covered_Bases": conserved_bases,
                "Conserved_Region_Overlap_Fraction": fmt_float(conserved_frac),
                "Conserved_Region_Overlap": yes_no(conserved_overlap),
                "Gene_Overlap_Count": feature_summary["Gene_Overlap_Count"],
                "Exon_Overlap_Count": feature_summary["Exon_Overlap_Count"],
                "CDS_Overlap_Count": feature_summary["CDS_Overlap_Count"],
                "Promoter_2kb_Overlap_Count": feature_summary["Promoter_2kb_Overlap_Count"],
                "Gene_IDs": feature_summary["Gene_Gene_IDs"],
                "Exon_Gene_IDs": feature_summary["Exon_Gene_IDs"],
                "CDS_Gene_IDs": feature_summary["CDS_Gene_IDs"],
                "Promoter_2kb_Gene_IDs": feature_summary["Promoter_2kb_Gene_IDs"],
                "Gene_Impact_Class": impact_class(gene_ids, exon_ids, cds_ids, promoter_ids),
                "Primary_Outgroup": PRIMARY_OUTGROUP,
            }
        )

        ins_evidence = evaluate_insertion_evidence(base, paf_index, args) if svtype == "INS" else None
        call = classify_alt_polarity_v5(base, phoenix_state, ins_evidence)
        flag = is_formal_dsv_v5(base, call.polarity, conserved_overlap, args.rare_frequency_threshold)
        ev_class = evidence_class(base, conserved_overlap) if flag else "Other"
        base.update(
            {
                "Phoenix_INS_Local_State": ins_evidence.state if ins_evidence else "NA",
                "Phoenix_INS_Spanning_Alignment_Count": ins_evidence.spanning_alignment_count if ins_evidence else "NA",
                "Phoenix_INS_Max_Query_Insertion": ins_evidence.max_query_insertion if ins_evidence else "NA",
                "Outgroup_Allele_State": call.outgroup_state,
                "ALT_Polarity_v5": call.polarity,
                "ALT_Polarity_v5_Evidence": call.evidence,
                "dSV_v5_Flag": yes_no(flag),
                "dSV_v5_Evidence_Class": ev_class,
            }
        )

        all_rows.append(base)
        if flag:
            candidate_rows.append(base)

        evidence_layer = base["Evidence_Layer"]
        summary["Total_Input_Clusters"] += 1
        summary[f"SVTYPE_{svtype}"] += 1
        summary[f"Evidence_Layer_{evidence_layer}"] += 1
        if svtype in CONSIDERED_SVTYPES:
            summary["Considered_SV_DEL_INS_DUP_INV"] += 1
        else:
            summary["Not_Considered_SV"] += 1
        if as_float(base["Frequency"], 1.0) <= args.rare_frequency_threshold:
            summary["Rare_05_SV"] += 1
        # Not_Considered_SV is already counted by SVTYPE above. Counting the
        # TRA polarity again would double this global summary category.
        if call.polarity != "Not_Considered":
            summary[f"{call.polarity}_SV"] += 1
        if flag:
            summary["dSV_v5_Total"] += 1
            summary[f"dSV_v5_{ev_class}"] += 1
        by_svtype[svtype]["Total_Count"] += 1
        by_svtype[svtype]["Rare_05_Count"] += int(as_float(base["Frequency"], 1.0) <= args.rare_frequency_threshold)
        by_svtype[svtype][f"{call.polarity}_Count"] += 1
        by_svtype[svtype]["dSV_v5_Count"] += int(flag)
        by_svtype[svtype][f"dSV_v5_{ev_class}_Count"] += int(flag)
        by_evidence[evidence_layer]["Total_Count"] += 1
        by_evidence[evidence_layer]["Rare_05_Count"] += int(as_float(base["Frequency"], 1.0) <= args.rare_frequency_threshold)
        by_evidence[evidence_layer]["ALT_Derived_Count"] += int(call.polarity == "ALT_Derived")
        by_evidence[evidence_layer]["dSV_v5_Count"] += int(flag)
        polarity_audit[(svtype, call.polarity, call.evidence)] += 1

    fields = [
        "SV_Key", "SV_ID", "Source_Cluster_ID", "Original_Chrom", "Chrom", "Pos", "Start", "End", "Original_Chr2", "Chr2",
        "SVTYPE", "SVTYPE_Group", "SVLEN", "Abs_SVLEN", "SVLEN_Min_bp", "SVLEN_Max_bp", "SUPP", "Call_Count",
        "Sample_Count", "Frequency", "Frequency_Bin", "Samples", "Dominant_Caller_Combo", "Evidence_Layer", "Is_Private",
        "Is_Rare_05", "Is_LowFreq_10", "Is_Intra_Chrom", "Has_Valid_Interval", "Interval_Start0",
        "Interval_End0", "Interval_Length_bp", "Repeat_Overlap_bp", "Repeat_Overlap_Fraction",
        "Dominant_Repeat_Group", "Dominant_Repeat_Superfamily", "Note", "Phoenix_Overlap_Count",
        "Phoenix_Covered_Bases", "SV_Bed_Length", "Phoenix_Ref_Coverage_Fraction", "Phoenix_RefCov_v5",
        "Phoenix_Reference_State_v5", "Phoenix_INS_Local_State", "Phoenix_INS_Spanning_Alignment_Count",
        "Phoenix_INS_Max_Query_Insertion", "Outgroup_Allele_State", "ALT_Polarity_v5",
        "ALT_Polarity_v5_Evidence", "Conserved_Region_Overlap", "Conserved_Region_Covered_Bases",
        "Conserved_Region_Overlap_Fraction", "Gene_Overlap_Count", "Exon_Overlap_Count", "CDS_Overlap_Count",
        "Promoter_2kb_Overlap_Count", "Gene_IDs", "Exon_Gene_IDs", "CDS_Gene_IDs",
        "Promoter_2kb_Gene_IDs", "Gene_Impact_Class", "dSV_v5_Flag", "dSV_v5_Evidence_Class",
    ]

    write_tsv(outdir / "sv_catalog.dsv_hap38.tsv", all_rows, fields)
    write_tsv(outdir / "dsv_hap38_candidates.tsv", candidate_rows, fields)

    summary_rows = [{"Metric": key, "Value": value} for key, value in sorted(summary.items())]
    summary_rows.extend(
        [
            {"Metric": "Unique_VCF_SV_Key", "Value": len(unique_sv_keys)},
            {"Metric": "Unique_Source_Cluster_ID", "Value": len(unique_source_clusters)},
            {"Metric": "Source_Cluster_ID_Collisions", "Value": len(catalog_rows) - len(unique_source_clusters)},
            {"Metric": "Input_Catalog_Rows", "Value": len(catalog_rows)},
            {"Metric": "Rare_Frequency_Threshold", "Value": args.rare_frequency_threshold},
            {"Metric": "Sample_Count_Denominator", "Value": args.sample_count},
        ]
    )
    write_simple_tsv(outdir / "dsv_hap38_summary.tsv", summary_rows, ["Metric", "Value"])

    svtype_fields = [
        "SVTYPE", "Total_Count", "Rare_05_Count", "ALT_Derived_Count", "ALT_Ancestral_Count",
        "Unknown_Count", "Not_Considered_Count", "dSV_v5_Count", "dSV_v5_CDS_Overlap_Count",
        "dSV_v5_Conserved_Region_Proxy_Count", "dSV_v5_CDS_And_Conserved_Region_Proxy_Count",
    ]
    svtype_rows = []
    for svtype in sorted(by_svtype, key=lambda x: ["DEL", "INS", "DUP", "INV", "TRA"].index(x) if x in {"DEL", "INS", "DUP", "INV", "TRA"} else 99):
        row = {"SVTYPE": svtype}
        row.update({field: by_svtype[svtype].get(field, 0) for field in svtype_fields if field != "SVTYPE"})
        svtype_rows.append(row)
    write_simple_tsv(outdir / "dsv_hap38_summary_by_svtype.tsv", svtype_rows, svtype_fields)

    evidence_fields = ["Evidence_Layer", "Total_Count", "Rare_05_Count", "ALT_Derived_Count", "dSV_v5_Count"]
    evidence_rows = []
    for evidence_layer in sorted(by_evidence):
        evidence_rows.append(
            {"Evidence_Layer": evidence_layer, **{field: by_evidence[evidence_layer].get(field, 0) for field in evidence_fields[1:]}}
        )
    write_simple_tsv(outdir / "dsv_hap38_summary_by_evidence.tsv", evidence_rows, evidence_fields)

    audit_rows = [
        {"SVTYPE": svtype, "ALT_Polarity_v5": polarity, "Evidence": evidence, "SV_Count": count}
        for (svtype, polarity, evidence), count in sorted(polarity_audit.items())
    ]
    write_simple_tsv(outdir / "dsv_hap38_polarity_audit.tsv", audit_rows, ["SVTYPE", "ALT_Polarity_v5", "Evidence", "SV_Count"])

    svtype_total = sum(counts.get("Total_Count", 0) for counts in by_svtype.values())
    considered_total = summary.get("Considered_SV_DEL_INS_DUP_INV", 0) + summary.get("Not_Considered_SV", 0)
    polarity_total = sum(summary.get(f"{polarity}_SV", 0) for polarity in ("ALT_Derived", "ALT_Ancestral", "Unknown", "Not_Considered"))
    evidence_total = sum(counts.get("Total_Count", 0) for counts in by_evidence.values())
    polarity_audit_total = sum(polarity_audit.values())
    dsv_summary_total = summary.get("dSV_v5_Total", 0)
    dsv_svtype_total = sum(counts.get("dSV_v5_Count", 0) for counts in by_svtype.values())

    integrity_rows = [
        {"Check_Name": "Input_Catalog_Rows", "Observed_Value": len(catalog_rows), "Expected_Value": args.expected_catalog_rows, "Status": "Pass" if len(catalog_rows) == args.expected_catalog_rows else "Fail", "Note": "Rows excluding header"},
        {"Check_Name": "Output_Catalog_Total_Conservation", "Observed_Value": len(all_rows), "Expected_Value": len(catalog_rows), "Status": "Pass" if len(all_rows) == len(catalog_rows) else "Fail", "Note": "Output catalog rows equal input catalog rows"},
        {"Check_Name": "SVTYPE_Total_Conservation", "Observed_Value": svtype_total, "Expected_Value": len(catalog_rows), "Status": "Pass" if svtype_total == len(catalog_rows) else "Fail", "Note": "Sum of per-SVTYPE Total_Count equals input rows"},
        {"Check_Name": "Considered_Total_Conservation", "Observed_Value": considered_total, "Expected_Value": len(catalog_rows), "Status": "Pass" if considered_total == len(catalog_rows) else "Fail", "Note": "Considered DEL/INS/DUP/INV plus not-considered TRA equals input rows"},
        {"Check_Name": "Polarity_Total_Conservation", "Observed_Value": polarity_total, "Expected_Value": len(catalog_rows), "Status": "Pass" if polarity_total == len(catalog_rows) else "Fail", "Note": "ALT_Derived plus ALT_Ancestral plus Unknown plus Not_Considered equals input rows"},
        {"Check_Name": "Evidence_Layer_Total_Conservation", "Observed_Value": evidence_total, "Expected_Value": len(catalog_rows), "Status": "Pass" if evidence_total == len(catalog_rows) else "Fail", "Note": "Sum of evidence-layer totals equals input rows"},
        {"Check_Name": "Polarity_Audit_Total_Conservation", "Observed_Value": polarity_audit_total, "Expected_Value": len(catalog_rows), "Status": "Pass" if polarity_audit_total == len(catalog_rows) else "Fail", "Note": "Sum of detailed polarity-audit categories equals input rows"},
        {"Check_Name": "dSV_Summary_Candidate_Conservation", "Observed_Value": dsv_summary_total, "Expected_Value": len(candidate_rows), "Status": "Pass" if dsv_summary_total == len(candidate_rows) else "Fail", "Note": "Global dSV summary count equals candidate-table rows"},
        {"Check_Name": "dSV_SVTYPE_Candidate_Conservation", "Observed_Value": dsv_svtype_total, "Expected_Value": len(candidate_rows), "Status": "Pass" if dsv_svtype_total == len(candidate_rows) else "Fail", "Note": "Sum of per-SVTYPE dSV counts equals candidate-table rows"},
        {"Check_Name": "Unique_VCF_SV_Key", "Observed_Value": len(unique_sv_keys), "Expected_Value": len(catalog_rows), "Status": "Pass" if len(unique_sv_keys) == len(catalog_rows) else "Fail", "Note": "Use OP39SV IDs as output keys"},
        {"Check_Name": "Source_Cluster_ID_Collisions", "Observed_Value": len(catalog_rows) - len(unique_source_clusters), "Expected_Value": "Recorded", "Status": "Pass", "Note": "Expected namespace reuse across mixed evidence layers; retained only for provenance"},
        {"Check_Name": "VCF_TSV_One_To_One", "Observed_Value": vcf_tsv_mismatch_count, "Expected_Value": 0, "Status": "Pass" if vcf_tsv_mismatch_count == 0 else "Fail", "Note": "Row-order coordinate/type/count/evidence agreement"},
        {"Check_Name": "Sample_Count_Denominator", "Observed_Value": args.sample_count, "Expected_Value": args.sample_count, "Status": "Pass", "Note": "Frequency denominator"},
        {"Check_Name": "Frequency_Consistency", "Observed_Value": frequency_mismatch_count, "Expected_Value": 0, "Status": "Pass" if frequency_mismatch_count == 0 else "Fail", "Note": "Frequency equals Sample_Count/sample_count"},
        {"Check_Name": "Carrier_Count_Consistency", "Observed_Value": carrier_mismatch_count, "Expected_Value": 0, "Status": "Pass" if carrier_mismatch_count == 0 else "Fail", "Note": "Unique Samples count equals Sample_Count"},
        {"Check_Name": "Sample_Count_Range", "Observed_Value": sample_count_range_bad, "Expected_Value": 0, "Status": "Pass" if sample_count_range_bad == 0 else "Fail", "Note": "Sample_Count is within 1..denominator"},
        {"Check_Name": "Allowed_SVTYPES", "Observed_Value": sum(count for svtype, count in svtype_counts.items() if svtype not in {"DEL", "INS", "DUP", "INV", "TRA"}), "Expected_Value": 0, "Status": "Pass" if set(svtype_counts) <= {"DEL", "INS", "DUP", "INV", "TRA"} else "Fail", "Note": "Expected mixed-catalog SVTYPE set"},
    ]
    for svtype in ("DEL", "INS", "INV", "DUP", "TRA"):
        integrity_rows.append({"Check_Name": f"SVTYPE_{svtype}", "Observed_Value": svtype_counts.get(svtype, 0), "Expected_Value": "Recorded", "Status": "Pass", "Note": "Observed input count; source checksum is the count contract"})
    bad_candidates = [
        row["SV_Key"] for row in candidate_rows
        if row["SVTYPE"] not in CONSIDERED_SVTYPES
        or as_int(row["Sample_Count"]) != 1
        or row["ALT_Polarity_v5"] != "ALT_Derived"
        or (as_int(row["CDS_Overlap_Count"]) <= 0 and row["Conserved_Region_Overlap"] != "Yes")
    ]
    integrity_rows.append({"Check_Name": "dSV_v5_Formal_Filter", "Observed_Value": len(bad_candidates), "Expected_Value": 0, "Status": "Pass" if not bad_candidates else "Fail", "Note": "Bad candidate rows"})
    write_simple_tsv(outdir / "dsv_hap38_integrity_summary.tsv", integrity_rows, ["Check_Name", "Observed_Value", "Expected_Value", "Status", "Note"])

    metadata = {
        "Analysis": "Hap38 mixed-catalog dSV calls using the previous dSV-v5 polarity and functional filters",
        "Source_VCF": str(args.source_vcf),
        "Carrier_Aware_SV_Catalog": str(args.sv_catalog),
        "Definition": "rare singleton DEL/INS/DUP/INV + Phoenix-supported ALT-derived + CDS or conserved-region proxy overlap",
        "Rare_Frequency_Threshold": args.rare_frequency_threshold,
        "Sample_Count": args.sample_count,
        "Expected_Catalog_Rows": args.expected_catalog_rows,
        "Primary_Outgroup": PRIMARY_OUTGROUP,
        "Conserved_Region_Label": "multi-outgroup conserved/alignable-region proxy; not GERP, phyloP, or phastCons",
        "Identifier_Note": "SV_Key/SV_ID use unique OP39SV IDs from the declared VCF; Source_Cluster_ID preserves the layer-local TSV identifier.",
        "Evidence_Layer_Note": "The selected mixed catalog contains Tier1 read-plus-assembly DEL/INS and SyRI large-rearrangement INV/DUP/TRA layers.",
        "Compatibility_Note": "SUPP is retained as a compatibility alias for carrier Sample_Count, not caller support.",
    }
    (outdir / "notes").mkdir(exist_ok=True)
    (outdir / "notes" / "analysis_metadata.json").write_text(json.dumps(metadata, indent=2, sort_keys=True) + "\n")
    (outdir / "notes" / "analysis_notes.md").write_text(
        "# Hap38 dSV notes\n\n"
        "This analysis uses the user-selected 178,314-record mixed high-confidence SV catalog and "
        "reproduces the previous dSV-v5 rarity, Phoenix polarity, and CDS/conserved-proxy filters. "
        "Formal dSVs consider DEL, INS, DUP, and INV; TRA is retained only in input and polarity audits. "
        "Evidence_Layer must remain attached to every record because the mixed catalog does not apply "
        "one identical caller-support rule to all SV types. The conserved-region track is a multi-outgroup "
        "conserved/alignable-region proxy, not GERP, phyloP, or phastCons.\n"
    )
    failures = [row for row in integrity_rows if row["Status"] == "Fail"]
    if failures:
        raise SystemExit(f"[ERROR] Integrity checks failed: {failures[:3]}")
    print(f"[OK] Wrote hap38 dSV catalog to {outdir}")
    print(f"[OK] Total input clusters: {len(all_rows)}; dSV candidates: {len(candidate_rows)}")


if __name__ == "__main__":
    main()
