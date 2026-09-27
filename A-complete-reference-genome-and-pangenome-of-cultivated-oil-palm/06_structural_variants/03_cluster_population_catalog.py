#!/usr/bin/env python3

import argparse
import csv
import json
import os
import random
import statistics
from collections import Counter, defaultdict


SVTYPE_ORDER = ["DEL", "INS", "INV", "DUP", "TRA", "HDR", "INVDP", "TDM"]
SVTYPE_RANK = {name: idx for idx, name in enumerate(SVTYPE_ORDER)}
PLOT_GROUPS = ["All", "DEL", "INS", "INV", "DUP", "Other"]
PLOT_TYPE_GROUPS = ["DEL", "INS", "INV", "DUP", "Other"]
SPAN_TYPES = {"DEL", "DUP", "INV", "HDR"}


class SVCall:
    __slots__ = (
        "idx",
        "sample",
        "chrom",
        "start",
        "end",
        "svtype",
        "svlen",
        "svlen_raw",
        "supp",
        "caller_combo",
        "size_bin",
    )

    def __init__(
        self,
        idx,
        sample,
        chrom,
        start,
        end,
        svtype,
        svlen,
        svlen_raw,
        supp,
        caller_combo,
        size_bin,
    ):
        self.idx = idx
        self.sample = sample
        self.chrom = chrom
        self.start = start
        self.end = end
        self.svtype = svtype
        self.svlen = svlen
        self.svlen_raw = svlen_raw
        self.supp = supp
        self.caller_combo = caller_combo
        self.size_bin = size_bin


class Cluster:
    __slots__ = (
        "internal_id",
        "chrom",
        "svtype",
        "starts",
        "ends",
        "svlens",
        "samples",
        "call_count",
        "supp_values",
        "caller_combos",
        "size_bins",
        "rep_start",
        "rep_end",
        "rep_len",
        "min_start",
        "max_start",
        "min_end",
        "max_end",
        "interval_start0",
        "interval_end0",
        "interval_len",
        "repeat_overlap_bp",
        "repeat_overlap_fraction",
        "top_repeat_group",
        "top_repeat_superfamily",
        "top_repeat_bp",
    )

    def __init__(self, internal_id, call):
        self.internal_id = internal_id
        self.chrom = call.chrom
        self.svtype = call.svtype
        self.starts = []
        self.ends = []
        self.svlens = []
        self.samples = set()
        self.call_count = 0
        self.supp_values = []
        self.caller_combos = Counter()
        self.size_bins = Counter()
        self.rep_start = call.start
        self.rep_end = normalized_end(call)
        self.rep_len = max(1, call.svlen)
        self.min_start = call.start
        self.max_start = call.start
        self.min_end = normalized_end(call)
        self.max_end = normalized_end(call)
        self.interval_start0 = 0
        self.interval_end0 = 1
        self.interval_len = 1
        self.repeat_overlap_bp = 0
        self.repeat_overlap_fraction = 0.0
        self.top_repeat_group = "No repeat"
        self.top_repeat_superfamily = "NA"
        self.top_repeat_bp = 0
        self.add_call(call)

    def add_call(self, call):
        end = normalized_end(call)
        self.starts.append(call.start)
        self.ends.append(end)
        self.svlens.append(max(1, call.svlen))
        self.samples.add(call.sample)
        self.call_count += 1
        self.supp_values.append(call.supp)
        if call.caller_combo:
            self.caller_combos[call.caller_combo] += 1
        if call.size_bin:
            self.size_bins[call.size_bin] += 1
        self.rep_start = int(round(sum(self.starts) / len(self.starts)))
        self.rep_end = int(round(sum(self.ends) / len(self.ends)))
        self.rep_len = median_int(self.svlens)
        self.min_start = min(self.min_start, call.start)
        self.max_start = max(self.max_start, call.start)
        self.min_end = min(self.min_end, end)
        self.max_end = max(self.max_end, end)


class RepeatRecord:
    __slots__ = ("start", "end", "group", "superfamily")

    def __init__(self, start, end, group, superfamily):
        self.start = start
        self.end = end
        self.group = group
        self.superfamily = superfamily


def parse_args():
    parser = argparse.ArgumentParser(
        description="Build a cross-sample SUPP>=2 SV catalog, frequency/saturation tables, and repeat-overlap summaries."
    )
    parser.add_argument(
        "--records",
        default="results/highconf_sv.records.tsv",
        help="SUPP>=2 per-sample SV record TSV.",
    )
    parser.add_argument(
        "--repeat-bed",
        default="FL_Hap2.fa.mod.EDTA.TEanno.bed",
        help="EDTA TE annotation BED. Use 'none' to skip repeat overlap.",
    )
    parser.add_argument(
        "--fai",
        default="FL_Hap2.fa.fai",
        help="Reference FASTA index for chromosome ordering.",
    )
    parser.add_argument(
        "--out-dir",
        default="results/population_repeat",
        help="Output directory.",
    )
    parser.add_argument("--merge-distance", type=int, default=1000)
    parser.add_argument("--ins-distance", type=int, default=500)
    parser.add_argument("--min-reciprocal-overlap", type=float, default=0.50)
    parser.add_argument("--min-size-similarity", type=float, default=0.50)
    parser.add_argument(
        "--ins-repeat-window",
        type=int,
        default=50,
        help="Reference breakpoint window on each side for INS repeat context.",
    )
    parser.add_argument("--saturation-iterations", type=int, default=200)
    parser.add_argument("--seed", type=int, default=20260625)
    parser.add_argument("--force", action="store_true", help="Overwrite existing outputs.")
    return parser.parse_args()


def safe_int(value, default=0):
    try:
        if value is None or value == "" or value == "NA":
            return default
        return int(float(value))
    except ValueError:
        return default


def median_int(values):
    if not values:
        return 0
    return int(round(statistics.median(values)))


def normalized_end(call):
    if call.svtype == "INS":
        return call.start
    return max(call.start + 1, call.end)


def call_span(call):
    if call.svtype == "INS":
        return call.start, call.start + 1
    end = normalized_end(call)
    return call.start, max(call.start + 1, end)


def cluster_span(cluster):
    if cluster.svtype == "INS":
        return cluster.rep_start, cluster.rep_start + 1
    return cluster.rep_start, max(cluster.rep_start + 1, cluster.rep_end)


def size_similarity(a, b):
    if a <= 0 or b <= 0:
        return 1.0
    return min(a, b) / max(a, b)


def reciprocal_overlap(s1, e1, s2, e2):
    len1 = max(1, e1 - s1)
    len2 = max(1, e2 - s2)
    overlap = max(0, min(e1, e2) - max(s1, s2))
    if overlap <= 0:
        return 0.0
    return min(overlap / len1, overlap / len2)


def match_cluster(call, cluster, args):
    if call.svtype != cluster.svtype or call.chrom != cluster.chrom:
        return False, float("inf")

    call_len = max(1, call.svlen)
    size_sim = size_similarity(call_len, cluster.rep_len)
    if size_sim < args.min_size_similarity:
        return False, float("inf")

    if call.svtype == "INS":
        delta = abs(call.start - cluster.rep_start)
        if delta <= args.ins_distance:
            len_delta = abs(call_len - cluster.rep_len) / max(call_len, cluster.rep_len)
            return True, delta + len_delta
        return False, float("inf")

    if call.svtype in SPAN_TYPES:
        s1, e1 = call_span(call)
        s2, e2 = cluster_span(cluster)
        endpoint_delta = abs(s1 - s2) + abs(e1 - e2)
        endpoints_close = abs(s1 - s2) <= args.merge_distance and abs(e1 - e2) <= args.merge_distance
        ro = reciprocal_overlap(s1, e1, s2, e2)
        if endpoints_close or ro >= args.min_reciprocal_overlap:
            return True, endpoint_delta - (ro * 1000.0)
        return False, float("inf")

    delta = abs(call.start - cluster.rep_start)
    if delta <= args.merge_distance:
        return True, delta
    return False, float("inf")


def read_chrom_order(fai_path):
    order = {}
    lengths = {}
    if fai_path and os.path.exists(fai_path):
        with open(fai_path, "r", encoding="utf-8") as handle:
            for idx, line in enumerate(handle):
                fields = line.rstrip("\n").split("\t")
                if len(fields) >= 2:
                    order[fields[0]] = idx
                    lengths[fields[0]] = safe_int(fields[1])
    return order, lengths


def chrom_rank(chrom, chrom_order):
    return chrom_order.get(chrom, len(chrom_order) + 1000)


def read_records(path):
    calls = []
    with open(path, "r", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"Sample", "Chrom", "Start", "End", "SVTYPE", "SVLEN_bp", "SUPP", "Caller_Combo", "Size_Bin"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise SystemExit(f"Missing columns in records TSV: {', '.join(sorted(missing))}")
        for idx, row in enumerate(reader):
            svtype = row["SVTYPE"]
            if svtype not in SVTYPE_RANK:
                continue
            start = safe_int(row["Start"])
            end = safe_int(row["End"], start)
            svlen_raw = safe_int(row.get("SVLEN_raw", row["SVLEN_bp"]))
            svlen = abs(safe_int(row["SVLEN_bp"]))
            if svlen <= 0:
                svlen = max(1, abs(end - start))
            calls.append(
                SVCall(
                    idx=idx,
                    sample=row["Sample"],
                    chrom=row["Chrom"],
                    start=start,
                    end=end,
                    svtype=svtype,
                    svlen=svlen,
                    svlen_raw=svlen_raw,
                    supp=safe_int(row["SUPP"]),
                    caller_combo=row["Caller_Combo"],
                    size_bin=row["Size_Bin"],
                )
            )
    return calls


def cluster_calls(calls, chrom_order, args):
    sorted_calls = sorted(
        calls,
        key=lambda c: (chrom_rank(c.chrom, chrom_order), c.chrom, SVTYPE_RANK[c.svtype], c.start, c.end, c.sample),
    )
    assignments = {}
    clusters = []
    active = []
    current_key = None
    next_internal_id = 1

    for call in sorted_calls:
        key = (call.chrom, call.svtype)
        if key != current_key:
            clusters.extend(active)
            active = []
            current_key = key

        if call.svtype == "INS":
            keep_active = []
            for cluster in active:
                if cluster.rep_start + args.ins_distance >= call.start:
                    keep_active.append(cluster)
                else:
                    clusters.append(cluster)
            active = keep_active
        else:
            keep_active = []
            for cluster in active:
                if cluster.max_end + args.merge_distance >= call.start:
                    keep_active.append(cluster)
                else:
                    clusters.append(cluster)
            active = keep_active

        best_cluster = None
        best_score = float("inf")
        for cluster in active:
            ok, score = match_cluster(call, cluster, args)
            if ok and score < best_score:
                best_cluster = cluster
                best_score = score

        if best_cluster is None:
            best_cluster = Cluster(next_internal_id, call)
            next_internal_id += 1
            active.append(best_cluster)
        else:
            best_cluster.add_call(call)

        assignments[call.idx] = best_cluster.internal_id

    clusters.extend(active)
    clusters.sort(key=lambda c: (chrom_rank(c.chrom, chrom_order), c.chrom, c.rep_start, c.rep_end, SVTYPE_RANK[c.svtype]))
    return clusters, assignments


def size_bin_from_len(length_bp):
    if length_bp < 1000:
        return "50bp-1kb"
    if length_bp < 10000:
        return "1-10kb"
    if length_bp < 100000:
        return "10-100kb"
    if length_bp < 1000000:
        return "100kb-1Mb"
    return ">=1Mb"


def svtype_group(svtype):
    return svtype if svtype in {"DEL", "INS", "INV", "DUP"} else "Other"


def frequency_bin(sample_count, total_samples):
    # Bins adapted for the 33-accession oil-palm panel.
    if sample_count <= 1:
        return "Private (1)"
    if sample_count <= 5:
        return "Low (2-5)"
    if sample_count <= 16:
        return "Intermediate (6-16)"
    if sample_count < total_samples:
        return f"High (17-{total_samples - 1})"
    return f"Core ({total_samples})"


def normalize_repeat_group(raw_name, order, superfamily):
    text = "/".join([raw_name or "", order or "", superfamily or ""]).lower()
    if "gypsy" in text:
        return "LTR/Gypsy"
    if "copia" in text:
        return "LTR/Copia"
    if "ltr" in text:
        return "Other LTR"
    if "helitron" in text:
        return "Helitron"
    if "line" in text:
        return "LINE"
    if "sine" in text:
        return "SINE"
    if any(token in text for token in ("tir", "cacta", "mutator", "tc1", "mariner", "hat", "pile")):
        return "TIR/other DNA"
    return "Other TE"


def read_repeat_bed(path):
    repeats_by_chrom = defaultdict(list)
    if not path or path.lower() == "none":
        return repeats_by_chrom
    with open(path, "r", encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                continue
            chrom = fields[0]
            start = safe_int(fields[1])
            end = safe_int(fields[2])
            if end <= start:
                continue
            raw_name = fields[4] if len(fields) > 4 else "NA"
            order = fields[11] if len(fields) > 11 else "NA"
            superfamily = fields[12] if len(fields) > 12 else raw_name
            repeat_group = normalize_repeat_group(raw_name, order, superfamily)
            repeats_by_chrom[chrom].append(RepeatRecord(start, end, repeat_group, superfamily))
    for chrom in repeats_by_chrom:
        repeats_by_chrom[chrom].sort(key=lambda rec: (rec.start, rec.end))
    return repeats_by_chrom


def cluster_reference_interval(cluster, ins_window):
    if cluster.svtype == "INS":
        start0 = max(0, cluster.rep_start - 1 - ins_window)
        end0 = max(start0 + 1, cluster.rep_start + ins_window)
        return start0, end0
    start = min(cluster.rep_start, cluster.rep_end)
    end = max(cluster.rep_start, cluster.rep_end)
    start0 = max(0, start - 1)
    end0 = max(start0 + 1, end)
    return start0, end0


def merge_segments(segments):
    if not segments:
        return 0
    segments.sort()
    total = 0
    cur_s, cur_e = segments[0]
    for start, end in segments[1:]:
        if start <= cur_e:
            cur_e = max(cur_e, end)
        else:
            total += cur_e - cur_s
            cur_s, cur_e = start, end
    total += cur_e - cur_s
    return total


def annotate_repeat_overlap(clusters, repeats_by_chrom, chrom_order, ins_window):
    clusters_by_chrom = defaultdict(list)
    for cluster in clusters:
        start0, end0 = cluster_reference_interval(cluster, ins_window)
        cluster.interval_start0 = start0
        cluster.interval_end0 = end0
        cluster.interval_len = max(1, end0 - start0)
        clusters_by_chrom[cluster.chrom].append(cluster)

    for chrom in clusters_by_chrom:
        clusters_by_chrom[chrom].sort(key=lambda c: (c.interval_start0, c.interval_end0))

    for chrom in sorted(clusters_by_chrom, key=lambda c: (chrom_rank(c, chrom_order), c)):
        repeats = repeats_by_chrom.get(chrom, [])
        pointer = 0
        active = []
        for cluster in clusters_by_chrom[chrom]:
            start0 = cluster.interval_start0
            end0 = cluster.interval_end0
            while pointer < len(repeats) and repeats[pointer].start < end0:
                active.append(repeats[pointer])
                pointer += 1
            if active:
                active = [rec for rec in active if rec.end > start0]
            segments = []
            group_bp = Counter()
            superfamily_bp = Counter()
            for rec in active:
                if rec.start >= end0 or rec.end <= start0:
                    continue
                overlap = min(end0, rec.end) - max(start0, rec.start)
                if overlap <= 0:
                    continue
                segments.append((max(start0, rec.start), min(end0, rec.end)))
                group_bp[rec.group] += overlap
                superfamily_bp[rec.superfamily] += overlap
            unique_bp = merge_segments(segments)
            cluster.repeat_overlap_bp = unique_bp
            cluster.repeat_overlap_fraction = min(1.0, unique_bp / cluster.interval_len)
            if group_bp:
                cluster.top_repeat_group, cluster.top_repeat_bp = group_bp.most_common(1)[0]
                cluster.top_repeat_superfamily = superfamily_bp.most_common(1)[0][0] if superfamily_bp else "NA"


def quantile(values, prob):
    if not values:
        return 0.0
    ordered = sorted(values)
    pos = prob * (len(ordered) - 1)
    lower = int(pos)
    upper = min(lower + 1, len(ordered) - 1)
    fraction = pos - lower
    return ordered[lower] * (1 - fraction) + ordered[upper] * fraction


def build_saturation(clusters, samples, iterations, seed):
    per_sample_group = {
        group: {sample: set() for sample in samples}
        for group in PLOT_GROUPS
    }
    for idx, cluster in enumerate(clusters):
        group = svtype_group(cluster.svtype)
        for sample in cluster.samples:
            per_sample_group["All"][sample].add(idx)
            per_sample_group[group][sample].add(idx)

    rng = random.Random(seed)
    values = {group: {n: [] for n in range(1, len(samples) + 1)} for group in PLOT_GROUPS}
    for _ in range(iterations):
        order = samples[:]
        rng.shuffle(order)
        for group in PLOT_GROUPS:
            seen = set()
            for n, sample in enumerate(order, 1):
                seen.update(per_sample_group[group][sample])
                values[group][n].append(len(seen))

    rows = []
    for group in PLOT_GROUPS:
        for n in range(1, len(samples) + 1):
            vals = values[group][n]
            rows.append(
                {
                    "SVTYPE_Group": group,
                    "N_Samples": n,
                    "Iterations": iterations,
                    "Mean_Clusters": sum(vals) / len(vals) if vals else 0,
                    "Median_Clusters": statistics.median(vals) if vals else 0,
                    "Lower_95": quantile(vals, 0.025),
                    "Upper_95": quantile(vals, 0.975),
                    "Min_Clusters": min(vals) if vals else 0,
                    "Max_Clusters": max(vals) if vals else 0,
                }
            )
    return rows


def write_tsv(path, rows, fieldnames):
    with open(path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({name: row.get(name, "") for name in fieldnames})


def check_outputs(paths, force):
    existing = [path for path in paths if os.path.exists(path)]
    if existing and not force:
        raise SystemExit("Refusing to overwrite existing outputs:\n" + "\n".join(existing))


def main():
    args = parse_args()
    os.makedirs(args.out_dir, exist_ok=True)
    output_paths = {
        "catalog": os.path.join(args.out_dir, "sv_population_catalog.tsv"),
        "membership": os.path.join(args.out_dir, "sv_cluster_membership.tsv"),
        "frequency": os.path.join(args.out_dir, "sv_frequency_spectrum.tsv"),
        "saturation": os.path.join(args.out_dir, "sv_saturation_curve.tsv"),
        "repeat_type": os.path.join(args.out_dir, "sv_repeat_by_svtype.tsv"),
        "repeat_frequency": os.path.join(args.out_dir, "sv_repeat_by_frequency.tsv"),
        "repeat_group": os.path.join(args.out_dir, "sv_repeat_group_by_svtype.tsv"),
        "qa": os.path.join(args.out_dir, "population_repeat_catalog_qa.md"),
        "metadata": os.path.join(args.out_dir, "population_repeat_catalog_metadata.json"),
    }
    check_outputs(output_paths.values(), args.force)

    chrom_order, chrom_lengths = read_chrom_order(args.fai)
    calls = read_records(args.records)
    if not calls:
        raise SystemExit("No valid SV calls were loaded.")
    samples = sorted({call.sample for call in calls})

    clusters, assignments = cluster_calls(calls, chrom_order, args)
    internal_to_final = {cluster.internal_id: f"SVCL{idx:07d}" for idx, cluster in enumerate(clusters, 1)}
    final_id_to_cluster = {internal_to_final[cluster.internal_id]: cluster for cluster in clusters}

    repeats_by_chrom = read_repeat_bed(args.repeat_bed)
    if repeats_by_chrom:
        annotate_repeat_overlap(clusters, repeats_by_chrom, chrom_order, args.ins_repeat_window)
    else:
        for cluster in clusters:
            start0, end0 = cluster_reference_interval(cluster, args.ins_repeat_window)
            cluster.interval_start0 = start0
            cluster.interval_end0 = end0
            cluster.interval_len = max(1, end0 - start0)

    catalog_rows = []
    for cluster in clusters:
        cluster_id = internal_to_final[cluster.internal_id]
        sample_count = len(cluster.samples)
        group = svtype_group(cluster.svtype)
        caller_combo = cluster.caller_combos.most_common(1)[0][0] if cluster.caller_combos else "NA"
        size_bin = size_bin_from_len(cluster.rep_len)
        catalog_rows.append(
            {
                "Cluster_ID": cluster_id,
                "Chrom": cluster.chrom,
                "Start": cluster.rep_start,
                "End": cluster.rep_end,
                "SVTYPE": cluster.svtype,
                "SVTYPE_Group": group,
                "SVLEN_Median_bp": median_int(cluster.svlens),
                "SVLEN_Min_bp": min(cluster.svlens),
                "SVLEN_Max_bp": max(cluster.svlens),
                "Start_Min": cluster.min_start,
                "Start_Max": cluster.max_start,
                "End_Min": cluster.min_end,
                "End_Max": cluster.max_end,
                "Call_Count": cluster.call_count,
                "Sample_Count": sample_count,
                "Frequency": sample_count / len(samples),
                "Frequency_Bin": frequency_bin(sample_count, len(samples)),
                "Samples": ";".join(sorted(cluster.samples)),
                "Dominant_Caller_Combo": caller_combo,
                "Size_Bin": size_bin,
                "Interval_Start0": cluster.interval_start0,
                "Interval_End0": cluster.interval_end0,
                "Interval_Length_bp": cluster.interval_len,
                "Repeat_Overlap_bp": cluster.repeat_overlap_bp,
                "Repeat_Overlap_Fraction": f"{cluster.repeat_overlap_fraction:.6f}",
                "Dominant_Repeat_Group": cluster.top_repeat_group,
                "Dominant_Repeat_Superfamily": cluster.top_repeat_superfamily,
                "Dominant_Repeat_bp": cluster.top_repeat_bp,
            }
        )

    catalog_fields = [
        "Cluster_ID",
        "Chrom",
        "Start",
        "End",
        "SVTYPE",
        "SVTYPE_Group",
        "SVLEN_Median_bp",
        "SVLEN_Min_bp",
        "SVLEN_Max_bp",
        "Start_Min",
        "Start_Max",
        "End_Min",
        "End_Max",
        "Call_Count",
        "Sample_Count",
        "Frequency",
        "Frequency_Bin",
        "Samples",
        "Dominant_Caller_Combo",
        "Size_Bin",
        "Interval_Start0",
        "Interval_End0",
        "Interval_Length_bp",
        "Repeat_Overlap_bp",
        "Repeat_Overlap_Fraction",
        "Dominant_Repeat_Group",
        "Dominant_Repeat_Superfamily",
        "Dominant_Repeat_bp",
    ]
    write_tsv(output_paths["catalog"], catalog_rows, catalog_fields)

    cluster_sample_counts = {
        internal_to_final[cluster.internal_id]: len(cluster.samples)
        for cluster in clusters
    }
    membership_rows = []
    for call in sorted(calls, key=lambda c: c.idx):
        cluster_id = internal_to_final[assignments[call.idx]]
        membership_rows.append(
            {
                "Cluster_ID": cluster_id,
                "Cluster_Sample_Count": cluster_sample_counts[cluster_id],
                "Sample": call.sample,
                "Chrom": call.chrom,
                "Start": call.start,
                "End": call.end,
                "SVTYPE": call.svtype,
                "SVLEN_bp": call.svlen,
                "SVLEN_raw": call.svlen_raw,
                "SUPP": call.supp,
                "Caller_Combo": call.caller_combo,
                "Size_Bin": call.size_bin,
            }
        )
    write_tsv(
        output_paths["membership"],
        membership_rows,
        [
            "Cluster_ID",
            "Cluster_Sample_Count",
            "Sample",
            "Chrom",
            "Start",
            "End",
            "SVTYPE",
            "SVLEN_bp",
            "SVLEN_raw",
            "SUPP",
            "Caller_Combo",
            "Size_Bin",
        ],
    )

    frequency_counts = {group: Counter() for group in PLOT_GROUPS}
    for cluster in clusters:
        sample_count = len(cluster.samples)
        frequency_counts["All"][sample_count] += 1
        frequency_counts[svtype_group(cluster.svtype)][sample_count] += 1
    frequency_rows = []
    for group in PLOT_GROUPS:
        for n in range(1, len(samples) + 1):
            frequency_rows.append(
                {
                    "SVTYPE_Group": group,
                    "Sample_Count": n,
                    "Frequency": n / len(samples),
                    "Cluster_Count": frequency_counts[group][n],
                }
            )
    write_tsv(output_paths["frequency"], frequency_rows, ["SVTYPE_Group", "Sample_Count", "Frequency", "Cluster_Count"])

    saturation_rows = build_saturation(clusters, samples, args.saturation_iterations, args.seed)
    write_tsv(
        output_paths["saturation"],
        saturation_rows,
        [
            "SVTYPE_Group",
            "N_Samples",
            "Iterations",
            "Mean_Clusters",
            "Median_Clusters",
            "Lower_95",
            "Upper_95",
            "Min_Clusters",
            "Max_Clusters",
        ],
    )

    repeat_by_type = []
    repeat_by_frequency = []
    repeat_group_counts = []
    for group in PLOT_TYPE_GROUPS:
        group_clusters = [cluster for cluster in clusters if svtype_group(cluster.svtype) == group]
        repeat_clusters = [cluster for cluster in group_clusters if cluster.repeat_overlap_bp > 0]
        repeat_by_type.append(
            {
                "SVTYPE_Group": group,
                "Total_Clusters": len(group_clusters),
                "Repeat_Overlapping_Clusters": len(repeat_clusters),
                "Repeat_Overlap_Rate": len(repeat_clusters) / len(group_clusters) if group_clusters else 0,
                "Median_Repeat_Overlap_Fraction": statistics.median(
                    [cluster.repeat_overlap_fraction for cluster in group_clusters]
                )
                if group_clusters
                else 0,
            }
        )
        dominant_counts = Counter(
            cluster.top_repeat_group if cluster.repeat_overlap_bp > 0 else "No repeat"
            for cluster in group_clusters
        )
        for repeat_group, count in sorted(dominant_counts.items()):
            repeat_group_counts.append(
                {
                    "SVTYPE_Group": group,
                    "Repeat_Group": repeat_group,
                    "Cluster_Count": count,
                    "Fraction": count / len(group_clusters) if group_clusters else 0,
                }
            )

        bin_order = ["Private (1)", "Low (2-5)", "Intermediate (6-16)", f"High (17-{len(samples) - 1})", f"Core ({len(samples)})"]
        for freq_bin in bin_order:
            bin_clusters = [
                cluster
                for cluster in group_clusters
                if frequency_bin(len(cluster.samples), len(samples)) == freq_bin
            ]
            bin_repeat = [cluster for cluster in bin_clusters if cluster.repeat_overlap_bp > 0]
            repeat_by_frequency.append(
                {
                    "SVTYPE_Group": group,
                    "Frequency_Bin": freq_bin,
                    "Total_Clusters": len(bin_clusters),
                    "Repeat_Overlapping_Clusters": len(bin_repeat),
                    "Repeat_Overlap_Rate": len(bin_repeat) / len(bin_clusters) if bin_clusters else 0,
                }
            )

    write_tsv(
        output_paths["repeat_type"],
        repeat_by_type,
        [
            "SVTYPE_Group",
            "Total_Clusters",
            "Repeat_Overlapping_Clusters",
            "Repeat_Overlap_Rate",
            "Median_Repeat_Overlap_Fraction",
        ],
    )
    write_tsv(
        output_paths["repeat_frequency"],
        repeat_by_frequency,
        [
            "SVTYPE_Group",
            "Frequency_Bin",
            "Total_Clusters",
            "Repeat_Overlapping_Clusters",
            "Repeat_Overlap_Rate",
        ],
    )
    write_tsv(
        output_paths["repeat_group"],
        repeat_group_counts,
        ["SVTYPE_Group", "Repeat_Group", "Cluster_Count", "Fraction"],
    )

    svtype_cluster_counts = Counter(cluster.svtype for cluster in clusters)
    repeat_total = sum(1 for cluster in clusters if cluster.repeat_overlap_bp > 0)
    metadata = {
        "records": os.path.abspath(args.records),
        "repeat_bed": os.path.abspath(args.repeat_bed) if args.repeat_bed.lower() != "none" else "none",
        "fai": os.path.abspath(args.fai) if args.fai else "",
        "out_dir": os.path.abspath(args.out_dir),
        "input_call_count": len(calls),
        "cluster_count": len(clusters),
        "sample_count": len(samples),
        "samples": samples,
        "svtype_cluster_counts": dict(svtype_cluster_counts),
        "repeat_overlapping_clusters": repeat_total,
        "repeat_overlap_rate": repeat_total / len(clusters) if clusters else 0,
        "merge_distance": args.merge_distance,
        "ins_distance": args.ins_distance,
        "min_reciprocal_overlap": args.min_reciprocal_overlap,
        "min_size_similarity": args.min_size_similarity,
        "ins_repeat_window": args.ins_repeat_window,
        "saturation_iterations": args.saturation_iterations,
        "seed": args.seed,
        "chrom_lengths_loaded": len(chrom_lengths),
    }
    with open(output_paths["metadata"], "w", encoding="utf-8") as handle:
        json.dump(metadata, handle, indent=2, sort_keys=True)
        handle.write("\n")

    qa_lines = [
        "# SUPP>=2 Population SV Catalog QA",
        "",
        f"- Input records: {len(calls):,}",
        f"- Samples: {len(samples)}",
        f"- Nonredundant SV clusters: {len(clusters):,}",
        f"- Repeat-overlapping clusters: {repeat_total:,} ({metadata['repeat_overlap_rate']:.2%})",
        f"- Merge distance: {args.merge_distance} bp for span SVs; {args.ins_distance} bp for INS",
        f"- Minimum reciprocal overlap: {args.min_reciprocal_overlap}",
        f"- Minimum size similarity: {args.min_size_similarity}",
        f"- INS repeat context window: +/-{args.ins_repeat_window} bp around reference breakpoint",
        "",
        "## Cluster Counts By SVTYPE",
    ]
    for svtype in SVTYPE_ORDER:
        qa_lines.append(f"- {svtype}: {svtype_cluster_counts.get(svtype, 0):,}")
    qa_lines.extend(
        [
            "",
            "## Caveats",
            "- This catalog is a coordinate-based cross-sample clustering of per-sample SUPP>=2 calls; it is intended for population-frequency and accumulation summaries.",
            "- INS repeat overlap uses the reference breakpoint window, not the inserted sequence itself.",
            "- TRA/HDR/INVDP/TDM are grouped as Other in the plotting summaries because they are sparse and length semantics are not directly comparable.",
        ]
    )
    with open(output_paths["qa"], "w", encoding="utf-8") as handle:
        handle.write("\n".join(qa_lines) + "\n")


if __name__ == "__main__":
    main()
