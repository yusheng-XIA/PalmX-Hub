#!/usr/bin/env python3
"""Build All38/African35 load-only ideal parental haplotypes by exact DP."""

import argparse
import csv
import itertools
import math
import os
import re
from collections import Counter, defaultdict
from pathlib import Path


YES_VALUES = {"yes", "true", "1"}
NA_VALUES = {"", "na", "nan", "none"}


def parse_args():
    parser = argparse.ArgumentParser(
        description="Build W=0 frequency-inclusive load-only ideal parental haplotypes."
    )
    parser.add_argument("--dsnp", required=True)
    parser.add_argument("--dsv-all38", required=True)
    parser.add_argument("--dsv-african35", required=True)
    parser.add_argument("--sample-manifest", required=True)
    parser.add_argument("--reference-fai", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--window-size", type=int, default=500_000)
    parser.add_argument("--primary-penalty", type=int, default=15)
    parser.add_argument("--penalties", default="0,5,15,40,100,300")
    parser.add_argument("--expected-all38-samples", type=int, default=38)
    parser.add_argument("--expected-african35-samples", type=int, default=35)
    parser.add_argument("--expected-all38-dsv", type=int, default=1480)
    parser.add_argument("--expected-african35-dsv", type=int, default=924)
    return parser.parse_args()


def as_int(value, default=None):
    try:
        if value is None or str(value).strip().lower() in NA_VALUES:
            return default
        return int(float(str(value)))
    except (TypeError, ValueError):
        return default


def parse_samples(value):
    if value is None:
        return []
    samples = []
    for token in re.split(r"[;,|]", str(value)):
        sample = token.strip()
        if sample and sample.lower() not in NA_VALUES:
            samples.append(sample)
    return sorted(set(samples))


def chrom_sort_key(chrom):
    match = re.fullmatch(r"chr(\d+)B", chrom or "")
    if match:
        return int(match.group(1))
    return 10**9


def read_tsv(path):
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = [dict(row) for row in reader]
        fields = reader.fieldnames or []
    if not fields:
        raise ValueError(f"Missing TSV header: {path}")
    return fields, rows


def require_fields(fields, required, path):
    missing = [field for field in required if field not in fields]
    if missing:
        raise ValueError(f"Missing columns in {path}: {','.join(missing)}")


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


def write_text(path, text):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        raise FileExistsError(f"Refusing to overwrite existing output: {path}")
    tmp = path.with_name(path.name + f".tmp.{os.getpid()}")
    with open(tmp, "w") as handle:
        handle.write(text)
    os.replace(tmp, path)


def read_manifest(path):
    fields, rows = read_tsv(path)
    require_fields(fields, ["Sample_ID", "Include_All38", "Include_African35", "Status"], path)
    if len(rows) != len({row["Sample_ID"] for row in rows}):
        raise ValueError("Duplicate Sample_ID in manifest")
    metadata = {}
    for row in rows:
        sample = row["Sample_ID"].strip()
        if row.get("Status", "").strip().upper() != "PASS":
            raise ValueError(f"Manifest status is not PASS: {sample}")
        metadata[sample] = {
            "Include_All38": row["Include_All38"].strip().lower() in YES_VALUES,
            "Include_African35": row["Include_African35"].strip().lower() in YES_VALUES,
            "Population_Scope": row.get("Population_Scope", "NA"),
        }
    return metadata


def read_fai(path, window_size):
    lengths = {}
    with open(path) as handle:
        for line in handle:
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            length = as_int(fields[1] if len(fields) > 1 else None)
            if length is None or length <= 0:
                raise ValueError(f"Invalid FAI line: {line.rstrip()}")
            lengths[fields[0]] = length
    expected = {f"chr{i:02d}B" for i in range(1, 17)}
    if set(lengths) != expected:
        raise ValueError("Reference chromosome set mismatch")
    chromosomes = sorted(lengths, key=chrom_sort_key)
    windows = {
        chrom: (lengths[chrom] + window_size - 1) // window_size for chrom in chromosomes
    }
    return lengths, chromosomes, windows


def dsnp_position0(row):
    pos1 = as_int(row.get("Pos"))
    return None if pos1 is None else pos1 - 1


def dsv_position0(row):
    start0 = as_int(row.get("Interval_Start0"))
    end0 = as_int(row.get("Interval_End0"))
    if start0 is not None and end0 is not None and end0 > start0:
        return (start0 + end0 - 1) // 2
    start1 = as_int(row.get("Start"))
    end1 = as_int(row.get("End"))
    if start1 is not None and end1 is not None:
        return max(0, ((start1 + end1) // 2) - 1)
    pos1 = as_int(row.get("Pos"))
    return None if pos1 is None else max(0, pos1 - 1)


def load_dsnp(path, all_samples, chrom_lengths):
    fields, rows = read_tsv(path)
    require_fields(
        fields,
        ["SNP_ID", "Chrom", "Pos", "Sample_Count", "Samples", "ALT_Polarity_Phoenix"],
        path,
    )
    events = []
    errors = []
    for row_number, row in enumerate(rows, 2):
        if row.get("ALT_Polarity_Phoenix") != "ALT_Derived":
            errors.append(f"row{row_number}:polarity={row.get('ALT_Polarity_Phoenix')}")
            continue
        chrom = row.get("Chrom", "")
        pos0 = dsnp_position0(row)
        carriers = parse_samples(row.get("Samples"))
        declared = as_int(row.get("Sample_Count"))
        if (
            chrom not in chrom_lengths
            or pos0 is None
            or not (0 <= pos0 < chrom_lengths.get(chrom, 0))
            or not carriers
            or declared != len(carriers)
            or set(carriers) - all_samples
        ):
            errors.append(
                f"row{row_number}:{chrom}:{pos0}:declared={declared}:parsed={len(carriers)}"
            )
            continue
        events.append((row["SNP_ID"], chrom, pos0, tuple(carriers)))
    if errors:
        raise ValueError("Invalid ALT-derived dSNP rows: " + ";".join(errors[:10]))
    return events


def load_dsv(path, carrier_column, panel_samples, chrom_lengths):
    fields, rows = read_tsv(path)
    require_fields(fields, ["SV_Key", "Chrom", carrier_column], path)
    events = []
    errors = []
    for row_number, row in enumerate(rows, 2):
        chrom = row.get("Chrom", "")
        pos0 = dsv_position0(row)
        carriers = parse_samples(row.get(carrier_column))
        if (
            len(carriers) != 1
            or carriers[0] not in panel_samples
            or chrom not in chrom_lengths
            or pos0 is None
            or not (0 <= pos0 < chrom_lengths.get(chrom, 0))
        ):
            errors.append(f"row{row_number}:{chrom}:{pos0}:carriers={','.join(carriers)}")
            continue
        events.append((row["SV_Key"], chrom, pos0, carriers[0]))
    if errors:
        raise ValueError(f"Invalid dSV rows in {path}: " + ";".join(errors[:10]))
    return events


def exact_dp(local_load, donors, penalty):
    """Return exact donor-index path and objective for windows x donors integer costs."""
    if not local_load or not donors:
        raise ValueError("DP requires non-empty windows and donors")
    n_windows = len(local_load)
    n_donors = len(donors)
    if any(len(row) != n_donors for row in local_load):
        raise ValueError("DP matrix width mismatch")
    previous = [int(value) for value in local_load[0]]
    back = [[-1] * n_donors for _ in range(n_windows)]
    for window in range(1, n_windows):
        current = [0] * n_donors
        for donor in range(n_donors):
            predecessor = min(
                range(n_donors),
                key=lambda prev: (
                    previous[prev] + (0 if prev == donor else penalty),
                    0 if prev == donor else 1,
                    donors[prev],
                ),
            )
            current[donor] = (
                previous[predecessor]
                + (0 if predecessor == donor else penalty)
                + int(local_load[window][donor])
            )
            back[window][donor] = predecessor
        previous = current
    end = min(range(n_donors), key=lambda donor: (previous[donor], donors[donor]))
    path = [0] * n_windows
    path[-1] = end
    for window in range(n_windows - 1, 0, -1):
        path[window - 1] = back[window][path[window]]
    return path, previous[end]


def path_objective(local_load, path, penalty):
    objective = 0
    for window, donor in enumerate(path):
        objective += int(local_load[window][donor])
        if window and donor != path[window - 1]:
            objective += penalty
    return objective


def run_dp_self_tests():
    tests = []
    donors = ["A", "B", "C"]
    matrices = [
        [[1, 3, 2], [2, 0, 4], [3, 1, 0], [0, 4, 2]],
        [[0, 0, 1], [0, 0, 1], [2, 1, 0]],
    ]
    for matrix_index, matrix in enumerate(matrices, 1):
        for penalty in (0, 1, 2, 5):
            path, objective = exact_dp(matrix, donors, penalty)
            brute = min(
                path_objective(matrix, candidate, penalty)
                for candidate in itertools.product(range(len(donors)), repeat=len(matrix))
            )
            passed = objective == brute and objective == path_objective(matrix, path, penalty)
            tests.append(
                {
                    "Test": f"Synthetic_{matrix_index}_P{penalty}",
                    "Status": "PASS" if passed else "FAIL",
                    "DP_Objective": objective,
                    "Exhaustive_Objective": brute,
                }
            )
    return tests


def initialize_counts(samples, chromosomes, window_counts):
    return {
        sample: {
            chrom: [0] * window_counts[chrom]
            for chrom in chromosomes
        }
        for sample in samples
    }


def reconstruct_segments(path, donors, chrom, chrom_length, window_size, dsv, dsnp):
    segments = []
    start = 0
    current = path[0]
    for window in range(1, len(path) + 1):
        changed = window == len(path) or path[window] != current
        if not changed:
            continue
        donor = donors[current]
        end_bp = min(window * window_size, chrom_length)
        dsv_count = sum(dsv[donor][chrom][start:window])
        dsnp_count = sum(dsnp[donor][chrom][start:window])
        segments.append(
            {
                "Chrom": chrom,
                "Start_Window_Index": start,
                "End_Window_Index_Exclusive": window,
                "Segment_Start_0based": start * window_size,
                "Segment_End_0based": end_bp,
                "Donor_ID": donor,
                "Window_Count": window - start,
                "Segment_Length_bp": end_bp - start * window_size,
                "DSV_Count": dsv_count,
                "DSNP_Count": dsnp_count,
                "Total_Load": dsv_count + dsnp_count,
            }
        )
        if window < len(path):
            start = window
            current = path[window]
    return segments


def add_check(checks, check, passed, observed, expected):
    checks.append(
        {
            "Check": check,
            "Status": "PASS" if passed else "FAIL",
            "Observed": observed,
            "Expected": expected,
        }
    )


def build_panel(
    panel_label,
    samples,
    dsnp_events,
    dsv_events,
    expected_dsv,
    chrom_lengths,
    chromosomes,
    window_counts,
    window_size,
    penalties,
    primary_penalty,
    outroot,
    self_tests,
):
    panel_key = panel_label.lower()
    panel_dir = Path(outroot) / panel_key
    panel_dir.mkdir(parents=True, exist_ok=False)
    qa_dir = panel_dir / "qa"
    qa_dir.mkdir()
    sample_set = set(samples)

    dsnp_counts = initialize_counts(samples, chromosomes, window_counts)
    dsv_counts = initialize_counts(samples, chromosomes, window_counts)
    selected_dsnp_sites = 0
    fixed_dsnp_sites = 0
    no_panel_carrier = 0
    for _, chrom, pos0, carriers in dsnp_events:
        selected = sorted(set(carriers) & sample_set)
        if not selected:
            no_panel_carrier += 1
            continue
        selected_dsnp_sites += 1
        if set(selected) == sample_set:
            fixed_dsnp_sites += 1
        window = min(pos0 // window_size, window_counts[chrom] - 1)
        for sample in selected:
            dsnp_counts[sample][chrom][window] += 1
    for _, chrom, pos0, sample in dsv_events:
        window = min(pos0 // window_size, window_counts[chrom] - 1)
        dsv_counts[sample][chrom][window] += 1

    matrix_rows = []
    for sample in samples:
        for chrom in chromosomes:
            for window in range(window_counts[chrom]):
                start = window * window_size
                end = min(start + window_size, chrom_lengths[chrom])
                dsv_value = dsv_counts[sample][chrom][window]
                dsnp_value = dsnp_counts[sample][chrom][window]
                matrix_rows.append(
                    {
                        "Panel": panel_label,
                        "Sample_ID": sample,
                        "Chrom": chrom,
                        "Window_Index": window,
                        "Window_Start_0based": start,
                        "Window_End_0based": end,
                        "DSV_Count": dsv_value,
                        "DSNP_Count": dsnp_value,
                        "Total_Load": dsv_value + dsnp_value,
                    }
                )
    matrix_fields = [
        "Panel", "Sample_ID", "Chrom", "Window_Index", "Window_Start_0based",
        "Window_End_0based", "DSV_Count", "DSNP_Count", "Total_Load",
    ]
    write_tsv(panel_dir / "load_matrix_500kb.tsv", matrix_rows, matrix_fields)

    sweep_rows = []
    sweep_paths = {}
    sweep_chrom_stats = {}
    for penalty in penalties:
        total_dsv = total_dsnp = total_breakpoints = total_segments = total_objective = 0
        used_donors = set()
        paths = {}
        chrom_stats = {}
        for chrom in chromosomes:
            local_dsv = [
                [dsv_counts[sample][chrom][window] for sample in samples]
                for window in range(window_counts[chrom])
            ]
            local_dsnp = [
                [dsnp_counts[sample][chrom][window] for sample in samples]
                for window in range(window_counts[chrom])
            ]
            local_total = [
                [local_dsv[w][d] + local_dsnp[w][d] for d in range(len(samples))]
                for w in range(window_counts[chrom])
            ]
            path, objective = exact_dp(local_total, samples, penalty)
            if objective != path_objective(local_total, path, penalty):
                raise RuntimeError(f"DP objective recount mismatch: {panel_label} {chrom} P={penalty}")
            residual_dsv = sum(local_dsv[w][path[w]] for w in range(len(path)))
            residual_dsnp = sum(local_dsnp[w][path[w]] for w in range(len(path)))
            breakpoints = sum(path[w] != path[w - 1] for w in range(1, len(path)))
            segments = breakpoints + 1
            donor_names = {samples[index] for index in path}
            paths[chrom] = path
            chrom_stats[chrom] = {
                "Residual_DSV": residual_dsv,
                "Residual_DSNP": residual_dsnp,
                "Breakpoint_Count": breakpoints,
                "Segment_Count": segments,
                "Donor_Count": len(donor_names),
                "Objective": objective,
            }
            total_dsv += residual_dsv
            total_dsnp += residual_dsnp
            total_breakpoints += breakpoints
            total_segments += segments
            total_objective += objective
            used_donors.update(donor_names)
        sweep_rows.append(
            {
                "Panel": panel_label,
                "Breakpoint_Penalty_P": penalty,
                "Residual_DSV": total_dsv,
                "Residual_DSNP": total_dsnp,
                "Residual_Total_Load": total_dsv + total_dsnp,
                "Breakpoint_Count": total_breakpoints,
                "Segment_Count": total_segments,
                "Donor_Count": len(used_donors),
                "Objective": total_objective,
                "GWAS_Weight": 0,
            }
        )
        sweep_paths[penalty] = paths
        sweep_chrom_stats[penalty] = chrom_stats
    sweep_fields = [
        "Panel", "Breakpoint_Penalty_P", "Residual_DSV", "Residual_DSNP",
        "Residual_Total_Load", "Breakpoint_Count", "Segment_Count", "Donor_Count",
        "Objective", "GWAS_Weight",
    ]
    write_tsv(panel_dir / "penalty_sweep.tsv", sweep_rows, sweep_fields)

    primary_paths = sweep_paths[primary_penalty]
    primary_stats = sweep_chrom_stats[primary_penalty]
    path_rows = []
    segment_rows = []
    by_chrom_rows = []
    donor_contribution = defaultdict(Counter)
    cumulative_genome = 0
    for chrom in chromosomes:
        path = primary_paths[chrom]
        cumulative_chrom = 0
        previous = None
        for window, donor_index in enumerate(path):
            donor = samples[donor_index]
            switched = previous is not None and donor_index != previous
            dsv_value = dsv_counts[donor][chrom][window]
            dsnp_value = dsnp_counts[donor][chrom][window]
            local = dsv_value + dsnp_value
            cumulative_chrom += local + (primary_penalty if switched else 0)
            cumulative_genome += local + (primary_penalty if switched else 0)
            start = window * window_size
            end = min(start + window_size, chrom_lengths[chrom])
            path_rows.append(
                {
                    "Panel": panel_label,
                    "Chrom": chrom,
                    "Window_Index": window,
                    "Window_Start_0based": start,
                    "Window_End_0based": end,
                    "Donor_ID": donor,
                    "DSV_Count": dsv_value,
                    "DSNP_Count": dsnp_value,
                    "Total_Load": local,
                    "Donor_Switch_From_Previous": "Yes" if switched else "No",
                    "Cumulative_Chrom_Objective": cumulative_chrom,
                    "Cumulative_Genome_Objective": cumulative_genome,
                    "Breakpoint_Penalty_P": primary_penalty,
                    "GWAS_Weight": 0,
                }
            )
            donor_contribution[donor]["Window_Count"] += 1
            donor_contribution[donor]["Selected_Length_bp"] += end - start
            donor_contribution[donor]["Selected_DSV"] += dsv_value
            donor_contribution[donor]["Selected_DSNP"] += dsnp_value
            previous = donor_index
        segments = reconstruct_segments(
            path, samples, chrom, chrom_lengths[chrom], window_size, dsv_counts, dsnp_counts
        )
        for segment_number, row in enumerate(segments, 1):
            row.update(
                {
                    "Panel": panel_label,
                    "Segment_ID": f"{chrom}_SEG{segment_number:04d}",
                    "Breakpoint_Penalty_P": primary_penalty,
                    "GWAS_Weight": 0,
                }
            )
            segment_rows.append(row)
            donor_contribution[row["Donor_ID"]]["Segment_Count"] += 1
        stat = primary_stats[chrom]
        by_chrom_rows.append(
            {
                "Panel": panel_label,
                "Chrom": chrom,
                "Window_Count": window_counts[chrom],
                **stat,
                "Residual_Total_Load": stat["Residual_DSV"] + stat["Residual_DSNP"],
                "Breakpoint_Penalty_P": primary_penalty,
                "GWAS_Weight": 0,
            }
        )

    total_stat = {
        key: sum(int(row[key]) for row in by_chrom_rows)
        for key in (
            "Window_Count", "Residual_DSV", "Residual_DSNP", "Residual_Total_Load",
            "Breakpoint_Count", "Segment_Count", "Objective",
        )
    }
    total_stat["Donor_Count"] = len({row["Donor_ID"] for row in path_rows})
    by_chrom_rows.append(
        {
            "Panel": panel_label,
            "Chrom": "TOTAL",
            **total_stat,
            "Breakpoint_Penalty_P": primary_penalty,
            "GWAS_Weight": 0,
        }
    )
    by_chrom_fields = [
        "Panel", "Chrom", "Window_Count", "Residual_DSV", "Residual_DSNP",
        "Residual_Total_Load", "Breakpoint_Count", "Segment_Count", "Donor_Count",
        "Objective", "Breakpoint_Penalty_P", "GWAS_Weight",
    ]
    write_tsv(panel_dir / "ideal_loadonly_by_chrom.tsv", by_chrom_rows, by_chrom_fields)
    path_fields = [
        "Panel", "Chrom", "Window_Index", "Window_Start_0based", "Window_End_0based",
        "Donor_ID", "DSV_Count", "DSNP_Count", "Total_Load",
        "Donor_Switch_From_Previous", "Cumulative_Chrom_Objective",
        "Cumulative_Genome_Objective", "Breakpoint_Penalty_P", "GWAS_Weight",
    ]
    write_tsv(panel_dir / "ideal_loadonly_path.tsv", path_rows, path_fields)
    segment_fields = [
        "Panel", "Segment_ID", "Chrom", "Start_Window_Index", "End_Window_Index_Exclusive",
        "Segment_Start_0based", "Segment_End_0based", "Donor_ID", "Window_Count",
        "Segment_Length_bp", "DSV_Count", "DSNP_Count", "Total_Load",
        "Breakpoint_Penalty_P", "GWAS_Weight",
    ]
    write_tsv(panel_dir / "ideal_loadonly_segments.tsv", segment_rows, segment_fields)

    genome_bp = sum(chrom_lengths.values())
    donor_rows = []
    for donor, counts in donor_contribution.items():
        donor_rows.append(
            {
                "Panel": panel_label,
                "Donor_ID": donor,
                "Window_Count": counts["Window_Count"],
                "Segment_Count": counts["Segment_Count"],
                "Selected_Length_bp": counts["Selected_Length_bp"],
                "Selected_Genome_Pct": f"{100 * counts['Selected_Length_bp'] / genome_bp:.6f}",
                "Selected_DSV": counts["Selected_DSV"],
                "Selected_DSNP": counts["Selected_DSNP"],
                "Selected_Total_Load": counts["Selected_DSV"] + counts["Selected_DSNP"],
            }
        )
    donor_rows.sort(key=lambda row: (-int(row["Window_Count"]), row["Donor_ID"]))
    write_tsv(
        panel_dir / "donor_contribution.tsv",
        donor_rows,
        [
            "Panel", "Donor_ID", "Window_Count", "Segment_Count", "Selected_Length_bp",
            "Selected_Genome_Pct", "Selected_DSV", "Selected_DSNP", "Selected_Total_Load",
        ],
    )

    fixed_row = {
        "Panel": panel_label,
        "Sample_Count": len(samples),
        "Source_ALT_Derived_Site_Count": len(dsnp_events),
        "Panel_Derived_Site_Count": selected_dsnp_sites,
        "Variable_Derived_Site_Count": selected_dsnp_sites - fixed_dsnp_sites,
        "Fixed_Unavoidable_Derived_Site_Count": fixed_dsnp_sites,
        "No_Panel_Carrier_Site_Count": no_panel_carrier,
    }
    write_tsv(
        panel_dir / "fixed_unavoidable_summary.tsv",
        [fixed_row],
        list(fixed_row),
    )

    donor_total_load = {}
    for sample in samples:
        donor_total_load[sample] = sum(
            sum(dsnp_counts[sample][chrom]) + sum(dsv_counts[sample][chrom])
            for chrom in chromosomes
        )
    best_donor = min(samples, key=lambda sample: (donor_total_load[sample], sample))
    best_donor_load = donor_total_load[best_donor]
    residual_total = total_stat["Residual_Total_Load"]
    reduction = best_donor_load - residual_total
    reduction_pct = 100 * reduction / best_donor_load if best_donor_load else math.nan
    informative_chrom = min(
        chromosomes,
        key=lambda chrom: (
            -int(primary_stats[chrom]["Breakpoint_Count"]),
            -(int(primary_stats[chrom]["Residual_DSV"]) + int(primary_stats[chrom]["Residual_DSNP"])),
            chrom_sort_key(chrom),
        ),
    )
    panel_summary = {
        "Panel": panel_label,
        "Sample_Count": len(samples),
        "Chromosome_Count": len(chromosomes),
        "Window_Count": sum(window_counts.values()),
        "Frequency_Inclusive_ALT_Derived_Site_Count": selected_dsnp_sites,
        "Fixed_Unavoidable_Derived_Site_Count": fixed_dsnp_sites,
        "Formal_DSV_Count": len(dsv_events),
        "Primary_Breakpoint_Penalty_P": primary_penalty,
        "GWAS_Weight": 0,
        "Residual_DSV": total_stat["Residual_DSV"],
        "Residual_DSNP": total_stat["Residual_DSNP"],
        "Residual_Total_Load": residual_total,
        "Breakpoint_Count": total_stat["Breakpoint_Count"],
        "Segment_Count": total_stat["Segment_Count"],
        "Donor_Count": total_stat["Donor_Count"],
        "Best_Single_Donor": best_donor,
        "Best_Single_Donor_Load": best_donor_load,
        "Mosaic_Load_Reduction": reduction,
        "Mosaic_Load_Reduction_Pct": "NA" if not math.isfinite(reduction_pct) else f"{reduction_pct:.6f}",
        "Informative_Zoom_Chrom": informative_chrom,
    }
    write_tsv(panel_dir / "panel_summary.tsv", [panel_summary], list(panel_summary))

    expanded = []
    for segment in segment_rows:
        expanded.extend(
            (segment["Chrom"], window, segment["Donor_ID"])
            for window in range(
                int(segment["Start_Window_Index"]), int(segment["End_Window_Index_Exclusive"])
            )
        )
    expected_expanded = [
        (row["Chrom"], int(row["Window_Index"]), row["Donor_ID"])
        for row in path_rows
    ]
    p0_row = next(row for row in sweep_rows if int(row["Breakpoint_Penalty_P"]) == 0)
    p0_minimum = 0
    for chrom in chromosomes:
        for window in range(window_counts[chrom]):
            p0_minimum += min(
                dsv_counts[sample][chrom][window] + dsnp_counts[sample][chrom][window]
                for sample in samples
            )
    recount_dsv = sum(int(row["DSV_Count"]) for row in path_rows)
    recount_dsnp = sum(int(row["DSNP_Count"]) for row in path_rows)
    recount_breaks = sum(row["Donor_Switch_From_Previous"] == "Yes" for row in path_rows)
    checks = []
    add_check(checks, "DP_synthetic_exhaustive_tests", all(t["Status"] == "PASS" for t in self_tests), sum(t["Status"] == "PASS" for t in self_tests), len(self_tests))
    add_check(checks, "Panel_sample_count", len(samples) in (38, 35), len(samples), "38_or_35")
    add_check(checks, "Reference_chromosome_count", len(chromosomes) == 16, len(chromosomes), 16)
    add_check(checks, "Reference_window_count", sum(window_counts.values()) == 3486, sum(window_counts.values()), 3486)
    add_check(checks, "Formal_DSV_count", len(dsv_events) == expected_dsv, len(dsv_events), expected_dsv)
    add_check(checks, "Frequency_inclusive_dSNP_nonempty", selected_dsnp_sites > 0, selected_dsnp_sites, ">0")
    add_check(checks, "Load_matrix_row_count", len(matrix_rows) == len(samples) * 3486, len(matrix_rows), len(samples) * 3486)
    add_check(checks, "Primary_path_row_count", len(path_rows) == 3486, len(path_rows), 3486)
    add_check(checks, "Path_donor_membership", all(row["Donor_ID"] in sample_set for row in path_rows), True, True)
    add_check(checks, "Segment_path_reconstruction", expanded == expected_expanded, True, True)
    add_check(checks, "Residual_DSV_recount", recount_dsv == total_stat["Residual_DSV"], recount_dsv, total_stat["Residual_DSV"])
    add_check(checks, "Residual_DSNP_recount", recount_dsnp == total_stat["Residual_DSNP"], recount_dsnp, total_stat["Residual_DSNP"])
    add_check(checks, "Breakpoint_recount", recount_breaks == total_stat["Breakpoint_Count"], recount_breaks, total_stat["Breakpoint_Count"])
    add_check(checks, "P0_equals_per_window_minimum", int(p0_row["Residual_Total_Load"]) == p0_minimum, p0_row["Residual_Total_Load"], p0_minimum)
    add_check(checks, "Residual_DSNP_above_fixed_floor", total_stat["Residual_DSNP"] >= fixed_dsnp_sites, total_stat["Residual_DSNP"], f">={fixed_dsnp_sites}")
    add_check(checks, "GWAS_weight_zero", all(int(row["GWAS_Weight"]) == 0 for row in sweep_rows + path_rows + by_chrom_rows), 0, 0)
    add_check(checks, "Penalty_set_exact", [int(row["Breakpoint_Penalty_P"]) for row in sweep_rows] == penalties, ",".join(str(row["Breakpoint_Penalty_P"]) for row in sweep_rows), ",".join(map(str, penalties)))
    write_tsv(
        qa_dir / "final_integrity_summary.tsv",
        checks,
        ["Check", "Status", "Observed", "Expected"],
    )
    if any(row["Status"] == "FAIL" for row in checks):
        raise RuntimeError(f"Panel integrity failure: {panel_label}")
    return panel_summary, checks, path_rows


def main():
    args = parse_args()
    penalties = [int(value) for value in args.penalties.split(",") if value.strip()]
    if penalties != [0, 5, 15, 40, 100, 300]:
        raise ValueError("Reviewed formal penalty set must be exactly 0,5,15,40,100,300")
    if args.primary_penalty != 15 or args.primary_penalty not in penalties:
        raise ValueError("Reviewed primary breakpoint penalty must be 15")
    outroot = Path(args.outdir)
    for target in (outroot / "all38", outroot / "african35"):
        if target.exists():
            raise FileExistsError(f"Refusing to overwrite existing panel output: {target}")

    metadata = read_manifest(args.sample_manifest)
    all38 = sorted(sample for sample, row in metadata.items() if row["Include_All38"])
    african35 = sorted(sample for sample, row in metadata.items() if row["Include_African35"])
    if len(all38) != args.expected_all38_samples or len(african35) != args.expected_african35_samples:
        raise ValueError(f"Panel sample mismatch: All38={len(all38)} African35={len(african35)}")
    chrom_lengths, chromosomes, window_counts = read_fai(args.reference_fai, args.window_size)
    dsnp_events = load_dsnp(args.dsnp, set(all38), chrom_lengths)
    dsv_all38 = load_dsv(args.dsv_all38, "Samples", set(all38), chrom_lengths)
    dsv_african35 = load_dsv(
        args.dsv_african35, "Samples_African35", set(african35), chrom_lengths
    )
    self_tests = run_dp_self_tests()
    if any(row["Status"] != "PASS" for row in self_tests):
        raise RuntimeError("DP synthetic exhaustive tests failed")
    write_tsv(
        outroot / "dp_unit_tests.tsv",
        self_tests,
        ["Test", "Status", "DP_Objective", "Exhaustive_Objective"],
    )

    all_summary, all_checks, all_path = build_panel(
        "All38", all38, dsnp_events, dsv_all38, args.expected_all38_dsv,
        chrom_lengths, chromosomes, window_counts, args.window_size, penalties,
        args.primary_penalty, outroot, self_tests,
    )
    afr_summary, afr_checks, afr_path = build_panel(
        "African35", african35, dsnp_events, dsv_african35, args.expected_african35_dsv,
        chrom_lengths, chromosomes, window_counts, args.window_size, penalties,
        args.primary_penalty, outroot, self_tests,
    )

    all_map = {(row["Chrom"], int(row["Window_Index"])): row["Donor_ID"] for row in all_path}
    afr_map = {(row["Chrom"], int(row["Window_Index"])): row["Donor_ID"] for row in afr_path}
    common = sorted(set(all_map) & set(afr_map))
    agreement = sum(all_map[key] == afr_map[key] for key in common)
    comparison = {
        "Primary_Panel": "All38",
        "Sensitivity_Panel": "African35",
        "Window_Count": len(common),
        "Same_Donor_Window_Count": agreement,
        "Same_Donor_Window_Pct": f"{100 * agreement / len(common):.6f}",
        "All38_Residual_DSV": all_summary["Residual_DSV"],
        "African35_Residual_DSV": afr_summary["Residual_DSV"],
        "All38_Residual_DSNP": all_summary["Residual_DSNP"],
        "African35_Residual_DSNP": afr_summary["Residual_DSNP"],
        "All38_Breakpoint_Count": all_summary["Breakpoint_Count"],
        "African35_Breakpoint_Count": afr_summary["Breakpoint_Count"],
        "GWAS_Weight": 0,
    }
    write_tsv(outroot / "panel_comparison.tsv", [comparison], list(comparison))
    write_tsv(outroot / "panel_summary.tsv", [all_summary, afr_summary], list(all_summary))

    combined_checks = []
    for panel, checks in (("All38", all_checks), ("African35", afr_checks)):
        combined_checks.extend({"Panel": panel, **row} for row in checks)
    write_tsv(
        outroot / "final_integrity_summary.tsv",
        combined_checks,
        ["Panel", "Check", "Status", "Observed", "Expected"],
    )
    parameter_rows = [
        {"Parameter": "Analysis_Mode", "Value": "Load_only_IPH", "Status": "Locked"},
        {"Parameter": "GWAS_Weight", "Value": 0, "Status": "Locked"},
        {"Parameter": "Window_Size_bp", "Value": args.window_size, "Status": "Locked"},
        {"Parameter": "Primary_Breakpoint_Penalty_P", "Value": args.primary_penalty, "Status": "Locked"},
        {"Parameter": "Penalty_Sweep", "Value": ",".join(map(str, penalties)), "Status": "Locked"},
        {"Parameter": "DSNP_Event_Weight", "Value": 1, "Status": "Locked"},
        {"Parameter": "DSV_Event_Weight", "Value": 1, "Status": "Locked"},
    ]
    write_tsv(outroot / "parameter_manifest.tsv", parameter_rows, ["Parameter", "Value", "Status"])
    readme = """# Stage 09 hap38 load-only ideal parental haplotypes

This directory contains exact dynamic-programming IPHs with GWAS weight fixed at zero.
All38 is the primary panel and African35 is the sensitivity panel. The optimizer uses
frequency-inclusive Phoenix ALT-derived conserved SNPs, formal panel dSVs, 500-kb windows,
equal per-event dSNP/dSV weights and a primary breakpoint penalty P=15. The complete
P={0,5,15,40,100,300} sensitivity is retained.

The output is a theoretical load-minimizing donor mosaic, not an existing genome, a proven
superior parent, an F1 prediction or a reconstruction of the potato paper's original algorithm.
Conserved-region evidence is a proxy, dSNP numerically dominates the equal-event objective,
and P=15 is an inherited comparability setting rather than an empirical recombination rate.
"""
    write_text(outroot / "README.md", readme)
    if any(row["Status"] == "FAIL" for row in combined_checks):
        raise RuntimeError("Combined Stage 09 integrity failure")
    print(
        "[OK] Load-only IPH complete | "
        f"All38={all_summary['Residual_DSV']}/{all_summary['Residual_DSNP']} "
        f"breaks={all_summary['Breakpoint_Count']} donors={all_summary['Donor_Count']} | "
        f"African35={afr_summary['Residual_DSV']}/{afr_summary['Residual_DSNP']} "
        f"breaks={afr_summary['Breakpoint_Count']} donors={afr_summary['Donor_Count']} | GWAS=0"
    )


if __name__ == "__main__":
    main()
