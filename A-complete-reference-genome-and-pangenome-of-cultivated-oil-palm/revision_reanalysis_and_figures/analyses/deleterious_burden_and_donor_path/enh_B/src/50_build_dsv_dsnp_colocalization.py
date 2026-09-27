#!/usr/bin/env python3
"""Build panel-aware hap38 dSV-dSNP colocalization tables for Fig. e/f."""

import argparse
import bisect
import csv
import json
import math
import os
import re
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

try:
    from scipy import stats
except Exception:  # pragma: no cover - scipy is expected on the server
    stats = None


BASE = Path("${ANALYSIS_DIR}/21_MS/06_result/dSVs")
BIN_SIZE = 10_000
FLANK_BP = 1_000_000


def parse_args():
    parser = argparse.ArgumentParser(
        description="Build Fig. e/f dSV-dSNP correlation and flank enrichment tables."
    )
    parser.add_argument("--workdir", default=str(BASE))
    parser.add_argument("--dsv", default=str(BASE / "results/01_core_dsv/dsv_v5_candidates.tsv"))
    parser.add_argument(
        "--dsnp",
        default=str(BASE / "results/06_dSNP_phoenix_polarity_v2/polarity/dsnp_v2_phoenix_alt_derived_candidates.tsv"),
    )
    parser.add_argument("--fai", default=str(BASE / "input/Africa_hap2.fa.fai"))
    parser.add_argument("--sample-manifest", required=True)
    parser.add_argument("--sample-column", default="Sample_ID")
    parser.add_argument("--include-column", required=True)
    parser.add_argument("--carrier-column", default="Samples")
    parser.add_argument("--panel-label", required=True)
    parser.add_argument("--expected-samples", type=int, required=True)
    parser.add_argument("--expected-dsv", type=int, required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--bootstrap", type=int, default=1000)
    parser.add_argument("--background-replicates", type=int, default=20)
    parser.add_argument("--seed", type=int, default=20260630)
    parser.add_argument("--force", action="store_true")
    return parser.parse_args()


def chrom_sort_key(chrom):
    match = re.match(r"chr(\d+)([A-Za-z]*)$", chrom or "")
    if match:
        return (0, int(match.group(1)), match.group(2))
    return (1, chrom or "")


def parse_samples(value):
    if value is None:
        return []
    samples = []
    for part in re.split(r"[;,|]", str(value)):
        sample = part.strip()
        if sample and sample.upper() != "NA":
            samples.append(sample)
    return samples


def read_panel_samples(path, sample_column, include_column):
    fields, rows = read_rows(path)
    missing = [name for name in (sample_column, include_column) if name not in fields]
    if missing:
        raise ValueError(f"Sample manifest missing columns: {','.join(missing)}")
    samples = sorted(
        row[sample_column].strip()
        for row in rows
        if row.get(include_column, "").strip().lower() in {"yes", "true", "1"}
    )
    if not samples or len(samples) != len(set(samples)):
        raise ValueError("Selected sample IDs are empty or duplicated")
    return samples


def as_int(value, default=None):
    try:
        if value is None or value == "" or str(value).upper() == "NA":
            return default
        return int(float(str(value)))
    except ValueError:
        return default


def read_rows(path):
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = [dict(row) for row in reader]
        fieldnames = reader.fieldnames or []
    return fieldnames, rows


def write_tsv(path, rows, fieldnames, force=False):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists() and not force:
        raise FileExistsError(f"Refusing to overwrite existing file without --force: {path}")
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "NA") for field in fieldnames})


def write_text(path, text, force=False):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists() and not force:
        raise FileExistsError(f"Refusing to overwrite existing file without --force: {path}")
    with open(path, "w") as handle:
        handle.write(text)


def read_fai(path):
    chrom_lengths = {}
    with open(path) as handle:
        for line in handle:
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 2:
                chrom_lengths[fields[0]] = as_int(fields[1], 0)
    return chrom_lengths


def focal_position(row):
    start0 = as_int(row.get("Interval_Start0"))
    end0 = as_int(row.get("Interval_End0"))
    if start0 is not None and end0 is not None and end0 > start0:
        return int((start0 + end0) / 2)
    pos = as_int(row.get("Pos"))
    if pos is None:
        return None
    return pos


def pearson_summary(rows, label):
    x = np.array([float(row["dSV_Count"]) for row in rows], dtype=float)
    y = np.array([float(row["dSNP_Count"]) for row in rows], dtype=float)
    n = int(len(rows))
    if n < 3 or np.std(x) == 0 or np.std(y) == 0:
        r = math.nan
        p = math.nan
    elif stats is not None:
        res = stats.pearsonr(x, y)
        r = float(res.statistic)
        p = float(res.pvalue)
    else:
        r = float(np.corrcoef(x, y)[0, 1])
        t = abs(r) * math.sqrt((n - 2) / max(1e-15, 1 - r * r))
        p = math.erfc(t / math.sqrt(2))
    return {
        "Analysis": label,
        "Pearson_r": f"{r:.6g}",
        "P_value": f"{p:.6g}",
        "N": n,
        "dSV_Total": int(x.sum()),
        "dSNP_Total": int(y.sum()),
    }


def build_sample_chrom_counts(dsv_rows, dsnp_rows, chromosomes, samples, label):
    dsv_counts = Counter()
    dsnp_counts = Counter()
    for row in dsv_rows:
        chrom = row["Chrom"]
        if chrom not in chromosomes:
            continue
        for sample in parse_samples(row.get("Samples")):
            if sample in samples:
                dsv_counts[(sample, chrom)] += 1
    for row in dsnp_rows:
        chrom = row["Chrom"]
        if chrom not in chromosomes:
            continue
        for sample in parse_samples(row.get("Samples")):
            if sample in samples:
                dsnp_counts[(sample, chrom)] += 1
    out = []
    for sample in samples:
        for chrom in chromosomes:
            out.append(
                {
                    "Analysis": label,
                    "Sample": sample,
                    "Chrom": chrom,
                    "dSV_Count": dsv_counts[(sample, chrom)],
                    "dSNP_Count": dsnp_counts[(sample, chrom)],
                }
            )
    return out


def prepare_dsnp_index(dsnp_rows, chromosomes, main_samples):
    index = defaultdict(list)
    for row in dsnp_rows:
        chrom = row.get("Chrom")
        if chrom not in chromosomes:
            continue
        pos = as_int(row.get("Pos"))
        if pos is None:
            continue
        carriers = tuple(sorted(s for s in parse_samples(row.get("Samples")) if s in main_samples))
        if not carriers:
            continue
        index[chrom].append((pos, carriers, row.get("SNP_ID", "NA")))
    for chrom in index:
        index[chrom].sort(key=lambda item: item[0])
    return index


def functional_class(row):
    evidence = row.get("dSV_v5_Evidence_Class", "NA")
    if "CDS" in evidence:
        return "CDS_supported"
    if evidence == "Conserved_Region_Proxy":
        return "Conserved_region_proxy_only"
    gene_impact = row.get("Gene_Impact_Class", "NA")
    if gene_impact and gene_impact != "NA":
        return gene_impact
    return "Other"


def prepare_focal_infos(dsv_rows, chromosomes, main_samples):
    focal_infos = []
    for row in dsv_rows:
        chrom = row.get("Chrom")
        if chrom not in chromosomes:
            continue
        carriers = tuple(sorted(s for s in parse_samples(row.get("Samples")) if s in main_samples))
        if not carriers:
            continue
        center = focal_position(row)
        if center is None:
            continue
        focal_infos.append(
            {
                "Focal_DSV_ID": row.get("SV_Key") or row.get("SV_ID") or f"{chrom}:{center}",
                "Chrom": chrom,
                "Focal_Pos": center,
                "Focal_Carriers": carriers,
                "Focal_Carrier_Count": len(carriers),
                "Functional_Class": functional_class(row),
                "SVTYPE": row.get("SVTYPE", "NA") or "NA",
            }
        )
    return focal_infos


def count_carrier_bins(chrom_items, center, carriers, n_bins, positions=None):
    arr = np.zeros(n_bins, dtype=float)
    if not chrom_items or not carriers:
        return arr
    carrier_set = set(carriers)
    carrier_n = len(carrier_set)
    if positions is None:
        positions = [item[0] for item in chrom_items]
    left = center - FLANK_BP
    right = center + FLANK_BP
    lo = bisect.bisect_left(positions, left)
    hi = bisect.bisect_right(positions, right)
    for pos, dsnp_carriers, _snp_id in chrom_items[lo:hi]:
        rel = pos - center
        if rel < -FLANK_BP or rel > FLANK_BP:
            continue
        bin_idx = int((rel + FLANK_BP) // BIN_SIZE)
        if bin_idx >= n_bins:
            bin_idx = n_bins - 1
        hits = len(carrier_set & set(dsnp_carriers))
        if hits:
            arr[bin_idx] += hits / carrier_n
    return arr


def random_center(chrom_length, rng):
    if chrom_length <= 1:
        return 1
    return int(rng.integers(1, chrom_length + 1))


def build_window_counts(dsv_rows, dsnp_index, chromosomes, main_samples):
    n_bins = (2 * FLANK_BP) // BIN_SIZE
    focal_rows = []
    by_focal_counts = []
    for row in dsv_rows:
        chrom = row.get("Chrom")
        if chrom not in chromosomes:
            continue
        carriers = tuple(sorted(s for s in parse_samples(row.get("Samples")) if s in main_samples))
        if not carriers:
            continue
        center = focal_position(row)
        if center is None:
            continue
        chrom_items = dsnp_index.get(chrom, [])
        if not chrom_items:
            continue
        positions = [item[0] for item in chrom_items]
        left = center - FLANK_BP
        right = center + FLANK_BP
        lo = bisect.bisect_left(positions, left)
        hi = bisect.bisect_right(positions, right)
        coupling = np.zeros(n_bins, dtype=float)
        repulsion = np.zeros(n_bins, dtype=float)
        focal_set = set(carriers)
        focal_n = len(focal_set)
        noncarrier_n = len(main_samples) - focal_n
        if noncarrier_n <= 0:
            continue
        for pos, dsnp_carriers, _snp_id in chrom_items[lo:hi]:
            rel = pos - center
            if rel < -FLANK_BP or rel > FLANK_BP:
                continue
            bin_idx = int((rel + FLANK_BP) // BIN_SIZE)
            if bin_idx >= n_bins:
                bin_idx = n_bins - 1
            carrier_set = set(dsnp_carriers)
            coupling_hits = len(focal_set & carrier_set)
            repulsion_hits = len(carrier_set - focal_set)
            if coupling_hits:
                coupling[bin_idx] += coupling_hits / focal_n
            if repulsion_hits:
                repulsion[bin_idx] += repulsion_hits / noncarrier_n
        focal_id = row.get("SV_Key") or row.get("SV_ID") or f"{chrom}:{center}"
        for phase, arr in (("coupling", coupling), ("repulsion", repulsion)):
            for bin_idx, value in enumerate(arr):
                start = -FLANK_BP + bin_idx * BIN_SIZE
                end = start + BIN_SIZE
                focal_rows.append(
                    {
                        "Focal_DSV_ID": focal_id,
                        "Chrom": chrom,
                        "Focal_Pos": center,
                        "Focal_Carriers": ",".join(carriers),
                        "Focal_Carrier_Count": focal_n,
                        "Noncarrier_Sample_Count": noncarrier_n,
                        "Phase": phase,
                        "Distance_Start_bp": start,
                        "Distance_End_bp": end,
                        "Distance_Midpoint_kb": f"{(start + end) / 2000:.1f}",
                        "Accumulated_dSNP_Count": f"{value:.8g}",
                    }
                )
                by_focal_counts.append((phase, bin_idx, value))
    return focal_rows, by_focal_counts, n_bins


def summarize_windows(by_focal_counts, n_bins, bootstrap, seed):
    grouped = defaultdict(list)
    for phase, bin_idx, value in by_focal_counts:
        grouped[(phase, bin_idx)].append(float(value))
    rng = np.random.default_rng(seed)
    summary = []
    for phase in ("coupling", "repulsion"):
        for bin_idx in range(n_bins):
            values = np.asarray(grouped.get((phase, bin_idx), []), dtype=float)
            start = -FLANK_BP + bin_idx * BIN_SIZE
            end = start + BIN_SIZE
            if len(values) == 0:
                mean = ci_low = ci_high = math.nan
                n_eff = 0
            else:
                mean = float(values.mean())
                n_eff = int(len(values))
                if bootstrap > 0 and n_eff > 1:
                    sample_idx = rng.integers(0, n_eff, size=(bootstrap, n_eff))
                    boot_means = values[sample_idx].mean(axis=1)
                    ci_low, ci_high = np.quantile(boot_means, [0.025, 0.975])
                    ci_low = float(ci_low)
                    ci_high = float(ci_high)
                else:
                    ci_low = mean
                    ci_high = mean
            summary.append(
                {
                    "Phase": phase,
                    "Distance_Start_bp": start,
                    "Distance_End_bp": end,
                    "Distance_Start_kb": f"{start / 1000:.0f}",
                    "Distance_End_kb": f"{end / 1000:.0f}",
                    "Distance_Midpoint_kb": f"{(start + end) / 2000:.1f}",
                    "Mean_Accumulated_dSNP_Count": "NA" if math.isnan(mean) else f"{mean:.8g}",
                    "Bootstrap_CI_Low": "NA" if math.isnan(ci_low) else f"{ci_low:.8g}",
                    "Bootstrap_CI_High": "NA" if math.isnan(ci_high) else f"{ci_high:.8g}",
                    "Effective_Focal_DSV_Count": n_eff,
                    "Bootstrap_Replicates": bootstrap,
                }
            )
    return summary


def summarize_accumulated_counts(focal_rows, n_bins, bootstrap, seed, repulsion_mode):
    """Summarize total dSNP counts accumulated across focal dSVs per 10-kb bin."""
    grouped = defaultdict(list)
    for row in focal_rows:
        phase = row["Phase"]
        bin_idx = int((as_int(row["Distance_Start_bp"], 0) + FLANK_BP) // BIN_SIZE)
        value = float(row["Accumulated_dSNP_Count"])
        if phase == "repulsion" and repulsion_mode == "raw_noncarrier_total":
            value *= float(row["Noncarrier_Sample_Count"])
        grouped[(phase, bin_idx)].append(value)

    rng = np.random.default_rng(seed)
    summary = []
    for phase in ("coupling", "repulsion"):
        for bin_idx in range(n_bins):
            values = np.asarray(grouped.get((phase, bin_idx), []), dtype=float)
            start = -FLANK_BP + bin_idx * BIN_SIZE
            end = start + BIN_SIZE
            if len(values) == 0:
                total = mean = ci_low = ci_high = math.nan
                n_eff = 0
            else:
                total = float(values.sum())
                mean = float(values.mean())
                n_eff = int(len(values))
                if bootstrap > 0 and n_eff > 1:
                    sample_idx = rng.integers(0, n_eff, size=(bootstrap, n_eff))
                    boot_totals = values[sample_idx].sum(axis=1)
                    ci_low, ci_high = np.quantile(boot_totals, [0.025, 0.975])
                    ci_low = float(ci_low)
                    ci_high = float(ci_high)
                else:
                    ci_low = total
                    ci_high = total
            summary.append(
                {
                    "Phase": phase,
                    "Repulsion_Mode": repulsion_mode,
                    "Distance_Start_bp": start,
                    "Distance_End_bp": end,
                    "Distance_Start_kb": f"{start / 1000:.0f}",
                    "Distance_End_kb": f"{end / 1000:.0f}",
                    "Distance_Midpoint_kb": f"{(start + end) / 2000:.1f}",
                    "Accumulated_dSNP_Count": "NA" if math.isnan(total) else f"{total:.8g}",
                    "Mean_Per_Focal_DSV": "NA" if math.isnan(mean) else f"{mean:.8g}",
                    "Bootstrap_CI_Low": "NA" if math.isnan(ci_low) else f"{ci_low:.8g}",
                    "Bootstrap_CI_High": "NA" if math.isnan(ci_high) else f"{ci_high:.8g}",
                    "Effective_Focal_DSV_Count": n_eff,
                    "Bootstrap_Replicates": bootstrap,
                }
            )
    return summary


def summarize_series_counts(series_counts, n_bins, bootstrap, seed):
    grouped = defaultdict(list)
    for series, bin_idx, value in series_counts:
        grouped[(series, bin_idx)].append(float(value))
    rng = np.random.default_rng(seed)
    series_order = sorted({series for series, _bin_idx, _value in series_counts})
    summary = []
    for series in series_order:
        for bin_idx in range(n_bins):
            values = np.asarray(grouped.get((series, bin_idx), []), dtype=float)
            start = -FLANK_BP + bin_idx * BIN_SIZE
            end = start + BIN_SIZE
            if len(values) == 0:
                total = mean = ci_low = ci_high = mean_ci_low = mean_ci_high = math.nan
                n_eff = 0
            else:
                total = float(values.sum())
                mean = float(values.mean())
                n_eff = int(len(values))
                if bootstrap > 0 and n_eff > 1:
                    sample_idx = rng.integers(0, n_eff, size=(bootstrap, n_eff))
                    boot_totals = values[sample_idx].sum(axis=1)
                    ci_low, ci_high = np.quantile(boot_totals, [0.025, 0.975])
                    mean_ci_low = float(ci_low) / n_eff
                    mean_ci_high = float(ci_high) / n_eff
                    ci_low = float(ci_low)
                    ci_high = float(ci_high)
                else:
                    ci_low = ci_high = total
                    mean_ci_low = mean_ci_high = mean
            summary.append(
                {
                    "Series": series,
                    "Distance_Start_bp": start,
                    "Distance_End_bp": end,
                    "Distance_Start_kb": f"{start / 1000:.0f}",
                    "Distance_End_kb": f"{end / 1000:.0f}",
                    "Distance_Midpoint_kb": f"{(start + end) / 2000:.1f}",
                    "Accumulated_dSNP_Count": "NA" if math.isnan(total) else f"{total:.8g}",
                    "Mean_Per_Focal_DSV": "NA" if math.isnan(mean) else f"{mean:.8g}",
                    "Bootstrap_CI_Low": "NA" if math.isnan(ci_low) else f"{ci_low:.8g}",
                    "Bootstrap_CI_High": "NA" if math.isnan(ci_high) else f"{ci_high:.8g}",
                    "Mean_CI_Low": "NA" if math.isnan(mean_ci_low) else f"{mean_ci_low:.8g}",
                    "Mean_CI_High": "NA" if math.isnan(mean_ci_high) else f"{mean_ci_high:.8g}",
                    "Effective_Focal_DSV_Count": n_eff,
                    "Bootstrap_Replicates": bootstrap,
                }
            )
    return summary


def build_observed_vs_background(focal_infos, dsnp_index, chrom_lengths, background_replicates, bootstrap, seed):
    n_bins = (2 * FLANK_BP) // BIN_SIZE
    rng = np.random.default_rng(seed)
    positions_by_chrom = {chrom: [item[0] for item in items] for chrom, items in dsnp_index.items()}
    detail_rows = []
    series_counts = []
    for focal in focal_infos:
        chrom = focal["Chrom"]
        chrom_items = dsnp_index.get(chrom, [])
        positions = positions_by_chrom.get(chrom, [])
        carriers = focal["Focal_Carriers"]
        observed = count_carrier_bins(chrom_items, focal["Focal_Pos"], carriers, n_bins, positions)
        background = np.zeros(n_bins, dtype=float)
        reps = max(1, int(background_replicates))
        for _rep in range(reps):
            center = random_center(chrom_lengths.get(chrom, 0), rng)
            background += count_carrier_bins(chrom_items, center, carriers, n_bins, positions)
        background = background / reps
        for series, arr in (("Observed_dSV_carrier", observed), ("Matched_random_background", background)):
            for bin_idx, value in enumerate(arr):
                start = -FLANK_BP + bin_idx * BIN_SIZE
                end = start + BIN_SIZE
                detail_rows.append(
                    {
                        "Focal_DSV_ID": focal["Focal_DSV_ID"],
                        "Chrom": chrom,
                        "Focal_Pos": focal["Focal_Pos"],
                        "Focal_Carriers": ",".join(carriers),
                        "Series": series,
                        "Background_Replicates": reps if series == "Matched_random_background" else 0,
                        "Distance_Start_bp": start,
                        "Distance_End_bp": end,
                        "Distance_Midpoint_kb": f"{(start + end) / 2000:.1f}",
                        "Accumulated_dSNP_Count": f"{value:.8g}",
                    }
                )
                series_counts.append((series, bin_idx, value))
    summary = summarize_series_counts(series_counts, n_bins, bootstrap, seed + 101)
    return detail_rows, summary


def build_functional_class_counts(focal_infos, dsnp_index, bootstrap, seed):
    n_bins = (2 * FLANK_BP) // BIN_SIZE
    positions_by_chrom = {chrom: [item[0] for item in items] for chrom, items in dsnp_index.items()}
    detail_rows = []
    series_counts = []
    for focal in focal_infos:
        chrom = focal["Chrom"]
        chrom_items = dsnp_index.get(chrom, [])
        positions = positions_by_chrom.get(chrom, [])
        carriers = focal["Focal_Carriers"]
        arr = count_carrier_bins(chrom_items, focal["Focal_Pos"], carriers, n_bins, positions)
        series = focal["Functional_Class"]
        for bin_idx, value in enumerate(arr):
            start = -FLANK_BP + bin_idx * BIN_SIZE
            end = start + BIN_SIZE
            detail_rows.append(
                {
                    "Focal_DSV_ID": focal["Focal_DSV_ID"],
                    "Chrom": chrom,
                    "Focal_Pos": focal["Focal_Pos"],
                    "Focal_Carriers": ",".join(carriers),
                    "Functional_Class": series,
                    "Distance_Start_bp": start,
                    "Distance_End_bp": end,
                    "Distance_Midpoint_kb": f"{(start + end) / 2000:.1f}",
                    "Accumulated_dSNP_Count": f"{value:.8g}",
                }
            )
            series_counts.append((series, bin_idx, value))
    summary = summarize_series_counts(series_counts, n_bins, bootstrap, seed + 211)
    return detail_rows, summary


def build_svtype_class_counts(focal_infos, dsnp_index, bootstrap, seed, keep_types=("DEL", "INS")):
    n_bins = (2 * FLANK_BP) // BIN_SIZE
    positions_by_chrom = {chrom: [item[0] for item in items] for chrom, items in dsnp_index.items()}
    detail_rows = []
    series_counts = []
    keep_types = set(keep_types)
    for focal in focal_infos:
        series = focal.get("SVTYPE", "NA")
        if series not in keep_types:
            continue
        chrom = focal["Chrom"]
        chrom_items = dsnp_index.get(chrom, [])
        positions = positions_by_chrom.get(chrom, [])
        carriers = focal["Focal_Carriers"]
        arr = count_carrier_bins(chrom_items, focal["Focal_Pos"], carriers, n_bins, positions)
        for bin_idx, value in enumerate(arr):
            start = -FLANK_BP + bin_idx * BIN_SIZE
            end = start + BIN_SIZE
            detail_rows.append(
                {
                    "Focal_DSV_ID": focal["Focal_DSV_ID"],
                    "Chrom": chrom,
                    "Focal_Pos": focal["Focal_Pos"],
                    "Focal_Carriers": ",".join(carriers),
                    "SVTYPE": series,
                    "Distance_Start_bp": start,
                    "Distance_End_bp": end,
                    "Distance_Midpoint_kb": f"{(start + end) / 2000:.1f}",
                    "Accumulated_dSNP_Count": f"{value:.8g}",
                }
            )
            series_counts.append((series, bin_idx, value))
    summary = summarize_series_counts(series_counts, n_bins, bootstrap, seed + 313)
    return detail_rows, summary


def integrity_rows(
    dsv_rows,
    dsnp_rows,
    chromosomes,
    panel_samples,
    expected_samples,
    expected_dsv,
    counts,
    corr,
    window_summary,
    focal_infos,
):
    def status(ok):
        return "PASS" if ok else "FAIL"

    expected_grid_rows = expected_samples * len(chromosomes)
    phases = {row["Phase"] for row in window_summary}
    ci_ok = all(
        row["Bootstrap_CI_Low"] != "NA"
        and row["Bootstrap_CI_High"] != "NA"
        and float(row["Bootstrap_CI_Low"]) >= 0
        and float(row["Bootstrap_CI_High"]) >= 0
        for row in window_summary
    )
    return [
        {"Check": "dSV_input_rows", "Status": status(len(dsv_rows) == expected_dsv), "Value": len(dsv_rows)},
        {"Check": "dSNP_input_rows", "Status": status(len(dsnp_rows) > 0), "Value": len(dsnp_rows)},
        {"Check": "chromosome_count", "Status": status(len(chromosomes) == 16), "Value": len(chromosomes)},
        {
            "Check": "panel_sample_count",
            "Status": status(len(panel_samples) == expected_samples),
            "Value": len(panel_samples),
        },
        {
            "Check": "focal_dSV_count",
            "Status": status(len(focal_infos) == expected_dsv),
            "Value": len(focal_infos),
        },
        {
            "Check": "focal_dSV_singleton_carriers",
            "Status": status(all(row["Focal_Carrier_Count"] == 1 for row in focal_infos)),
            "Value": all(row["Focal_Carrier_Count"] == 1 for row in focal_infos),
        },
        {
            "Check": "fig_e_panel_grid_rows",
            "Status": status(len(counts) == expected_grid_rows),
            "Value": len(counts),
        },
        {
            "Check": "fig_e_panel_pearson_N",
            "Status": status(corr["N"] == expected_grid_rows),
            "Value": corr["N"],
        },
        {"Check": "fig_f_bin_size_bp", "Status": status(BIN_SIZE == 10000), "Value": BIN_SIZE},
        {
            "Check": "fig_f_distance_range_kb",
            "Status": status(window_summary[0]["Distance_Start_kb"] == "-1000" and window_summary[-1]["Distance_End_kb"] == "1000"),
            "Value": f"{window_summary[0]['Distance_Start_kb']}..{window_summary[-1]['Distance_End_kb']}",
        },
        {"Check": "fig_f_phase_presence", "Status": status(phases == {"coupling", "repulsion"}), "Value": ",".join(sorted(phases))},
        {"Check": "fig_f_ci_nonnegative_non_na", "Status": status(ci_ok), "Value": ci_ok},
    ]


def main():
    args = parse_args()
    outdir = Path(args.outdir)
    if outdir.exists() and any(outdir.rglob("*")) and not args.force:
        raise SystemExit(f"[ERROR] {outdir} exists and is not empty; use --force only for an intentional rerun.")
    for subdir in ("figures", "notes", "qa"):
        (outdir / subdir).mkdir(parents=True, exist_ok=True)

    dsv_fields, dsv_rows = read_rows(args.dsv)
    dsnp_fields, dsnp_rows = read_rows(args.dsnp)
    if args.carrier_column not in dsv_fields:
        raise SystemExit(f"[ERROR] dSV table lacks carrier column: {args.carrier_column}")
    if "Samples" not in dsnp_fields:
        raise SystemExit("[ERROR] dSNP table lacks Samples column")
    if args.carrier_column != "Samples":
        for row in dsv_rows:
            row["Samples"] = row.get(args.carrier_column, "")

    panel_samples = read_panel_samples(args.sample_manifest, args.sample_column, args.include_column)
    if len(panel_samples) != args.expected_samples:
        raise SystemExit(
            f"[ERROR] Panel sample count {len(panel_samples)} != expected {args.expected_samples}"
        )
    panel_sample_set = set(panel_samples)
    chrom_lengths = read_fai(args.fai)
    chromosomes = sorted([chrom for chrom in chrom_lengths if re.match(r"chr\d+B$", chrom)], key=chrom_sort_key)
    dsv_chroms = {row.get("Chrom") for row in dsv_rows}
    dsnp_chroms = {row.get("Chrom") for row in dsnp_rows}
    chromosomes = [chrom for chrom in chromosomes if chrom in dsv_chroms or chrom in dsnp_chroms]

    panel_counts = build_sample_chrom_counts(
        dsv_rows, dsnp_rows, chromosomes, panel_samples, args.panel_label
    )
    panel_corr = pearson_summary(panel_counts, args.panel_label)

    dsnp_index = prepare_dsnp_index(dsnp_rows, chromosomes, panel_sample_set)
    focal_infos = prepare_focal_infos(dsv_rows, set(chromosomes), panel_sample_set)
    focal_rows, by_focal_counts, n_bins = build_window_counts(
        dsv_rows, dsnp_index, set(chromosomes), panel_sample_set
    )
    window_summary = summarize_windows(by_focal_counts, n_bins, args.bootstrap, args.seed)
    accumulated_summary = summarize_accumulated_counts(
        focal_rows, n_bins, args.bootstrap, args.seed + 11, "single_noncarrier_equivalent"
    )
    accumulated_raw_summary = summarize_accumulated_counts(
        focal_rows, n_bins, args.bootstrap, args.seed + 17, "raw_noncarrier_total"
    )
    obs_bg_detail, obs_bg_summary = build_observed_vs_background(
        focal_infos, dsnp_index, chrom_lengths, args.background_replicates, args.bootstrap, args.seed + 31
    )
    functional_detail, functional_summary = build_functional_class_counts(
        focal_infos, dsnp_index, args.bootstrap, args.seed + 41
    )
    svtype_detail, svtype_summary = build_svtype_class_counts(
        focal_infos, dsnp_index, args.bootstrap, args.seed + 53
    )

    input_summary = [
        {"Metric": "dSV_Input", "Value": args.dsv},
        {"Metric": "dSV_Rows", "Value": len(dsv_rows)},
        {"Metric": "dSNP_Input", "Value": args.dsnp},
        {"Metric": "dSNP_Rows", "Value": len(dsnp_rows)},
        {"Metric": "FAI", "Value": args.fai},
        {"Metric": "Sample_Manifest", "Value": args.sample_manifest},
        {"Metric": "Include_Column", "Value": args.include_column},
        {"Metric": "dSV_Carrier_Column", "Value": args.carrier_column},
        {"Metric": "Panel", "Value": args.panel_label},
        {"Metric": "Chromosomes", "Value": ",".join(chromosomes)},
        {"Metric": "Panel_Sample_Count", "Value": len(panel_samples)},
        {"Metric": "Panel_Samples", "Value": ",".join(panel_samples)},
        {"Metric": "Focal_dSV_Count", "Value": len(focal_infos)},
        {"Metric": "Bootstrap_Replicates", "Value": args.bootstrap},
        {"Metric": "Matched_Background_Replicates_Per_Focal_DSV", "Value": args.background_replicates},
        {"Metric": "Random_Seed", "Value": args.seed},
    ]

    write_tsv(
        outdir / "fig_e_sample_chrom_counts.tsv",
        panel_counts,
        ["Analysis", "Sample", "Chrom", "dSV_Count", "dSNP_Count"],
        args.force,
    )
    write_tsv(
        outdir / "fig_e_correlation_summary.tsv",
        [panel_corr],
        ["Analysis", "Pearson_r", "P_value", "N", "dSV_Total", "dSNP_Total"],
        args.force,
    )
    write_tsv(
        outdir / "fig_f_window_counts_by_focal_dsv.tsv",
        focal_rows,
        [
            "Focal_DSV_ID",
            "Chrom",
            "Focal_Pos",
            "Focal_Carriers",
            "Focal_Carrier_Count",
            "Noncarrier_Sample_Count",
            "Phase",
            "Distance_Start_bp",
            "Distance_End_bp",
            "Distance_Midpoint_kb",
            "Accumulated_dSNP_Count",
        ],
        args.force,
    )
    write_tsv(
        outdir / "fig_f_window_summary.tsv",
        window_summary,
        [
            "Phase",
            "Distance_Start_bp",
            "Distance_End_bp",
            "Distance_Start_kb",
            "Distance_End_kb",
            "Distance_Midpoint_kb",
            "Mean_Accumulated_dSNP_Count",
            "Bootstrap_CI_Low",
            "Bootstrap_CI_High",
            "Effective_Focal_DSV_Count",
            "Bootstrap_Replicates",
        ],
        args.force,
    )
    accumulated_fields = [
        "Phase",
        "Repulsion_Mode",
        "Distance_Start_bp",
        "Distance_End_bp",
        "Distance_Start_kb",
        "Distance_End_kb",
        "Distance_Midpoint_kb",
        "Accumulated_dSNP_Count",
        "Mean_Per_Focal_DSV",
        "Bootstrap_CI_Low",
        "Bootstrap_CI_High",
        "Effective_Focal_DSV_Count",
        "Bootstrap_Replicates",
    ]
    write_tsv(
        outdir / "fig_f_accumulated_summary.tsv",
        accumulated_summary,
        accumulated_fields,
        args.force,
    )
    write_tsv(
        outdir / "fig_f_accumulated_raw_noncarrier_summary.tsv",
        accumulated_raw_summary,
        accumulated_fields,
        args.force,
    )
    series_summary_fields = [
        "Series",
        "Distance_Start_bp",
        "Distance_End_bp",
        "Distance_Start_kb",
        "Distance_End_kb",
        "Distance_Midpoint_kb",
        "Accumulated_dSNP_Count",
        "Mean_Per_Focal_DSV",
        "Bootstrap_CI_Low",
        "Bootstrap_CI_High",
        "Mean_CI_Low",
        "Mean_CI_High",
        "Effective_Focal_DSV_Count",
        "Bootstrap_Replicates",
    ]
    write_tsv(
        outdir / "fig_f_observed_vs_matched_background_by_focal_dsv.tsv",
        obs_bg_detail,
        [
            "Focal_DSV_ID",
            "Chrom",
            "Focal_Pos",
            "Focal_Carriers",
            "Series",
            "Background_Replicates",
            "Distance_Start_bp",
            "Distance_End_bp",
            "Distance_Midpoint_kb",
            "Accumulated_dSNP_Count",
        ],
        args.force,
    )
    write_tsv(
        outdir / "fig_f_observed_vs_matched_background_summary.tsv",
        obs_bg_summary,
        series_summary_fields,
        args.force,
    )
    write_tsv(
        outdir / "fig_f_functional_class_by_focal_dsv.tsv",
        functional_detail,
        [
            "Focal_DSV_ID",
            "Chrom",
            "Focal_Pos",
            "Focal_Carriers",
            "Functional_Class",
            "Distance_Start_bp",
            "Distance_End_bp",
            "Distance_Midpoint_kb",
            "Accumulated_dSNP_Count",
        ],
        args.force,
    )
    write_tsv(
        outdir / "fig_f_functional_class_summary.tsv",
        functional_summary,
        series_summary_fields,
        args.force,
    )
    write_tsv(
        outdir / "fig_f_svtype_class_by_focal_dsv.tsv",
        svtype_detail,
        [
            "Focal_DSV_ID",
            "Chrom",
            "Focal_Pos",
            "Focal_Carriers",
            "SVTYPE",
            "Distance_Start_bp",
            "Distance_End_bp",
            "Distance_Midpoint_kb",
            "Accumulated_dSNP_Count",
        ],
        args.force,
    )
    write_tsv(
        outdir / "fig_f_svtype_class_summary.tsv",
        svtype_summary,
        series_summary_fields,
        args.force,
    )
    write_tsv(outdir / "qa/input_summary.tsv", input_summary, ["Metric", "Value"], args.force)
    integrity = integrity_rows(
        dsv_rows,
        dsnp_rows,
        chromosomes,
        panel_samples,
        args.expected_samples,
        args.expected_dsv,
        panel_counts,
        panel_corr,
        window_summary,
        focal_infos,
    )
    write_tsv(outdir / "qa/final_integrity_summary.tsv", integrity, ["Check", "Status", "Value"], args.force)

    type_counts = Counter(row["SVTYPE"] for row in focal_infos)
    metadata = {
        "Panel": args.panel_label,
        "Sample_Count": len(panel_samples),
        "Focal_dSV_Count": len(focal_infos),
        "dSNP_Input_Count": len(dsnp_rows),
        "Chromosomes": chromosomes,
        "Figure_Families": {
            "Fig_e": "scatter_regression",
            "Fig_f": "ordered_distance_profile",
            "Combined": "multi_panel_manuscript_layout",
        },
        "Pattern_Documents": [
            "scatter-regression-marginal.md",
            "multi-panel-manuscript-layout.md",
        ],
        "Bootstrap_Replicates": args.bootstrap,
        "Background_Replicates_Per_Focal_dSV": args.background_replicates,
        "Seed": args.seed,
        "Reference": "Africa_hap2",
        "Coordinate_System": "1-based SNP positions; dSV midpoint from 0-based half-open interval",
        "Limitations": [
            "dSNP catalogue retains the locked legacy reference-match selection rule",
            "conserved-region evidence is a proxy and not experimental constraint evidence",
        ],
    }
    write_text(
        outdir / "notes/analysis_metadata.json",
        json.dumps(metadata, indent=2, ensure_ascii=False) + "\n",
        args.force,
    )

    legend = f"""# Figure e/f legend draft

Panel scope: **{args.panel_label}**.

**e,** Pearson correlation between Phoenix-polarized derived SNP (dSNP) and derived structural variant (dSV) counts across sample-by-chromosome bins. Each point represents one sample on one chromosome ({len(panel_samples)} samples x 16 chromosomes). The Pearson correlation is r = {float(panel_corr['Pearson_r']):.2f}, P = {float(panel_corr['P_value']):.2g} (N = {panel_corr['N']}).

**f,** Accumulated dSNP counts around {len(focal_infos)} focal dSVs within +/-1 Mb, summarized in 10-kb non-overlapping windows. The main comparison is the observed focal-dSV carrier background versus 20 random centers matched by carrier sample and chromosome. Curves are split at 0 kb so bins on opposite sides of the focal dSV are not connected. Bootstrap intervals in the source tables are based on {args.bootstrap} focal-dSV resamples.

Supplementary stratifications compare functional evidence classes and DEL versus INS. DEL n = {type_counts.get('DEL', 0)}; INS n = {type_counts.get('INS', 0)}; INV and DUP are not shown in the SVTYPE curve because their focal counts are too small for a stable two-class comparison.

Interpretation limits: the dSNP catalogue is a legacy-compatible proxy set because the locked reference-match selection rule excludes reference-mismatch rows. Conserved-region support is a multi-outgroup conserved/alignable-region proxy, not GERP, phyloP, phastCons, or experimental validation.
"""
    write_text(outdir / "notes/figure_ef_legend.md", legend, args.force)

    if any(row["Status"] != "PASS" for row in integrity):
        raise SystemExit("[ERROR] One or more final integrity checks failed.")
    print(f"[OK] Wrote Fig. e/f analysis tables to {outdir}")


if __name__ == "__main__":
    main()
