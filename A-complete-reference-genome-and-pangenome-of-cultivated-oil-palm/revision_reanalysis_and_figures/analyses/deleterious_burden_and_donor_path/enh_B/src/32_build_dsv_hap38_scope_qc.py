#!/usr/bin/env python3
"""All38/African35 scope sensitivity and old-versus-hap38 dSV comparison."""
import argparse
import bisect
import csv
import json
import math
import re
from collections import Counter, defaultdict
from pathlib import Path

CONSIDERED_SVTYPES = {"DEL", "INS", "DUP", "INV"}
RARE_THRESHOLD = 0.05
EXPECTED_ALL38 = 38
EXPECTED_AFRICAN35 = 35
EXPECTED_NON_AFRICAN = 3
EXPECTED_FORMAL_DSV = 1480


def parse_args():
    base = Path("${ANALYSIS_DIR}/21_MS/06_result/dSVs")
    root = base / "results-8.9"
    p = argparse.ArgumentParser(description="Build All38/African35 dSV scope sensitivity analysis.")
    p.add_argument("--catalog", default=str(root / "01_core_dsv/sv_catalog.dsv_hap38.tsv"))
    p.add_argument("--candidates", default=str(root / "01_core_dsv/dsv_hap38_candidates.tsv"))
    p.add_argument("--manifest", default=str(root / "config/Sample_Manifest.tsv"))
    p.add_argument("--sample-burden", default=str(root / "02_summary/dsv_sample_burden.tsv"))
    p.add_argument("--old-candidates", default=str(base / "results/01_core_dsv/dsv_v5_candidates.tsv"))
    p.add_argument("--outdir", default=str(root / "04_species_qc"))
    return p.parse_args()


def read_tsv(path):
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return list(reader), list(reader.fieldnames or [])


def write_tsv(path, rows, fields):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        raise FileExistsError(f"Refusing to overwrite: {path}")
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=list(fields), extrasaction="ignore", lineterminator="\n")
        writer.writeheader()
        writer.writerows({field: row.get(field, "NA") for field in fields} for row in rows)


def split_samples(value):
    text = str(value or "").strip()
    if not text or text.upper() == "NA":
        return []
    return [x.strip() for x in re.split(r"[;,|]+", text) if x.strip() and x.strip().upper() != "NA"]


def as_int(value, default=0):
    try:
        return int(float(str(value)))
    except (TypeError, ValueError):
        return default


def midpoint(row):
    pos = as_int(row.get("Pos"), -1)
    if pos >= 0:
        return pos
    return (as_int(row.get("Start")) + as_int(row.get("End"))) // 2


def functional_supported(row):
    return (
        as_int(row.get("CDS_Overlap_Count")) > 0
        or row.get("Conserved_Region_Overlap") == "Yes"
        or as_int(row.get("Conserved_Region_Covered_Bases")) > 0
    )


def scope_of_carriers(carriers, african, non_african):
    carrier_set = set(carriers)
    if carrier_set and carrier_set <= african:
        return "African35"
    if carrier_set and carrier_set <= non_african:
        return "American_or_Hybrid3"
    if carrier_set & african and carrier_set & non_african:
        return "Cross_Scope"
    return "Unknown"


def counter_rows(counter, key_name, total=None, extra=None):
    rows = []
    for key in sorted(counter):
        row = {key_name: key, "Count": counter[key]}
        if total is not None:
            row["Fraction"] = f"{counter[key] / total:.8f}" if total else "NA"
        if extra:
            row.update(extra)
        rows.append(row)
    return rows


def nearest_matches(old_rows, new_rows, distance=500):
    index = defaultdict(list)
    for row in new_rows:
        index[(row.get("Chrom"), row.get("SVTYPE"))].append((midpoint(row), row.get("SV_Key"), row.get("SV_ID")))
    positions = {}
    for key in index:
        index[key].sort()
        positions[key] = [x[0] for x in index[key]]
    detail = []
    matched = 0
    for row in old_rows:
        key = (row.get("Chrom"), row.get("SVTYPE"))
        op = midpoint(row)
        arr = index.get(key, [])
        pos = positions.get(key, [])
        i = bisect.bisect_left(pos, op)
        candidates = []
        for j in (i - 1, i):
            if 0 <= j < len(arr):
                candidates.append(arr[j])
        best = min(candidates, key=lambda x: abs(x[0] - op)) if candidates else None
        is_match = best is not None and abs(best[0] - op) <= distance
        matched += int(is_match)
        detail.append({
            "Old_SV_Key": row.get("SV_Key", "NA"),
            "Old_SV_ID": row.get("SV_ID", "NA"),
            "Chrom": row.get("Chrom", "NA"),
            "SVTYPE": row.get("SVTYPE", "NA"),
            "Old_Position_bp": op,
            "Nearest_New_SV_Key": best[1] if best else "NA",
            "Nearest_New_SV_ID": best[2] if best else "NA",
            "Nearest_New_Position_bp": best[0] if best else "NA",
            "Distance_bp": abs(best[0] - op) if best else "NA",
            "Matched_Same_Chrom_Type_Within_500bp": "Yes" if is_match else "No",
        })
    return detail, matched


def main():
    args = parse_args()
    outdir = Path(args.outdir)
    if outdir.exists() and any(path.is_file() for path in outdir.rglob("*")):
        raise FileExistsError(f"Refusing to write into non-empty output directory: {outdir}")
    outdir.mkdir(parents=True, exist_ok=True)
    (outdir / "figures").mkdir(exist_ok=True)
    (outdir / "notes").mkdir(exist_ok=True)
    (outdir / "qa").mkdir(exist_ok=True)

    catalog, catalog_fields = read_tsv(args.catalog)
    formal, formal_fields = read_tsv(args.candidates)
    manifest, manifest_fields = read_tsv(args.manifest)
    burden, _ = read_tsv(args.sample_burden)
    old, _ = read_tsv(args.old_candidates)

    all38 = {r["Sample_ID"] for r in manifest if r.get("Include_All38") == "Yes"}
    african = {r["Sample_ID"] for r in manifest if r.get("Include_African35") == "Yes"}
    non_african = all38 - african
    group_by_sample = {r["Sample_ID"]: r.get("Population_Scope", "NA") for r in manifest}
    burden_by_sample = {r["Sample_Name"]: r for r in burden}
    formal_keys = {r.get("SV_Key") for r in formal}

    scope_counts = Counter()
    scope_type_counts = Counter()
    scope_layer_counts = Counter()
    scope_function_counts = Counter()
    sample_counts = Counter()
    non_african_detail = []
    african_fixed_subset = []
    unknown_carriers = set()
    for row in formal:
        carriers = split_samples(row.get("Samples"))
        unknown_carriers.update(set(carriers) - all38)
        scope = scope_of_carriers(carriers, african, non_african)
        scope_counts[scope] += 1
        scope_type_counts[(scope, row.get("SVTYPE", "NA"))] += 1
        scope_layer_counts[(scope, row.get("Evidence_Layer", "NA"))] += 1
        scope_function_counts[(scope, row.get("dSV_v5_Evidence_Class", "NA"))] += 1
        for sample in carriers:
            sample_counts[sample] += 1
        if carriers and set(carriers) <= african:
            african_fixed_subset.append(row)
        if set(carriers) & non_african:
            out = dict(row)
            out["Carrier_Scope"] = scope
            non_african_detail.append(out)

    sample_rows = []
    for row in manifest:
        sample = row["Sample_ID"]
        b = burden_by_sample.get(sample, {})
        old_comparator = row.get("Old_Collapsed_ID")
        if not old_comparator or old_comparator == "NA":
            old_comparator = sample if row.get("New8_Mode") == "Reused30" else "NA"
        sample_rows.append({
            "Sample_ID": sample,
            "Population_Scope": row.get("Population_Scope", "NA"),
            "Include_All38": row.get("Include_All38", "NA"),
            "Include_African35": row.get("Include_African35", "NA"),
            "Formal_All38_dSV_Count": sample_counts.get(sample, 0),
            "Formal_All38_dSV_Total_Length_bp": b.get("DSV_Total_Length_bp", 0),
            "Old_Comparator": old_comparator,
            "Sample_Update_Mode": row.get("New8_Mode", "NA"),
        })

    scope_summary = []
    for scope in ["African35", "American_or_Hybrid3", "Cross_Scope", "Unknown"]:
        scope_summary.append({
            "Carrier_Scope": scope,
            "Formal_All38_dSV_Count": scope_counts.get(scope, 0),
            "Fraction_Of_All38_Formal": f"{scope_counts.get(scope, 0) / len(formal):.8f}" if formal else "NA",
        })

    by_type = [
        {"Carrier_Scope": scope, "SVTYPE": svtype, "Formal_All38_dSV_Count": count}
        for (scope, svtype), count in sorted(scope_type_counts.items())
    ]
    by_layer = [
        {"Carrier_Scope": scope, "Evidence_Layer": layer, "Formal_All38_dSV_Count": count}
        for (scope, layer), count in sorted(scope_layer_counts.items())
    ]
    by_function = [
        {"Carrier_Scope": scope, "Functional_Evidence_Class": evidence, "Formal_All38_dSV_Count": count}
        for (scope, evidence), count in sorted(scope_function_counts.items())
    ]

    # Re-identify with the same formal v5 rule after restricting carriers to African35.
    african_reidentified = []
    african_rule_violations = 0
    for row in catalog:
        carriers = [s for s in split_samples(row.get("Samples")) if s in african]
        if not carriers:
            continue
        count = len(carriers)
        freq = count / EXPECTED_AFRICAN35
        keep = (
            row.get("SVTYPE") in CONSIDERED_SVTYPES
            and freq <= RARE_THRESHOLD
            and row.get("ALT_Polarity_v5") == "ALT_Derived"
            and functional_supported(row)
        )
        if not keep:
            continue
        out = dict(row)
        out.update({
            "Samples_African35": ";".join(carriers),
            "Sample_Count_African35": count,
            "Frequency_African35": f"{freq:.8f}",
            "African35_dSV_Flag": "Yes",
            "African35_Change_Class": "Retained_from_All38_Formal" if row.get("SV_Key") in formal_keys else "Newly_Qualifying_African35",
        })
        african_reidentified.append(out)
        if count > 1 or freq > RARE_THRESHOLD or not functional_supported(row):
            african_rule_violations += 1

    reid_keys = {r.get("SV_Key") for r in african_reidentified}
    retained_keys = reid_keys & formal_keys
    newly_keys = reid_keys - formal_keys
    fixed_keys = {r.get("SV_Key") for r in african_fixed_subset}
    reid_summary = [
        {"Metric": "All38_Formal_dSV", "Value": len(formal), "Definition": "Formal main result"},
        {"Metric": "All38_Formal_With_African35_Carrier", "Value": len(african_fixed_subset), "Definition": "Fixed formal set; carriers restricted to African35"},
        {"Metric": "All38_Formal_With_American_or_Hybrid3_Carrier", "Value": len(non_african_detail), "Definition": "Fixed formal set; at least one carrier in non-African/hybrid3"},
        {"Metric": "African35_Reidentified_dSV", "Value": len(african_reidentified), "Definition": "Recomputed carrier count and 5% threshold within African35"},
        {"Metric": "African35_Retained_From_All38_Formal", "Value": len(retained_keys), "Definition": "Intersection of African35 reidentified and All38 formal"},
        {"Metric": "African35_Newly_Qualifying", "Value": len(newly_keys), "Definition": "Not formal in All38 but qualifies after scope restriction"},
        {"Metric": "All38_Formal_Not_In_African35_Reidentified", "Value": len(formal_keys - reid_keys), "Definition": "Primarily non-African/hybrid-carried formal dSVs"},
    ]
    reid_type = Counter((r["African35_Change_Class"], r.get("SVTYPE", "NA")) for r in african_reidentified)
    reid_type_rows = [
        {"African35_Change_Class": change, "SVTYPE": svtype, "dSV_Count": count}
        for (change, svtype), count in sorted(reid_type.items())
    ]
    reid_layer = Counter((r["African35_Change_Class"], r.get("Evidence_Layer", "NA")) for r in african_reidentified)
    reid_layer_rows = [
        {"African35_Change_Class": change, "Evidence_Layer": layer, "dSV_Count": count}
        for (change, layer), count in sorted(reid_layer.items())
    ]

    # Old-versus-hap38 comparison.
    old_type = Counter(r.get("SVTYPE", "NA") for r in old)
    new_type = Counter(r.get("SVTYPE", "NA") for r in formal)
    old_new_type_rows = []
    for svtype in sorted(set(old_type) | set(new_type)):
        old_new_type_rows.append({
            "SVTYPE": svtype,
            "Old_v5_dSV_Count": old_type.get(svtype, 0),
            "Hap38_dSV_Count": new_type.get(svtype, 0),
            "Delta_Hap38_Minus_Old": new_type.get(svtype, 0) - old_type.get(svtype, 0),
            "Retention_Fraction": f"{new_type.get(svtype, 0) / old_type.get(svtype, 0):.8f}" if old_type.get(svtype, 0) else "NA",
        })
    old_sample = Counter(s for r in old for s in split_samples(r.get("Samples")))
    old_new_sample_rows = []
    for row in sample_rows:
        sample = row["Sample_ID"]
        comparator = row["Old_Comparator"]
        directly_comparable = comparator != "NA" and row["Sample_Update_Mode"] == "Reused30"
        old_count = old_sample.get(comparator, 0) if comparator != "NA" else "NA"
        new_count = sample_counts.get(sample, 0)
        old_new_sample_rows.append({
            "Hap38_Sample": sample,
            "Old_Comparator": comparator,
            "Comparison_Status": "Direct_reused_haplotype" if directly_comparable else ("Collapsed_old_sample_not_one_to_one" if comparator != "NA" else "No_old_counterpart"),
            "Old_v5_dSV_Count": old_count,
            "Hap38_dSV_Count": new_count,
            "Delta_Hap38_Minus_Old": new_count - old_count if directly_comparable else "NA",
            "Population_Scope": group_by_sample.get(sample, "NA"),
        })
    nearest_detail, old_matched = nearest_matches(old, formal, 500)
    reverse_detail, new_matched = nearest_matches(formal, old, 500)
    spatial_summary = [
        {"Direction": "Old_to_Hap38", "Query_Count": len(old), "Matched_Count": old_matched, "Match_Fraction": f"{old_matched/len(old):.8f}", "Rule": "Same chromosome and SVTYPE; midpoint within 500 bp"},
        {"Direction": "Hap38_to_Old", "Query_Count": len(formal), "Matched_Count": new_matched, "Match_Fraction": f"{new_matched/len(formal):.8f}", "Rule": "Same chromosome and SVTYPE; midpoint within 500 bp"},
    ]

    non_african_sample_summary = []
    for sample in sorted(non_african):
        rows = [r for r in formal if sample in split_samples(r.get("Samples"))]
        non_african_sample_summary.append({
            "Sample_ID": sample,
            "Formal_All38_dSV_Count": len(rows),
            "DEL": sum(r.get("SVTYPE") == "DEL" for r in rows),
            "INS": sum(r.get("SVTYPE") == "INS" for r in rows),
            "DUP": sum(r.get("SVTYPE") == "DUP" for r in rows),
            "INV": sum(r.get("SVTYPE") == "INV" for r in rows),
            "Tier1_cuteSV_AND_assembly": sum(r.get("Evidence_Layer") == "Tier1_cuteSV_AND_assembly" for r in rows),
            "SyRI_large_rearrangement": sum(r.get("Evidence_Layer") == "SyRI_large_rearrangement" for r in rows),
        })

    # Integrity checks.
    checks = []
    def check(name, observed, expected, passed, note="OK"):
        checks.append({"Check_Name": name, "Observed_Value": observed, "Expected_Value": expected, "Status": "Pass" if passed else "Fail", "Note": note})
    check("All38_Manifest_Count", len(all38), EXPECTED_ALL38, len(all38) == EXPECTED_ALL38)
    check("African35_Manifest_Count", len(african), EXPECTED_AFRICAN35, len(african) == EXPECTED_AFRICAN35)
    check("Non_African_Manifest_Count", len(non_african), EXPECTED_NON_AFRICAN, len(non_african) == EXPECTED_NON_AFRICAN, ";".join(sorted(non_african)))
    check("Formal_dSV_Count", len(formal), EXPECTED_FORMAL_DSV, len(formal) == EXPECTED_FORMAL_DSV)
    check("Formal_Unique_SV_Key", len(formal_keys), len(formal), len(formal_keys) == len(formal))
    check("Unknown_Formal_Carriers", len(unknown_carriers), 0, not unknown_carriers, ";".join(sorted(unknown_carriers)) or "OK")
    check("Formal_Scope_Total", sum(scope_counts.values()), len(formal), sum(scope_counts.values()) == len(formal))
    check("Fixed_African_Plus_NonAfrican_Total", len(fixed_keys | {r.get('SV_Key') for r in non_african_detail}), len(formal), len(fixed_keys | {r.get('SV_Key') for r in non_african_detail}) == len(formal))
    check("African35_Reidentified_Unique_Key", len(reid_keys), len(african_reidentified), len(reid_keys) == len(african_reidentified))
    check("African35_Rule_Violations", african_rule_violations, 0, african_rule_violations == 0)
    observed_non_african = sum(any(s in non_african for s in split_samples(r.get("Samples_African35"))) for r in african_reidentified)
    check("African35_Reidentified_NonAfrican_Carriers", observed_non_african, 0, observed_non_african == 0)
    check("Evidence_Layer_Scope_Total", sum(scope_layer_counts.values()), len(formal), sum(scope_layer_counts.values()) == len(formal))
    if any(r["Status"] == "Fail" for r in checks):
        write_tsv(outdir / "qa/final_integrity_summary.tsv", checks, ["Check_Name", "Observed_Value", "Expected_Value", "Status", "Note"])
        raise SystemExit("Scope QC integrity gate failed")

    write_tsv(outdir / "panel_sample_manifest.tsv", sample_rows, ["Sample_ID", "Population_Scope", "Include_All38", "Include_African35", "Formal_All38_dSV_Count", "Formal_All38_dSV_Total_Length_bp", "Old_Comparator", "Sample_Update_Mode"])
    write_tsv(outdir / "all38_formal_dsv_by_scope.tsv", scope_summary, ["Carrier_Scope", "Formal_All38_dSV_Count", "Fraction_Of_All38_Formal"])
    write_tsv(outdir / "all38_formal_dsv_by_scope_svtype.tsv", by_type, ["Carrier_Scope", "SVTYPE", "Formal_All38_dSV_Count"])
    write_tsv(outdir / "all38_formal_dsv_by_scope_evidence_layer.tsv", by_layer, ["Carrier_Scope", "Evidence_Layer", "Formal_All38_dSV_Count"])
    write_tsv(outdir / "all38_formal_dsv_by_scope_functional_evidence.tsv", by_function, ["Carrier_Scope", "Functional_Evidence_Class", "Formal_All38_dSV_Count"])
    write_tsv(outdir / "all38_formal_african35_carrier_subset.tsv", african_fixed_subset, formal_fields)
    write_tsv(outdir / "non_african_hybrid_formal_dsv_detail.tsv", non_african_detail, formal_fields + ["Carrier_Scope"])
    write_tsv(outdir / "non_african_hybrid_summary_by_sample.tsv", non_african_sample_summary, ["Sample_ID", "Formal_All38_dSV_Count", "DEL", "INS", "DUP", "INV", "Tier1_cuteSV_AND_assembly", "SyRI_large_rearrangement"])
    write_tsv(outdir / "african35_reidentified_candidates.tsv", african_reidentified, catalog_fields + ["Samples_African35", "Sample_Count_African35", "Frequency_African35", "African35_dSV_Flag", "African35_Change_Class"])
    write_tsv(outdir / "african35_reidentification_summary.tsv", reid_summary, ["Metric", "Value", "Definition"])
    write_tsv(outdir / "african35_reidentified_by_svtype.tsv", reid_type_rows, ["African35_Change_Class", "SVTYPE", "dSV_Count"])
    write_tsv(outdir / "african35_reidentified_by_evidence_layer.tsv", reid_layer_rows, ["African35_Change_Class", "Evidence_Layer", "dSV_Count"])
    write_tsv(outdir / "old_vs_hap38_dsv_by_svtype.tsv", old_new_type_rows, ["SVTYPE", "Old_v5_dSV_Count", "Hap38_dSV_Count", "Delta_Hap38_Minus_Old", "Retention_Fraction"])
    write_tsv(outdir / "old_vs_hap38_sample_burden.tsv", old_new_sample_rows, ["Hap38_Sample", "Old_Comparator", "Comparison_Status", "Old_v5_dSV_Count", "Hap38_dSV_Count", "Delta_Hap38_Minus_Old", "Population_Scope"])
    write_tsv(outdir / "old_to_hap38_nearest_dsv_match.tsv", nearest_detail, ["Old_SV_Key", "Old_SV_ID", "Chrom", "SVTYPE", "Old_Position_bp", "Nearest_New_SV_Key", "Nearest_New_SV_ID", "Nearest_New_Position_bp", "Distance_bp", "Matched_Same_Chrom_Type_Within_500bp"])
    write_tsv(outdir / "old_vs_hap38_spatial_overlap_summary.tsv", spatial_summary, ["Direction", "Query_Count", "Matched_Count", "Match_Fraction", "Rule"])
    write_tsv(outdir / "qa/final_integrity_summary.tsv", checks, ["Check_Name", "Observed_Value", "Expected_Value", "Status", "Note"])

    metadata = {
        "Analysis": "All38 formal dSV scope audit and African35 sensitivity analysis",
        "Formal_Panel": "All38",
        "Sensitivity_Panel": "African35",
        "Non_African_or_Hybrid_Samples": sorted(non_african),
        "Rare_Threshold": RARE_THRESHOLD,
        "Rules": ["DEL/INS/DUP/INV", "Phoenix ALT_Derived", "CDS or conserved-region proxy", "panel frequency <= 0.05"],
        "Inputs": {"Catalog": args.catalog, "Candidates": args.candidates, "Manifest": args.manifest, "Old_Candidates": args.old_candidates},
        "Important_Distinction": {
            "Fixed_Subset": "All38 formal dSVs whose carriers are in African35",
            "Reidentified": "Full catalog re-filtered after carrier restriction and denominator change to 35",
        },
    }
    (outdir / "notes/analysis_metadata.json").write_text(json.dumps(metadata, indent=2, ensure_ascii=False) + "\n")
    notes = f"""# All38/African35 dSV scope analysis\n\nAll38 remains the formal panel ({len(formal):,} dSVs). African35 is a sensitivity analysis and does not replace the formal catalogue.\n\nTwo outputs answer different questions:\n\n1. `all38_formal_african35_carrier_subset.tsv` keeps the formal All38 definition fixed and selects formal dSVs carried within African35.\n2. `african35_reidentified_candidates.tsv` returns to the full catalogue, removes the three American-or-hybrid haplotypes, recalculates carrier count/frequency with denominator 35 and reapplies the same formal v5 rule.\n\nThe three non-African/hybrid haplotypes are {', '.join(sorted(non_african))}. They are not automatically labelled technical outliers. Scope differences may reflect panel composition, divergence and/or technical effects and require locus-level validation.\n\nEvidence layers remain separated. Conserved-region support is a proxy and must not be called GERP/phyloP constraint.\n"""
    (outdir / "notes/analysis_notes.md").write_text(notes)
    print(f"[OK] All38 formal={len(formal)}; fixed African35 subset={len(african_fixed_subset)}")
    print(f"[OK] African35 reidentified={len(african_reidentified)}; newly qualifying={len(newly_keys)}")
    print(f"[OK] Old-to-hap38 spatial matches within 500 bp={old_matched}/{len(old)}")


if __name__ == "__main__":
    main()
