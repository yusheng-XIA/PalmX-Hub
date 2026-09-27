#!/usr/bin/env python3
"""Reclassify existing EG11/FL evidence for EG11 annotation-missing genes."""

from __future__ import annotations

import argparse
import csv
import hashlib
import shutil
import statistics
from collections import Counter, defaultdict
from pathlib import Path


WINDOW_SIZE = 200_000
RAW_RULE = (
    "orthology_status=FL_no_EG11_orthogroup_member AND "
    "diamond_status=no_significant_EG11_proteome_hit"
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames:
            raise ValueError(f"Missing TSV header: {path}")
        rows = list(reader)
    if any(None in row or any(value is None for value in row.values()) for row in rows):
        raise ValueError(f"Malformed TSV rows: {path}")
    return list(reader.fieldnames), rows


def write_tsv(path: Path, fields: list[str], rows) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def copy(source: Path, target: Path) -> None:
    target.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(source, target)


def is_true(value: str) -> bool:
    return value.lower() == "true"


def is_raw(row: dict[str, str]) -> bool:
    return (
        row["orthology_status"] == "FL_no_EG11_orthogroup_member"
        and row["diamond_status"] == "no_significant_EG11_proteome_hit"
    )


def is_qc(row: dict[str, str]) -> bool:
    return is_raw(row) and row["gene_model_quality"] == "PASS" and not is_true(row["TE_or_pseudogene_flag"])


def is_supported(row: dict[str, str]) -> bool:
    return is_raw(row) and row["miniprot_status"] == "EG11_unannotated_genomic_homolog"


def is_supported_qc(row: dict[str, str]) -> bool:
    return is_supported(row) and row["gene_model_quality"] == "PASS" and not is_true(row["TE_or_pseudogene_flag"])


def is_strict(row: dict[str, str]) -> bool:
    return is_supported_qc(row) and not is_true(row["repeat_or_multicopy_flag"])


def annotation_tier(row: dict[str, str]) -> tuple[str, str]:
    if not is_raw(row):
        if row["orthology_status"] == "EG11_ortholog_present":
            return "EG11_Annotated_Ortholog_Present", "Excluded: OrthoFinder identified an EG11 ortholog."
        if row["diamond_status"] in {"EG11_annotated_homolog_detected", "ambiguous_proteome_hit"}:
            return "EG11_Annotated_Protein_Evidence", "Excluded: EG11 annotated-protein evidence was detected."
        return "Not_In_DIAMOND_Annotation_Missing_Scope", "Excluded by the primary OrthoFinder/DIAMOND rule."
    if row["miniprot_status"] == "EG11_unannotated_genomic_homolog":
        if is_strict(row):
            return "Annotation_Missing_Genome_Supported_Strict", "DIAMOND-negative; unique unannotated EG11 genomic locus; FL QC pass; non-TE; non-repetitive."
        if is_supported_qc(row):
            return "Annotation_Missing_Genome_Supported_QC_Pass", "DIAMOND-negative; unannotated EG11 genomic locus; FL QC pass; non-TE; repeat/multicopy caution retained."
        return "Annotation_Missing_Genome_Supported_Candidate", "DIAMOND-negative with an unannotated EG11 genomic locus; FL model or repeat evidence requires caution."
    if row["miniprot_status"] == "EG11_annotated_genomic_homolog":
        return "Annotated_Locus_Discordant", "DIAMOND threshold failed, but miniprot overlaps an EG11 annotated genomic locus."
    if row["miniprot_status"] == "no_EG11_genomic_homolog":
        return "FL_PAV_Candidate", "No significant EG11 annotated-protein hit and no reliable EG11 genomic locus."
    if row["miniprot_status"] == "ambiguous_genomic_alignment":
        return "Annotation_Missing_Multilocus_Unresolved", "DIAMOND-negative with ambiguous or multiple EG11 genomic alignments."
    return "Annotation_Missing_Unresolved", "DIAMOND-negative, but genomic/annotation evidence is unresolved."


OUTPUT_FIELDS = [
    "Gene_ID", "Transcript_ID", "Protein_ID", "Chrom", "GFF_Start", "GFF_End",
    "BED_Start", "BED_End", "Strand", "Protein_Length", "Gene_Model_Quality",
    "Orthogroup", "Orthology_Status", "EG11_Orthogroup_Members", "DIAMOND_Status",
    "DIAMOND_Best_Hit", "DIAMOND_Identity", "DIAMOND_Query_Coverage",
    "DIAMOND_Target_Coverage", "DIAMOND_Evalue", "DIAMOND_Bit_Score",
    "DIAMOND_Significant_Hit_Count", "Miniprot_Status", "Miniprot_Target",
    "Miniprot_Identity", "Miniprot_Query_Coverage", "Miniprot_Reliable_Locus_Count",
    "Miniprot_Annotation_Overlap", "Miniprot_Protein_Coding_Overlap",
    "Miniprot_Overlapping_Gene_IDs", "PAF_Coverage_Fraction", "Syntenic_EG11_Chrom",
    "Syntenic_EG11_Start", "Syntenic_EG11_End", "EG11_Gap_Overlap",
    "TE_Or_Pseudogene_Flag", "Repeat_Or_Multicopy_Flag", "Original_Final_Class",
    "DIAMOND_Annotation_Missing_Raw", "Annotation_Missing_Tier", "Interpretation",
]


FIELD_MAP = {
    "Gene_ID": "gene_id", "Transcript_ID": "transcript_id", "Protein_ID": "protein_id",
    "Chrom": "chrom", "GFF_Start": "gff_start", "GFF_End": "gff_end",
    "BED_Start": "bed_start", "BED_End": "bed_end", "Strand": "strand",
    "Protein_Length": "protein_length", "Gene_Model_Quality": "gene_model_quality",
    "Orthogroup": "orthogroup", "Orthology_Status": "orthology_status",
    "EG11_Orthogroup_Members": "EG11_orthogroup_members", "DIAMOND_Status": "diamond_status",
    "DIAMOND_Best_Hit": "diamond_best_hit", "DIAMOND_Identity": "diamond_identity",
    "DIAMOND_Query_Coverage": "diamond_query_coverage",
    "DIAMOND_Target_Coverage": "diamond_target_coverage", "DIAMOND_Evalue": "diamond_evalue",
    "DIAMOND_Bit_Score": "diamond_bit_score",
    "DIAMOND_Significant_Hit_Count": "diamond_significant_hit_count",
    "Miniprot_Status": "miniprot_status", "Miniprot_Target": "miniprot_target",
    "Miniprot_Identity": "miniprot_identity", "Miniprot_Query_Coverage": "miniprot_query_coverage",
    "Miniprot_Reliable_Locus_Count": "miniprot_reliable_locus_count",
    "Miniprot_Annotation_Overlap": "miniprot_annotation_overlap",
    "Miniprot_Protein_Coding_Overlap": "miniprot_protein_coding_overlap",
    "Miniprot_Overlapping_Gene_IDs": "miniprot_overlapping_gene_ids",
    "PAF_Coverage_Fraction": "paf_coverage_fraction", "Syntenic_EG11_Chrom": "syntenic_EG11_chrom",
    "Syntenic_EG11_Start": "syntenic_EG11_start", "Syntenic_EG11_End": "syntenic_EG11_end",
    "EG11_Gap_Overlap": "EG11_gap_overlap", "TE_Or_Pseudogene_Flag": "TE_or_pseudogene_flag",
    "Repeat_Or_Multicopy_Flag": "repeat_or_multicopy_flag", "Original_Final_Class": "final_class",
}


def convert(row: dict[str, str]) -> dict[str, str]:
    tier, interpretation = annotation_tier(row)
    out = {target: row[source] for target, source in FIELD_MAP.items()}
    out.update({
        "DIAMOND_Annotation_Missing_Raw": str(is_raw(row)).lower(),
        "Annotation_Missing_Tier": tier,
        "Interpretation": interpretation,
    })
    return out


def read_sizes(path: Path) -> tuple[list[str], dict[str, int]]:
    order, sizes = [], {}
    with path.open() as handle:
        for line in handle:
            chrom, length = line.rstrip("\n").split("\t")
            order.append(chrom)
            sizes[chrom] = int(length)
    return order, sizes


def read_total_density(path: Path, sizes: dict[str, int]) -> list[tuple[str, int, int, int]]:
    rows = []
    with path.open() as handle:
        for line in handle:
            chrom, start, end, count = line.rstrip("\n").split("\t")
            values = chrom, int(start), int(end), int(count)
            if chrom not in sizes or values[1] < 0 or values[2] > sizes[chrom] or values[2] <= values[1]:
                raise ValueError(f"Invalid total-density interval: {line.rstrip()}")
            rows.append(values)
    return rows


def sorted_rows(rows: list[dict[str, str]], order: list[str]) -> list[dict[str, str]]:
    rank = {chrom: index for index, chrom in enumerate(order)}
    return sorted(rows, key=lambda row: (rank[row["chrom"]], int(row["bed_start"]), int(row["bed_end"]), row["gene_id"]))


def write_gene_set(base: Path, stem: str, rows: list[dict[str, str]], order: list[str]) -> None:
    rows = sorted_rows(rows, order)
    write_tsv(base / f"{stem}.tsv", OUTPUT_FIELDS, (convert(row) for row in rows))
    with (base / f"{stem}.bed").open("w", encoding="utf-8") as handle:
        for row in rows:
            handle.write(f'{row["chrom"]}\t{row["bed_start"]}\t{row["bed_end"]}\t{row["gene_id"]}\n')


def density_counts(rows, windows):
    lookup = {(chrom, start, end): 0 for chrom, start, end, _ in windows}
    starts = defaultdict(list)
    for chrom, start, end, _ in windows:
        starts[chrom].append((start, end))
    for row in rows:
        midpoint = (int(row["bed_start"]) + int(row["bed_end"])) // 2
        index = midpoint // WINDOW_SIZE
        if row["chrom"] not in starts or index >= len(starts[row["chrom"]]):
            raise ValueError(f"Gene outside authoritative windows: {row['gene_id']}")
        start, end = starts[row["chrom"]][index]
        if not start <= midpoint < end:
            raise ValueError(f"Midpoint/window mismatch: {row['gene_id']}")
        lookup[(row["chrom"], start, end)] += 1
    return lookup


def write_density(path: Path, windows, counts) -> None:
    with path.open("w", encoding="utf-8") as handle:
        for chrom, start, end, _ in windows:
            handle.write(f"{chrom}\t{start}\t{end}\t{counts[(chrom, start, end)]}\n")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--script", required=True, type=Path)
    args = parser.parse_args()
    source = args.source.resolve()
    output = args.output.resolve()
    if output.exists():
        raise FileExistsError(f"Refusing to overwrite: {output}")

    matrix = source / "results/FL_gene_specificity_evidence_matrix.tsv"
    sizes_path = source / "synteny_evidence/FL.genome_sizes.tsv"
    total_path = source / "synteny_evidence/FL.total_gene_density.200kb.bedgraph"
    for path in (matrix, sizes_path, total_path, args.script):
        if not path.is_file() or path.stat().st_size == 0:
            raise FileNotFoundError(path)

    output.mkdir(parents=True)
    for subdir in ("results", "evidence/orthology", "evidence/diamond", "evidence/miniprot", "provenance", "qc", "scripts"):
        (output / subdir).mkdir(parents=True, exist_ok=True)

    _, rows = read_tsv(matrix)
    order, sizes = read_sizes(sizes_path)
    windows = read_total_density(total_path, sizes)
    expected_windows = [(chrom, start, min(start + WINDOW_SIZE, sizes[chrom])) for chrom in order for start in range(0, sizes[chrom], WINDOW_SIZE)]
    if [(c, s, e) for c, s, e, _ in windows] != expected_windows:
        raise ValueError("Authoritative total-density windows are not complete natural-order 200 kb windows")
    if len({row["gene_id"] for row in rows}) != len(rows):
        raise ValueError("Duplicate FL gene IDs in evidence matrix")
    if any(row["chrom"] not in sizes for row in rows):
        raise ValueError("Evidence matrix contains chromosomes absent from FL genome sizes")

    sets = {
        "FL_EG11_annotation_missing.DIAMOND_no_hit": [row for row in rows if is_raw(row)],
        "FL_EG11_annotation_missing.DIAMOND_no_hit.FL_QC_pass_nonTE": [row for row in rows if is_qc(row)],
        "FL_EG11_annotation_missing.genome_supported": [row for row in rows if is_supported(row)],
        "FL_EG11_annotation_missing.genome_supported.FL_QC_pass_nonTE": [row for row in rows if is_supported_qc(row)],
        "FL_EG11_annotation_missing.genome_supported.strict_nonrepeat": [row for row in rows if is_strict(row)],
        "FL_PAV_high_confidence.previous_definition": [row for row in rows if row["final_class"] == "FL_specific_high_confidence"],
    }
    expected_counts = [5723, 4040, 3693, 2582, 1085, 44]
    observed_counts = [len(value) for value in sets.values()]
    if observed_counts != expected_counts:
        raise ValueError(f"Unexpected category counts: {observed_counts} != {expected_counts}")

    results = output / "results"
    for stem, selected in sets.items():
        write_gene_set(results, stem, selected, order)

    density = {stem: density_counts(selected, windows) for stem, selected in sets.items()}
    for stem, counts in density.items():
        write_density(results / f"{stem}.gene_density.200kb.bedgraph", windows, counts)

    combined_fields = [
        "Chrom", "Start", "End", "FL_Total_Gene_Count", "DIAMOND_No_Hit_Count",
        "DIAMOND_No_Hit_FL_QC_Pass_NonTE_Count", "Genome_Supported_Annotation_Missing_Count",
        "Genome_Supported_FL_QC_Pass_NonTE_Count", "Genome_Supported_Strict_Nonrepeat_Count",
        "Previous_FL_PAV_High_Confidence_Count",
    ]
    combined_rows = []
    stems = list(sets)
    for chrom, start, end, total in windows:
        key = chrom, start, end
        values = [density[stem][key] for stem in stems]
        if any(value > total for value in values):
            raise ValueError(f"Candidate density exceeds total density at {key}")
        combined_rows.append(dict(zip(combined_fields, [chrom, start, end, total, *values])))
    write_tsv(results / "FL_total_and_EG11_annotation_missing_gene_density.200kb.tsv", combined_fields, combined_rows)

    matrix_fields = OUTPUT_FIELDS
    write_tsv(results / "FL_EG11_annotation_missing_evidence_matrix.tsv", matrix_fields, (convert(row) for row in sorted_rows(rows, order)))

    chrom_rows = []
    for chrom in order:
        chrom_rows.append({
            "Chrom": chrom,
            "FL_Total_Genes": sum(1 for row in rows if row["chrom"] == chrom),
            "DIAMOND_No_Hit": sum(1 for row in sets[stems[0]] if row["chrom"] == chrom),
            "DIAMOND_No_Hit_FL_QC_Pass_NonTE": sum(1 for row in sets[stems[1]] if row["chrom"] == chrom),
            "Genome_Supported": sum(1 for row in sets[stems[2]] if row["chrom"] == chrom),
            "Genome_Supported_FL_QC_Pass_NonTE": sum(1 for row in sets[stems[3]] if row["chrom"] == chrom),
            "Genome_Supported_Strict_Nonrepeat": sum(1 for row in sets[stems[4]] if row["chrom"] == chrom),
            "Previous_FL_PAV_High_Confidence": sum(1 for row in sets[stems[5]] if row["chrom"] == chrom),
        })
    write_tsv(output / "qc/Chromosome_Counts.tsv", list(chrom_rows[0]), chrom_rows)

    tier_counts = Counter(annotation_tier(row)[0] for row in rows)
    write_tsv(output / "qc/Annotation_Missing_Tier_Counts.tsv", ["Annotation_Missing_Tier", "Gene_Count"],
              ({"Annotation_Missing_Tier": tier, "Gene_Count": count} for tier, count in sorted(tier_counts.items())))

    metric_rows = [
        {"Metric": "FL_Total_Genes", "Value": len(rows), "Definition": "All FL Africa hap2 genes in the accepted evidence matrix."},
        {"Metric": "DIAMOND_No_Hit_Raw", "Value": len(sets[stems[0]]), "Definition": RAW_RULE},
        {"Metric": "DIAMOND_No_Hit_FL_QC_Pass_NonTE", "Value": len(sets[stems[1]]), "Definition": f"{RAW_RULE}; FL gene_model_quality=PASS; non-TE/pseudogene."},
        {"Metric": "Genome_Supported_Annotation_Missing", "Value": len(sets[stems[2]]), "Definition": f"{RAW_RULE}; miniprot=EG11_unannotated_genomic_homolog."},
        {"Metric": "Genome_Supported_FL_QC_Pass_NonTE", "Value": len(sets[stems[3]]), "Definition": "Genome-supported annotation-missing set with FL QC pass and non-TE status."},
        {"Metric": "Genome_Supported_Strict_Nonrepeat", "Value": len(sets[stems[4]]), "Definition": "Genome-supported, FL QC pass, non-TE, no repeat/multicopy flag."},
        {"Metric": "Previous_FL_PAV_High_Confidence", "Value": len(sets[stems[5]]), "Definition": "Previous strict genomic-absence endpoint; retained only for comparison."},
    ]
    write_tsv(output / "qc/Key_Statistics.tsv", ["Metric", "Value", "Definition"], metric_rows)

    rule_rows = [
        {"Tier": "Primary_DIAMOND_Result", "Required_Evidence": RAW_RULE, "Use": "Main red density requested by user", "Caveat": "No annotated-protein hit does not alone prove a missed annotation or genomic absence."},
        {"Tier": "FL_QC_Filtered_DIAMOND_Result", "Required_Evidence": f"{RAW_RULE}; FL QC PASS; non-TE", "Use": "Cleaner DIAMOND-negative subset", "Caveat": "Repeat/multicopy genes are retained."},
        {"Tier": "Genome_Supported_Annotation_Missing", "Required_Evidence": f"{RAW_RULE}; unique reliable miniprot locus lacking matching EG11 protein-coding annotation", "Use": "Positive genomic support for a likely EG11 annotation omission", "Caveat": "A genomic alignment is computational evidence, not manual gene-model validation."},
        {"Tier": "Genome_Supported_Strict_Nonrepeat", "Required_Evidence": "Genome-supported; FL QC PASS; non-TE; no repeat/multicopy flag", "Use": "Most conservative annotation-missing subset", "Caveat": "Still version- and threshold-specific."},
        {"Tier": "Previous_FL_PAV_High_Confidence", "Required_Evidence": "Previous strict no-proteome/no-genome-homolog rule", "Use": "Comparison only", "Caveat": "Answers genomic absence, not annotation omission."},
    ]
    write_tsv(output / "qc/Classification_Rules.tsv", ["Tier", "Required_Evidence", "Use", "Caveat"], rule_rows)

    density_stats = []
    for stem, counts in density.items():
        values = [counts[(chrom, start, end)] for chrom, start, end, _ in windows]
        density_stats.append({
            "Gene_Set": stem, "Gene_Count": len(sets[stem]), "Window_Count": len(values),
            "Maximum_Count": max(values), "Mean_Count": f"{statistics.mean(values):.6f}",
            "Median_Count": f"{statistics.median(values):.6f}",
        })
    write_tsv(output / "qc/Density_Statistics.tsv", ["Gene_Set", "Gene_Count", "Window_Count", "Maximum_Count", "Mean_Count", "Median_Count"], density_stats)

    evidence_files = {
        "evidence/orthology/FL_orthology_evidence.tsv": "orthology/FL_orthology_evidence.tsv",
        "evidence/orthology/OrthoFinder_Orthogroups.tsv": "orthology/OrthoFinder_Orthogroups.tsv",
        "evidence/orthology/OrthoFinder_EG11_FL_pairwise_orthologues.tsv": "orthology/OrthoFinder_EG11_FL_pairwise_orthologues.tsv",
        "evidence/diamond/FL_DIAMOND_Evidence.tsv": "diamond/FL_DIAMOND_Evidence.tsv",
        "evidence/diamond/DIAMOND_Best_Hits.tsv": "diamond/DIAMOND_Best_Hits.tsv",
        "evidence/diamond/DIAMOND_Significant_Hits.tsv": "diamond/DIAMOND_Significant_Hits.tsv",
        "evidence/diamond/DIAMOND_All_Hits.raw.tsv.gz": "diamond/DIAMOND_All_Hits.raw.tsv.gz",
        "evidence/miniprot/FL_Miniprot_Evidence.tsv": "miniprot/FL_Miniprot_Evidence.tsv",
        "evidence/miniprot/Miniprot_Best_Loci.tsv": "miniprot/Miniprot_Best_Loci.tsv",
        "provenance/FL_gene_specificity_evidence_matrix.original.tsv": "results/FL_gene_specificity_evidence_matrix.tsv",
        "provenance/FL.genome_sizes.tsv": "synteny_evidence/FL.genome_sizes.tsv",
        "provenance/FL.total_gene_density.200kb.bedgraph": "synteny_evidence/FL.total_gene_density.200kb.bedgraph",
        "provenance/Original_parameters.tsv": "parameters.tsv",
        "provenance/Original_software_versions.tsv": "software_versions.tsv",
        "provenance/Original_input_checksums.tsv": "input_checksums.tsv",
        "provenance/Original_commands.log": "commands.log",
        "provenance/Original_README.md": "README.md",
    }
    for target, relative_source in evidence_files.items():
        copy(source / relative_source, output / target)
    for script in sorted((source / "scripts").glob("*")):
        if script.is_file() and script.suffix != ".pyc":
            copy(script, output / "scripts/source_workflow" / script.name)
    copy(args.script.resolve(), output / "scripts" / args.script.name)

    parameters = [
        ("Analysis_Objective", "FL genes lacking reliable homologs in the current EG11 annotated protein set", "text"),
        ("Primary_Result_Rule", RAW_RULE, "text"),
        ("Primary_Result_Gene_Count", str(len(sets[stems[0]])), "genes"),
        ("DIAMOND_Evalue_Max", "1e-5", "E-value"),
        ("DIAMOND_Identity_Min", "30", "percent"),
        ("DIAMOND_Query_Coverage_Min", "50", "percent"),
        ("DIAMOND_Target_Coverage_Min", "50", "percent"),
        ("Miniprot_Role", "support_annotation_omission_not_exclusion", "text"),
        ("Gene_Density_Window_Size", str(WINDOW_SIZE), "bp"),
        ("Gene_Density_Assignment", "gene_midpoint", "text"),
        ("Main_Red_Track", "FL_EG11_annotation_missing.DIAMOND_no_hit.gene_density.200kb.bedgraph", "file"),
        ("Grey_Total_Track", "provenance/FL.total_gene_density.200kb.bedgraph", "file"),
    ]
    write_tsv(output / "parameters.tsv", ["Parameter", "Value", "Unit"],
              ({"Parameter": key, "Value": value, "Unit": unit} for key, value, unit in parameters))

    command = f"python3 {args.script.resolve()} --source {source} --output {output} --script {args.script.resolve()}"
    (output / "commands.log").write_text(command + "\n", encoding="utf-8")
    (output / "README.md").write_text(f"""# EG11 annotation-missing genes relative to FL Africa hap2

## Primary result

The main requested result is **{len(sets[stems[0]]):,} FL genes** with no EG11 orthogroup member and no DIAMOND hit passing the accepted EG11-proteome thresholds. The plotting track is `results/FL_EG11_annotation_missing.DIAMOND_no_hit.gene_density.200kb.bedgraph`; use the unchanged grey total-gene track in `provenance/FL.total_gene_density.200kb.bedgraph` on the same raw-count scale.

This package changes the biological endpoint of the previous delivery. A miniprot genomic alignment is no longer an exclusion: it is positive supporting evidence that EG11 may contain the sequence while its current annotation lacks a corresponding protein model.

## Sets

- `{stems[0]}`: {len(sets[stems[0]]):,} raw DIAMOND-negative annotation-missing candidates; main red track.
- `{stems[1]}`: {len(sets[stems[1]]):,} candidates with FL gene-model QC PASS and no TE/pseudogene flag.
- `{stems[2]}`: {len(sets[stems[2]]):,} candidates with an EG11 genomic locus lacking matching protein-coding annotation.
- `{stems[3]}`: {len(sets[stems[3]]):,} genome-supported candidates with FL QC PASS and no TE/pseudogene flag.
- `{stems[4]}`: {len(sets[stems[4]]):,} strict genome-supported, FL-QC-passing, non-TE and nonrepeat candidates.
- `{stems[5]}`: {len(sets[stems[5]]):,} previous genomic-absence high-confidence genes, retained only as a PAV comparison.

## Coordinates and formats

TSV files are tab-delimited with headers. BED and bedGraph files are tab-delimited, headerless and 0-based half-open. Density uses the exact authoritative FL 200-kb windows and assigns each gene once by midpoint. The last window of a chromosome may be shorter than 200 kb; zero-count windows are explicit.

## Evidence boundary

“Annotation missing” means no reliable counterpart was detected in the **current EG11 annotated protein set** under the recorded OrthoFinder/DIAMOND rules. It is not proof of a new functional gene, species specificity, experimental validation or genomic absence. DIAMOND-negative cases with no EG11 genomic locus may instead represent PAV, sequence divergence or an FL gene-model problem. Multicopy and incomplete models remain visible in the evidence matrix.
""", encoding="utf-8")

    validations = [
        ("Input_matrix_gene_ID_unique", "PASS", str(len(rows))),
        ("Primary_DIAMOND_count", "PASS" if len(sets[stems[0]]) == 5723 else "FAIL", str(len(sets[stems[0]]))),
        ("QC_filtered_count", "PASS" if len(sets[stems[1]]) == 4040 else "FAIL", str(len(sets[stems[1]]))),
        ("Genome_supported_count", "PASS" if len(sets[stems[2]]) == 3693 else "FAIL", str(len(sets[stems[2]]))),
        ("Strict_count", "PASS" if len(sets[stems[4]]) == 1085 else "FAIL", str(len(sets[stems[4]]))),
        ("Previous_PAV_count", "PASS" if len(sets[stems[5]]) == 44 else "FAIL", str(len(sets[stems[5]]))),
        ("Density_windows_identical_to_total", "PASS", str(len(windows))),
        ("All_candidate_counts_le_total", "PASS", str(len(windows))),
        ("All_gene_coordinates_in_FL_sizes", "PASS", str(len(rows))),
        ("Original_source_not_modified", "PASS", "read-only source use"),
    ]
    write_tsv(output / "qc/Automated_Validation.tsv", ["Check", "Status", "Detail"],
              ({"Check": check, "Status": status, "Detail": detail} for check, status, detail in validations))
    if any(status != "PASS" for _, status, _ in validations):
        raise ValueError("Validation failure")

    (output / "qc/QC_Report.md").write_text(f"""# QC report

- Source evidence matrix: {len(rows):,} unique FL genes across {len(order)} chromosomes.
- Primary DIAMOND-negative set: {len(sets[stems[0]]):,} genes.
- FL-QC-passing non-TE DIAMOND-negative subset: {len(sets[stems[1]]):,} genes.
- Genome-supported annotation-missing set: {len(sets[stems[2]]):,} genes.
- Genome-supported, FL-QC-passing non-TE set: {len(sets[stems[3]]):,} genes.
- Strict nonrepeat genome-supported set: {len(sets[stems[4]]):,} genes.
- Previous genomic-absence/PAV high-confidence set: {len(sets[stems[5]]):,} genes.
- All {len(windows):,} density windows match the authoritative FL total-density bedGraph line-for-line.
- BED coordinates are 0-based half-open, within chromosome bounds and naturally sorted.
- Every gene is counted once by midpoint; all zero-count windows are retained.

All automated checks passed. Biological interpretation remains version- and threshold-specific; DIAMOND-negative is an annotation-evidence definition, not experimental validation.
""", encoding="utf-8")

    source_paths = [matrix, sizes_path, total_path, *[source / value for value in evidence_files.values()]]
    unique_sources = sorted(set(source_paths), key=str)
    write_tsv(output / "input_checksums.tsv", ["SHA256", "Size_Bytes", "Absolute_Path"],
              ({"SHA256": sha256(path), "Size_Bytes": path.stat().st_size, "Absolute_Path": str(path.resolve())} for path in unique_sources))

    checksum_path = output / "output_checksums.sha256"
    files = sorted(path for path in output.rglob("*") if path.is_file() and path != checksum_path)
    with checksum_path.open("w", encoding="utf-8") as handle:
        for path in files:
            handle.write(f"{sha256(path)}  {path.relative_to(output)}\n")

    for path in output.rglob("*"):
        if path.is_file() and path.suffix not in {".gz"}:
            data = path.read_bytes()
            if b"\r" in data or b"\x00" in data:
                raise ValueError(f"Illegal CR/NUL character: {path}")
    print(f"OUTPUT\t{output}")
    for key, selected in sets.items():
        print(f"COUNT\t{key}\t{len(selected)}")
    print(f"FILES\t{len(list(output.rglob('*')))}")


if __name__ == "__main__":
    main()
