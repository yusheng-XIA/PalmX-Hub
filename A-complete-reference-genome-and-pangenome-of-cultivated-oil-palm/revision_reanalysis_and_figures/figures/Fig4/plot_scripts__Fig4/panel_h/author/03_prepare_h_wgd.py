#!/usr/bin/env python3
"""Prepare deterministic current-PAV WGD inputs and selected prefixed CDS."""

from __future__ import annotations

import csv
import json
import os
from pathlib import Path

os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
import numpy as np


RUN = Path("${ANALYSIS_DIR}/22_answer_reviews/00_ms/03_V3/04_figure4/Fig4_d_i_pan39_material33_redraw_20260808")
PAN = Path("${ANALYSIS_DIR}/14_pan_genome/11_new_pan")
OLD_CDS = Path("${ANALYSIS_DIR}/20_results/Figure3/02_pan_genome/03_ka_ks/cds_sequences")
ATTEMPT = os.environ.get("ATTEMPT_ID", "attempt1")
WORK = RUN / "work/h_wgd" / ATTEMPT
ORTHO = PAN / "results/orthofinder_attempt1/Results_Aug05/Orthogroups/Orthogroups.tsv"
FREQ = PAN / "results/postprocess_attempt1/Orthogroups.material33.freq_class.tsv"
MANIFEST = PAN / "config/sample_manifest.tsv"
ORDER = ("core", "softcore", "shell", "cloud")


def current_class(freq: int) -> str:
    if freq == 33:
        return "core"
    if freq == 32:
        return "softcore"
    if freq >= 2:
        return "shell"
    return "cloud"


def orthogroups_by_class() -> dict[str, list[str]]:
    grouped = {c: [] for c in ORDER}
    with FREQ.open() as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            grouped[current_class(int(row["material_frequency"]))].append(row["Orthogroup"])
    return grouped


def write_mcl(grouped: dict[str, list[str]]) -> tuple[set[str], dict[str, int], dict[str, int], dict[str, int]]:
    class_lookup = {og: cat for cat, ogs in grouped.items() for og in ogs}
    eligible: dict[str, list[list[str]]] = {c: [] for c in ORDER}
    excluded_over_cap = {c: 0 for c in ORDER}
    family_gene_cap = 700
    with ORTHO.open() as handle:
        reader = csv.reader(handle, delimiter="\t")
        next(reader)
        for row in reader:
            if not row or row[0] not in class_lookup:
                continue
            genes = []
            for cell in row[1:]:
                genes.extend(x.strip() for x in cell.split(",") if x.strip())
            category = class_lookup[row[0]]
            if len(genes) > family_gene_cap:
                excluded_over_cap[category] += 1
                continue
            if len(genes) >= 2:
                genes = [gene.replace("__", "_SAMPLESEP_") for gene in genes]
                eligible[category].append(genes)

    eligible_counts = {c: len(eligible[c]) for c in ORDER}
    handles = {c: (WORK / f"{c}.mcl").open("w") for c in ORDER}
    counts = {c: 0 for c in ORDER}
    needed: set[str] = set()
    try:
        for category in ORDER:
            values = eligible[category]
            if len(values) < 2000:
                raise RuntimeError(f"Only {len(values)} eligible families for {category}")
            rng = np.random.RandomState(42)
            idx = rng.choice(len(values), 2000, replace=False)
            for i in idx:
                genes = values[int(i)]
                handles[category].write("\t".join(genes) + "\n")
                counts[category] += 1
                needed.update(genes)
    finally:
        for handle in handles.values():
            handle.close()
    return needed, counts, eligible_counts, excluded_over_cap


def sources() -> list[tuple[str, Path]]:
    result = []
    with MANIFEST.open() as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            sample = row["sample_id"]
            if row["input_mode"] == "existing_protein":
                source = OLD_CDS / f"{sample}.cds.fa"
            else:
                source = PAN / "inputs/derived_new_haplotypes" / sample / "raw.cds.fa"
            if not source.is_file() or source.stat().st_size == 0:
                raise FileNotFoundError(f"Missing CDS for {sample}: {source}")
            result.append((sample, source))
    if len(result) != 39:
        raise AssertionError(f"Expected 39 CDS sources, found {len(result)}")
    return result


def extract_selected_cds(needed: set[str], source_list: list[tuple[str, Path]]) -> set[str]:
    out = WORK / "all_selected_cds.fa"
    found: set[str] = set()
    with out.open("w") as dest:
        for sample, source in source_list:
            keep = False
            final_id = ""
            with source.open() as handle:
                for line in handle:
                    if line.startswith(">"):
                        original = line[1:].split()[0]
                        # New gffread-derived CDS uses transcript IDs (evm.model),
                        # whereas the normalized Orthofinder input uses gene IDs
                        # (evm.TU). Retained CDS already uses evm.TU and is unchanged.
                        original = original.replace("evm.model", "evm.TU", 1)
                        final_id = f"{sample}_SAMPLESEP_{original}"
                        keep = final_id in needed
                        if keep:
                            found.add(final_id)
                            dest.write(f">{final_id}\n")
                    elif keep:
                        dest.write(line)
    return found


def main() -> None:
    WORK.mkdir(parents=True, exist_ok=True)
    grouped = orthogroups_by_class()
    needed, mcl_counts, eligible_counts, excluded_over_cap = write_mcl(grouped)
    source_list = sources()
    found = extract_selected_cds(needed, source_list)
    missing = sorted(needed - found)
    report = {
        "attempt": ATTEMPT,
        "class_total_orthogroups": {"core": 20736, "softcore": 2057, "shell": 23564, "cloud": 2563},
        "eligible_after_2_to_700_gene_filter": eligible_counts,
        "excluded_over_700_genes": excluded_over_cap,
        "mcl_families": mcl_counts,
        "required_gene_ids": len(needed),
        "found_cds_ids": len(found),
        "missing_cds_ids": len(missing),
        "unsafe_double_underscore_ids": sum("__" in gene for gene in needed),
        "wgd_pair_delimiter_escape": "__ -> _SAMPLESEP_",
        "cds_sources": len(source_list),
        "random_seed": 42,
        "max_orthogroups_per_class": 2000,
        "family_gene_count_filter": "2 <= genes per orthogroup <= 700",
        "historical_input_max_family_size": 621,
    }
    (WORK / "prepare_report.json").write_text(json.dumps(report, indent=2) + "\n")
    if missing:
        (WORK / "missing_cds_ids.txt").write_text("\n".join(missing) + "\n")
        raise RuntimeError(f"Missing {len(missing)} selected CDS IDs")
    if report["unsafe_double_underscore_ids"] != 0:
        raise RuntimeError("Unsafe double-underscore gene IDs remain for wgd --pairwise")
    if any(mcl_counts[c] == 0 for c in ORDER):
        raise RuntimeError(f"Empty MCL class: {mcl_counts}")
    (RUN / "provenance" / f"h_wgd_prepare.{ATTEMPT}.SUCCESS").touch()
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
