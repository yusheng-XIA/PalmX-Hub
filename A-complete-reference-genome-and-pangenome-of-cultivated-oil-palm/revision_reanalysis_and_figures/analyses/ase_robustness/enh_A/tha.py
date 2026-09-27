#!/usr/bin/env python3
"""Build trait-oriented ASE modules without making causal assignments."""

from __future__ import annotations

import argparse
from pathlib import Path
import re

import numpy as np
import pandas as pd


MODULE_ORDER = [
    "Oil biosynthesis & storage",
    "De-novo / saturated FA",
    "Unsaturated FA",
    "TAG assembly & oil body",
    "Lipid oxidation / antioxidant",
    "Shell / cell wall / lignin",
]


def parse_emapper(path: Path) -> pd.DataFrame:
    header = None
    rows = []
    with path.open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("#query"):
                header = line[1:].rstrip("\n").split("\t")
                continue
            if line.startswith("#") or not line.strip() or header is None:
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < len(header):
                fields += [""] * (len(header) - len(fields))
            r = dict(zip(header, fields))
            rows.append(
                {
                    "gene_africa": r.get("query", ""),
                    "preferred_name": r.get("Preferred_name", ""),
                    "description": r.get("Description", ""),
                    "kegg_ko": r.get("KEGG_ko", ""),
                    "kegg_pathway": r.get("KEGG_Pathway", ""),
                }
            )
    return pd.DataFrame(rows).drop_duplicates("gene_africa")


def to_tu(gene: str) -> str:
    gene = str(gene)
    return re.sub(r"^evm\.model\.", "evm.TU.", gene)


def build_catalog(fatty_path: Path, emapper_path: Path, shell_path: Path) -> pd.DataFrame:
    fatty = pd.read_csv(fatty_path, sep="\t")
    fatty = fatty[["Section", "Enzyme", "GeneID", "Preferred_name", "Description"]].drop_duplicates()
    fatty = fatty.rename(
        columns={
            "GeneID": "gene_africa",
            "Preferred_name": "preferred_name",
            "Description": "description",
        }
    )
    records = []
    saturated = {
        "ACCase", "ACP", "FATA/B", "FabD (MCAT)", "FabG (KAR)", "FabI (ENR)",
        "KAS I/II", "KASIII",
    }
    unsaturated = {"SAD", "FAD2", "FAD3", "FAD6", "FAD7/FAD8"}
    tag = {"LACS", "KCS", "DGAT", "GPAT", "LPAT", "PAP", "PDAT", "PDCT"}
    for row in fatty.itertuples(index=False):
        base = {
            "gene_africa": row.gene_africa,
            "family": row.Enzyme,
            "preferred_name": row.preferred_name,
            "description": row.description,
            "evidence_source": "integrated_fatty_acid_catalog",
        }
        records.append({**base, "trait_module": "Oil biosynthesis & storage"})
        if row.Enzyme in saturated:
            records.append({**base, "trait_module": "De-novo / saturated FA"})
        if row.Enzyme in unsaturated:
            records.append({**base, "trait_module": "Unsaturated FA"})
        if row.Enzyme in tag:
            records.append({**base, "trait_module": "TAG assembly & oil body"})

    annotations = parse_emapper(emapper_path)
    annotations["text"] = (
        annotations["preferred_name"].fillna("") + " " + annotations["description"].fillna("")
    ).str.lower()
    oxidation_pattern = re.compile(
        r"lipoxygenase|peroxidase|catalase|tocopherol|glutathione peroxidase|"
        r"oxidative stress|antioxidant|hydroperoxide"
    )
    shell_pattern = re.compile(
        r"lignin|cell wall|cellulose|pectin|phenylpropanoid|laccase|secondary wall|"
        r"suberin|seedstick|endocarp"
    )
    for module, pattern, family in (
        ("Lipid oxidation / antioxidant", oxidation_pattern, "oxidation/antioxidant"),
        ("Shell / cell wall / lignin", shell_pattern, "shell/cell-wall/lignin"),
    ):
        selected = annotations[annotations["text"].str.contains(pattern, regex=True, na=False)]
        for row in selected.itertuples(index=False):
            records.append(
                {
                    "gene_africa": to_tu(row.gene_africa),
                    "trait_module": module,
                    "family": row.preferred_name if row.preferred_name not in ("", "-") else family,
                    "preferred_name": row.preferred_name,
                    "description": row.description,
                    "evidence_source": "Africa_hap2_eggnog_keyword",
                }
            )

    shell = pd.read_csv(shell_path, sep="\t")
    shell = shell[shell["Genome"].isin(["Africa_hap2", "dura", "pisifera"])].copy()
    for row in shell.itertuples(index=False):
        if row.Genome == "Africa_hap2":
            records.append(
                {
                    "gene_africa": to_tu(row.Gene_ID),
                    "trait_module": "Shell / cell wall / lignin",
                    "family": "SEEDSTICK-like shell candidate",
                    "preferred_name": "SEEDSTICK-like",
                    "description": f"project shell candidate; identity={row.Identity}; coverage={row.Coverage}",
                    "evidence_source": "project_shell_candidate",
                }
            )
    catalog = pd.DataFrame(records).drop_duplicates(["gene_africa", "trait_module", "family"])
    catalog["trait_module"] = pd.Categorical(catalog["trait_module"], MODULE_ORDER, ordered=True)
    return catalog.sort_values(["trait_module", "gene_africa", "family"])


def summarize(detail: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    eligible = detail[detail["eligible"]].copy()
    eligible["robust_A"] = eligible["ase_call"] == "Allele_A_biased"
    eligible["robust_B"] = eligible["ase_call"] == "Allele_B_biased"
    eligible["robust_any"] = eligible["robust_A"] | eligible["robust_B"]
    module = (
        eligible.groupby(["analysis", "trait_module", "stage_group"], observed=True)
        .agg(
            tested_gene_stage_rows=("gene_id", "size"),
            tested_genes=("gene_id", "nunique"),
            allele_A_biased_rows=("robust_A", "sum"),
            allele_B_biased_rows=("robust_B", "sum"),
            robust_ASE_rows=("robust_any", "sum"),
            median_log2_ratio=("log2_allele_ratio", "median"),
        )
        .reset_index()
    )
    module["robust_ASE_percentage"] = (
        100 * module["robust_ASE_rows"] / module["tested_gene_stage_rows"].clip(lower=1)
    )
    robust_median = (
        eligible[eligible["robust_any"]]
        .groupby(["analysis", "trait_module", "stage_group"], observed=True)["log2_allele_ratio"]
        .median()
        .rename("robust_median_log2_ratio")
        .reset_index()
    )
    module = module.merge(
        robust_median, on=["analysis", "trait_module", "stage_group"], how="left"
    )

    family = (
        eligible.groupby(["analysis", "trait_module", "family", "stage_group"], observed=True)
        .agg(
            tested_rows=("gene_id", "size"),
            genes=("gene_id", "nunique"),
            robust_ASE_rows=("robust_any", "sum"),
            median_log2_ratio=("log2_allele_ratio", "median"),
        )
        .reset_index()
    )
    family["robust_ASE_percentage"] = 100 * family["robust_ASE_rows"] / family["tested_rows"].clip(lower=1)
    return module, family


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--gene-stage-ase", type=Path, required=True)
    ap.add_argument("--bridge", type=Path, required=True)
    ap.add_argument("--fatty-catalog", type=Path, required=True)
    ap.add_argument("--africa-emapper", type=Path, required=True)
    ap.add_argument("--shell-candidates", type=Path, required=True)
    ap.add_argument("--output-dir", type=Path, required=True)
    args = ap.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    ase = pd.read_csv(args.gene_stage_ase, sep="\t")
    bridge = pd.read_csv(args.bridge, sep="\t")
    catalog = build_catalog(args.fatty_catalog, args.africa_emapper, args.shell_candidates)
    catalog.to_csv(args.output_dir / "trait_gene_catalog.tsv", sep="\t", index=False)

    fl = ase[ase["analysis"] == "FL"].merge(
        catalog,
        left_on="gene_id",
        right_on="gene_africa",
        how="inner",
        validate="many_to_many",
    )
    tn_catalog = catalog.merge(
        bridge[["orthogroup", "gene_dura", "gene_africa"]],
        on="gene_africa",
        how="inner",
        validate="many_to_one",
    )
    tn = ase[ase["analysis"] == "TN"].merge(
        tn_catalog,
        left_on="gene_id",
        right_on="gene_dura",
        how="inner",
        validate="many_to_many",
    )

    # Add the project-curated Dura SEEDSTICK-like candidate directly.
    shell_dura = "evm.TU.chr01.1793"
    shell_rows = ase[(ase["analysis"] == "TN") & (ase["gene_id"] == shell_dura)].copy()
    if not shell_rows.empty:
        shell_rows["gene_africa"] = np.nan
        shell_rows["orthogroup"] = np.nan
        shell_rows["gene_dura"] = shell_dura
        shell_rows["trait_module"] = "Shell / cell wall / lignin"
        shell_rows["family"] = "SEEDSTICK-like shell candidate"
        shell_rows["preferred_name"] = "SEEDSTICK-like"
        shell_rows["description"] = "project Dura shell candidate"
        shell_rows["evidence_source"] = "project_shell_candidate"
        tn = pd.concat([tn, shell_rows], ignore_index=True, sort=False)

    detail = pd.concat([tn, fl], ignore_index=True, sort=False)
    detail["allele_A_label"] = np.where(
        detail["analysis"] == "TN", "Dura-like / TK-like", "Africa hap2"
    )
    detail["allele_B_label"] = np.where(
        detail["analysis"] == "TN", "Pisifera-like / NS-like", "American hap1"
    )
    detail["interpretation_guardrail"] = "haplotype_expression_support_not_causal_trait_assignment"
    detail = detail.drop_duplicates(
        ["analysis", "gene_id", "stage", "trait_module", "family"]
    ).sort_values(["analysis", "trait_module", "family", "stage_index", "gene_id"])
    detail.to_csv(
        args.output_dir / "trait_haplotype_ASE.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
        float_format="%.8g",
    )
    module, family = summarize(detail)
    module.to_csv(args.output_dir / "trait_module_ASE_summary.tsv", sep="\t", index=False, float_format="%.8g")
    family.to_csv(args.output_dir / "trait_family_ASE_summary.tsv", sep="\t", index=False, float_format="%.8g")
    print(f"trait_catalog_rows={len(catalog)}")
    print(f"trait_ASE_rows={len(detail)}")
    print(f"trait_modules={detail['trait_module'].nunique()}")


if __name__ == "__main__":
    main()

