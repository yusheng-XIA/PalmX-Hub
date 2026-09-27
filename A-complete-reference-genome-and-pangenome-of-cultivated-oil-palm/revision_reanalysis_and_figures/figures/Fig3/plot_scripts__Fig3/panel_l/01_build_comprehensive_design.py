#!/usr/bin/env python3
"""Build a full-gene primitive-ancestry breeding evidence catalog."""
from __future__ import annotations

import argparse
import csv
import math
import os
import re
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
import pandas as pd


ALLOWED = {"Dura": "D", "Pisifera": "P", "Meizhou4": "M", "Meizhou4_like": "M"}


# Ordered broad-recall rules. A match is a candidate invitation, not proof of direction.
RULES = [
    # Developmental timing, hormones, and transcriptional control.
    ("M01", "Fruit development and lipid-window timing", "WRI1", r"\bwrinkled\s*1\b|\bwri1\b", "timing"),
    ("M01", "Fruit development and lipid-window timing", "LEC/ABI3/FUS3", r"leafy cotyledon|\blec1\b|\blec2\b|abscisic acid insensitive 3|\babi3\b|\bfus3\b", "timing"),
    ("M01", "Fruit development and lipid-window timing", "AP2/ERF", r"ap2[- /]erf|ethylene[- ]responsive transcription factor|ethylene response factor", "timing"),
    ("M01", "Fruit development and lipid-window timing", "NAC", r"\bnac\b.*transcription|transcription factor nac", "timing"),
    ("M01", "Fruit development and lipid-window timing", "MADS-box", r"mads[- ]box|agamous[- ]like", "timing"),
    ("M01", "Fruit development and lipid-window timing", "bZIP", r"basic leucine zipper|\bbzip\b", "timing"),
    ("M01", "Fruit development and lipid-window timing", "MYB", r"\bmyb\b.*transcription|transcription factor myb", "timing"),
    ("M01", "Fruit development and lipid-window timing", "DOF", r"dna binding with one finger|dof transcription", "timing"),
    ("M01", "Fruit development and lipid-window timing", "NF-Y", r"nuclear transcription factor y|\bnf-y[abc]\b", "timing"),
    ("M01", "Fruit development and lipid-window timing", "ARF/Aux-IAA", r"auxin response factor|auxin-responsive|aux[/ -]?iaa", "timing"),
    ("M01", "Fruit development and lipid-window timing", "Ethylene signalling", r"ethylene insensitive|ein3|ein4|etr1|ethylene receptor|1-aminocyclopropane-1-carboxylate", "timing"),
    ("M01", "Fruit development and lipid-window timing", "ABA signalling", r"abscisic acid|aba receptor|pyrabactin resistance|snf1-related protein kinase 2|snrk2", "timing"),
    ("M01", "Fruit development and lipid-window timing", "Gibberellin signalling", r"gibberellin receptor|gibberellin response|gibberellin signalling", "timing"),
    ("M01", "Fruit development and lipid-window timing", "Sugar-energy signalling", r"snf1-related protein kinase 1|snrk1|target of rapamycin|\btor kinase\b|trehalose[- ]6[- ]phosphate", "timing"),
    ("M01", "Fruit development and lipid-window timing", "Ripening/circadian", r"ripening|circadian|timing of cab expression|pseudo-response regulator", "timing"),

    # Carbon import and precursor supply.
    ("M02", "Carbon import and acetyl-CoA supply", "SWEET/SUT", r"sugar will eventually be exported|\bsweet\b|sucrose transporter|sucrose-proton symporter|sucrose transport protein", "favorable"),
    ("M02", "Carbon import and acetyl-CoA supply", "Sucrose synthase", r"sucrose synthase|sucrose-cleaving enzyme", "favorable"),
    ("M02", "Carbon import and acetyl-CoA supply", "Invertase", r"invertase|beta-fructofuranosidase", "context"),
    ("M02", "Carbon import and acetyl-CoA supply", "Hexose transport", r"hexose transporter|glucose transporter|fructose transporter|monosaccharide transporter", "favorable"),
    ("M02", "Carbon import and acetyl-CoA supply", "HXK/FRK", r"hexokinase|fructokinase|galactokinase", "favorable"),
    ("M02", "Carbon import and acetyl-CoA supply", "Glycolysis", r"phosphofructokinase|fructose-bisphosphate aldolase|glyceraldehyde-3-phosphate dehydrogenase|phosphoglycerate kinase|phosphoglycerate mutase|enolase|pyruvate kinase", "favorable"),
    ("M02", "Carbon import and acetyl-CoA supply", "Glycolysis description", r"fructose 6-phosphate.*fructose 1,6-bisphosphate|first committing step of glycolysis|phosphopyruvate hydratase", "favorable"),
    ("M02", "Carbon import and acetyl-CoA supply", "Pyruvate dehydrogenase", r"pyruvate dehydrogenase|dihydrolipoyl transacetylase|dihydrolipoamide dehydrogenase", "favorable"),
    ("M02", "Carbon import and acetyl-CoA supply", "ATP-citrate lyase", r"atp[- ]citrate lyase|citrate lyase", "favorable"),
    ("M02", "Carbon import and acetyl-CoA supply", "Malic enzyme", r"malic enzyme|malate dehydrogenase.*decarboxylating", "favorable"),
    ("M02", "Carbon import and acetyl-CoA supply", "PPP", r"glucose-6-phosphate dehydrogenase|6-phosphogluconate dehydrogenase|transketolase|transaldolase", "favorable"),
    ("M02", "Carbon import and acetyl-CoA supply", "PPP description", r"oxidative pentose-phosphate pathway|oxidative pentose phosphate pathway", "favorable"),
    ("M02", "Carbon import and acetyl-CoA supply", "Plastid carbon transport", r"phosphoenolpyruvate.*translocator|triose phosphate.*translocator|pyruvate transporter|bile acid sodium symporter 2|\bbass2\b|plastidic nucleotide transporter", "favorable"),
    ("M02", "Carbon import and acetyl-CoA supply", "PEP carboxylation", r"phosphoenolpyruvate carboxylase|phosphoenolpyruvate carboxykinase", "favorable"),
    ("M02", "Carbon import and acetyl-CoA supply", "Starch supply", r"adp[- ]glucose pyrophosphorylase|starch synthase|starch branching enzyme|starch debranching enzyme|alpha-amylase|beta-amylase|synthesis of starch.*adp", "context"),
    ("M02", "Carbon import and acetyl-CoA supply", "Fruit carbon fixation", r"ribulose.*bisphosphate carboxylase|\brubisco\b|carbon dioxide fixation", "context"),
    ("M02", "Carbon import and acetyl-CoA supply", "Broad sugar transport", r"sugar \(and other\) transporter|nucleotide-diphospho-sugar transferase|udp-glucose 6-dehydrogenase", "context"),

    # De novo fatty-acid synthesis.
    ("M03", "De novo fatty-acid synthesis", "ACCase", r"acetyl[- ]coenzyme a carboxylase|acetyl-coa carboxylase|biotin carboxyl carrier protein|\bbiotin carboxylase\b|carboxyl transferase.*acetyl", "favorable"),
    ("M03", "De novo fatty-acid synthesis", "MCAT", r"malonyl[- ]coa.*acyl carrier protein transacylase|malonyl[- ]acp transacylase|\bmcat\b", "favorable"),
    ("M03", "De novo fatty-acid synthesis", "ACP", r"\bacyl carrier protein\b|carrier of the growing fatty acid chain", "favorable"),
    ("M03", "De novo fatty-acid synthesis", "KAS", r"3-oxoacyl[- ]acyl carrier protein synthase|beta-ketoacyl-acp synthase|\bkas\s*(i|ii|iii|1|2|3)\b", "favorable"),
    ("M03", "De novo fatty-acid synthesis", "KAR", r"3-oxoacyl[- ]acyl carrier protein reductase|beta-ketoacyl-acp reductase|\bkar\b", "favorable"),
    ("M03", "De novo fatty-acid synthesis", "HAD", r"3-hydroxyacyl[- ]acyl carrier protein dehydratase|hydroxyacyl-acp dehydratase|\bhad\b", "favorable"),
    ("M03", "De novo fatty-acid synthesis", "ENR", r"enoyl[- ]acyl carrier protein reductase|enoyl-acp reductase|\benr\b", "favorable"),
    ("M03", "De novo fatty-acid synthesis", "FATA/FATB", r"acyl[- ]acyl carrier protein thioesterase|\bfata\b|\bfatb\b", "composition"),

    # Chain length and unsaturation.
    ("M04", "Fatty-acid elongation and composition", "SAD/FAB2", r"stearoyl[- ]acyl carrier protein desaturase|stearoyl-acp desaturase|acyl[- ]+acyl-carrier-protein desaturase|introduction of a cis double bond.*acyl|\bfab2\b", "composition"),
    ("M04", "Fatty-acid elongation and composition", "FAD2", r"omega-6 fatty acid desaturase|\bfad2\b", "composition"),
    ("M04", "Fatty-acid elongation and composition", "FAD3/6/7/8", r"omega-3 fatty acid desaturase|\bfad3\b|\bfad6\b|\bfad7\b|\bfad8\b", "composition"),
    ("M04", "Fatty-acid elongation and composition", "KCS", r"3-ketoacyl-coa synthase|beta-ketoacyl-coa synthase|very-long-chain.*condensing|\bkcs\b", "composition"),
    ("M04", "Fatty-acid elongation and composition", "KCR/HCD/ECR", r"3-ketoacyl-coa reductase|3-hydroxyacyl-coa dehydratase|trans-2,3-enoyl-coa reductase|\bkcr\b|\becr\b", "composition"),
    ("M04", "Fatty-acid elongation and composition", "Other desaturase", r"fatty acid desaturase|acyl-lipid desaturase", "composition"),
    ("M04", "Fatty-acid elongation and composition", "Sterol desaturase", r"sterol desaturase", "context"),
    ("M04", "Fatty-acid elongation and composition", "Elongation description", r"long-chain fatty acids elongation cycle|very long-chain fatty acids.*per cycle", "composition"),

    # Acyl activation and TAG assembly.
    ("M05", "Acyl activation and TAG assembly", "LACS", r"long-chain acyl-coa synthetase|long-chain-fatty-acid--coa ligase|\blacs\b", "favorable"),
    ("M05", "Acyl activation and TAG assembly", "ACBP", r"acyl-coa-binding protein|acyl-coenzyme a binding", "favorable"),
    ("M05", "Acyl activation and TAG assembly", "GPAT", r"glycerol-3-phosphate acyltransferase|\bgpat\b", "favorable"),
    ("M05", "Acyl activation and TAG assembly", "LPAT", r"lysophosphatidic acid acyltransferase|1-acylglycerol-3-phosphate acyltransferase|\blpat\b", "favorable"),
    ("M05", "Acyl activation and TAG assembly", "PAP/PAH", r"phosphatidate phosphatase|phosphatidic acid phosphatase|lipid phosphate phosphatase|\bpah[12]?\b", "favorable"),
    ("M05", "Acyl activation and TAG assembly", "DGAT", r"diacylglycerol acyltransferase|\bdgat\b", "favorable"),
    ("M05", "Acyl activation and TAG assembly", "PDAT", r"phospholipid:diacylglycerol acyltransferase|\bpdat\b", "favorable"),
    ("M05", "Acyl activation and TAG assembly", "LPCAT", r"lysophosphatidylcholine acyltransferase|\blpcat\b", "composition"),
    ("M05", "Acyl activation and TAG assembly", "PDCT/CPT", r"phosphatidylcholine:diacylglycerol cholinephosphotransferase|phosphatidylcholine diacylglycerol cholinephosphotransferase|\bpdct\b|cholinephosphotransferase", "composition"),
    ("M05", "Acyl activation and TAG assembly", "Glycerol supply", r"glycerol kinase|glycerol-3-phosphate dehydrogenase", "favorable"),
    ("M05", "Acyl activation and TAG assembly", "Broad acyltransferase", r"membrane-bound acyltransferase|\bo-acyltransferase\b|^acyltransferase$|\| acyltransferase \|", "context"),
    ("M05", "Acyl activation and TAG assembly", "Phospholipid synthesis", r"phosphatidylserine decarboxylase|cdp-diacylglycerol|phosphatidylethanolamine|phosphatidylinositol.*phosphatidyltransferase", "context"),
    ("M05", "Acyl activation and TAG assembly", "Galactolipid synthesis", r"monogalactosyldiacylglycerol synthase|digalactosyldiacylglycerol synthase", "context"),

    # Oil bodies and lipid trafficking.
    ("M06", "Oil-body biogenesis and lipid storage", "Oleosin", r"oleosin|oil body protein", "favorable"),
    ("M06", "Oil-body biogenesis and lipid storage", "Caleosin/steroleosin", r"caleosin|steroleosin", "favorable"),
    ("M06", "Oil-body biogenesis and lipid storage", "SEIPIN", r"seipin", "favorable"),
    ("M06", "Oil-body biogenesis and lipid storage", "LDAP/LDIP", r"lipid droplet-associated protein|lipid droplet-interacting protein|\bldap\b|\bldip\b", "favorable"),
    ("M06", "Oil-body biogenesis and lipid storage", "Lipid transfer", r"lipid transfer protein|phospholipid transfer protein|lipid-transfer protein", "context"),
    ("M06", "Oil-body biogenesis and lipid storage", "Lipid transporter", r"abc.*lipid transporter|lipid transporter|fatty acid transporter", "context"),

    # Mesocarp, sink strength, and cell wall.
    ("M07", "Mesocarp growth, sink strength and cell wall", "Cellulose synthase", r"cellulose synthase|cellulose-synthase-like", "context"),
    ("M07", "Mesocarp growth, sink strength and cell wall", "Expansin", r"expansin", "favorable"),
    ("M07", "Mesocarp growth, sink strength and cell wall", "XTH", r"xyloglucan endotransglucosylase|xyloglucan endotransglycosylase|\bxth\b", "favorable"),
    ("M07", "Mesocarp growth, sink strength and cell wall", "XEH/XET description", r"xyloglucan endohydrolysis|endotransglycosylation.*xyloglucan|cleaves and religates xyloglucan", "favorable"),
    ("M07", "Mesocarp growth, sink strength and cell wall", "PME/PMEI", r"pectin methylesterase|pectinesterase|pectin methylesterase inhibitor", "context"),
    ("M07", "Mesocarp growth, sink strength and cell wall", "Polygalacturonase", r"polygalacturonase|pectate lyase", "context"),
    ("M07", "Mesocarp growth, sink strength and cell wall", "Pectin acetylesterase", r"homogalacturonan.*pectin|acetyl esters.*pectin|degree of acetylation of pectin", "context"),
    ("M07", "Mesocarp growth, sink strength and cell wall", "Glycosidase", r"beta-galactosidase|alpha-arabinofuranosidase|endoglucanase", "context"),
    ("M07", "Mesocarp growth, sink strength and cell wall", "Arabinogalactan", r"arabinogalactan protein", "context"),
    ("M07", "Mesocarp growth, sink strength and cell wall", "Lignin pathway", r"phenylalanine ammonia-lyase|cinnamate 4-hydroxylase|4-coumarate.*ligase|cinnamoyl-coa reductase|cinnamyl alcohol dehydrogenase|caffeic acid.*methyltransferase|caffeoyl-coa.*methyltransferase|hydroxycinnamoyl.*transferase|laccase", "context"),
    ("M07", "Mesocarp growth, sink strength and cell wall", "Xylan/cell-wall assembly", r"glucuronoxylan|plant-type cell wall assembly|cell wall construction|lignin degradation", "context"),
    ("M07", "Mesocarp growth, sink strength and cell wall", "Seedstick", r"seedstick", "context"),
    ("M07", "Mesocarp growth, sink strength and cell wall", "Aquaporin", r"aquaporin|plasma membrane intrinsic protein|tonoplast intrinsic protein", "favorable"),

    # Antioxidant and quality retention.
    ("M08", "Antioxidant and postharvest stability", "SOD", r"superoxide dismutase", "stability"),
    ("M08", "Antioxidant and postharvest stability", "Catalase", r"\bcatalase\b", "stability"),
    ("M08", "Antioxidant and postharvest stability", "APX/GPX", r"ascorbate peroxidase|glutathione peroxidase", "stability"),
    ("M08", "Antioxidant and postharvest stability", "Glutathione cycle", r"glutathione reductase|dehydroascorbate reductase|monodehydroascorbate reductase|glutathione synthetase|gamma-glutamylcysteine synthetase", "stability"),
    ("M08", "Antioxidant and postharvest stability", "GST", r"glutathione s-transferase", "stability"),
    ("M08", "Antioxidant and postharvest stability", "Thioredoxin/peroxiredoxin", r"thioredoxin|peroxiredoxin|glutaredoxin", "stability"),
    ("M08", "Antioxidant and postharvest stability", "Tocopherol/VTE", r"tocopherol|homogentisate phytyltransferase|tocopherol cyclase|\bvte[1-6]\b", "stability"),
    ("M08", "Antioxidant and postharvest stability", "Carotenoid", r"phytoene synthase|phytoene desaturase|zeta-carotene desaturase|lycopene.*cyclase|carotenoid cleavage", "stability"),
    ("M08", "Antioxidant and postharvest stability", "Class III peroxidase", r"class iii peroxidase|peroxidase family|plant peroxidase", "context"),
    ("M08", "Antioxidant and postharvest stability", "Other peroxidase", r"\bperoxidase\b|ahpc/tsa antioxidant|thiol-specific.*hydroperoxide|protect cells.*hydrogen peroxide", "context"),
    ("M08", "Antioxidant and postharvest stability", "Polyphenol oxidase", r"polyphenol oxidase", "context"),
    ("M08", "Antioxidant and postharvest stability", "Ferritin/iron redox", r"\bferritin\b|stores iron.*non-toxic", "context"),

    # Lipolysis and oxidative-rancidity risks.
    ("M09", "Lipolysis and oxidative-rancidity risk", "LOX", r"lipoxygenase|\blox[1-9]?\b", "risk"),
    ("M09", "Lipolysis and oxidative-rancidity risk", "TAG lipase", r"triacylglycerol lipase|diacylglycerol lipase|monoacylglycerol lipase|sugar-dependent 1|\bsdp1\b", "risk"),
    ("M09", "Lipolysis and oxidative-rancidity risk", "Patatin/lipase", r"patatin-like phospholipase|patatin.*lipase|lipase.*class 3|class 3.*lipase|gdxg lipase|gdsl esterase lipase|ab hydrolase.*lipase family|alpha/beta hydrolase family", "risk"),
    ("M09", "Lipolysis and oxidative-rancidity risk", "Phospholipase", r"\bphospholipase\b|hydrolyzes glycerol-phospholipids", "risk"),
    ("M09", "Lipolysis and oxidative-rancidity risk", "Oxylipin cleavage", r"hydroperoxide lyase|allene oxide synthase|alpha-dioxygenase", "risk"),
    ("M09", "Lipolysis and oxidative-rancidity risk", "Beta oxidation", r"acyl-coa oxidase|multifunctional protein.*beta-oxidation|3-ketoacyl-coa thiolase|peroxisomal.*beta-oxidation", "risk"),
    ("M09", "Lipolysis and oxidative-rancidity risk", "Acyl-CoA dehydrogenase", r"acyl-coa dehydrogenase", "context"),

    # Targeted upstream transport and regulation not captured above.
    ("M10", "Upstream lipid regulation and trafficking", "Acyl editing/regulation", r"phospholipid:diacylglycerol|fatty acid export|fatty acid amide hydrolase|acyl-coa oxidase regulator", "context"),
    ("M10", "Upstream lipid regulation and trafficking", "Plastid lipid export", r"trigalactosyldiacylglycerol|tgd[1-5]|fatty acid export protein|fax[1-9]", "context"),
    ("M10", "Upstream lipid regulation and trafficking", "Phospholipid signalling", r"phosphatidylinositol.*kinase|diacylglycerol kinase|phosphatidic acid binding", "context"),

    # Curated prior topic and trait-module fallbacks. These preserve earlier evidence without
    # elevating it automatically; molecular evidence and resolved D/P/M are still mandatory.
    ("M02", "Carbon import and acetyl-CoA supply", "Prior carbon-allocation topic", r"sugar_lipid_carbon_allocation", "context"),
    ("M03", "De novo fatty-acid synthesis", "Prior high-oil topic", r"high_oil_synthesis_storage|de-novo / saturated fa", "context"),
    ("M04", "Fatty-acid elongation and composition", "Prior composition topic", r"oleic_and_unsaturation|unsaturated fa", "composition"),
    ("M05", "Acyl activation and TAG assembly", "Prior oil-storage module", r"oil biosynthesis & storage|tag assembly & oil body", "context"),
    ("M07", "Mesocarp growth, sink strength and cell wall", "Prior mesocarp/cell-wall topic", r"fruit_mesocarp_thickening|shell / cell wall / lignin", "context"),
    ("M08", "Antioxidant and postharvest stability", "Prior redox topic", r"oxidative_rancidity_redox|lipid oxidation / antioxidant|phenylpropanoid_defence", "context"),
    ("M09", "Lipolysis and oxidative-rancidity risk", "Prior hydrolytic-rancidity topic", r"hydrolytic_rancidity", "risk"),
]


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    for name in ("catalog", "mechanism", "protein", "heterosis", "ase", "bins", "gtf"):
        p.add_argument(f"--{name}", required=True, type=Path)
    p.add_argument("--output-dir", required=True, type=Path)
    p.add_argument("--checkpoint-dir", required=True, type=Path)
    return p.parse_args()


def read_tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False, low_memory=False)


def n(value, default=0.0) -> float:
    try:
        x = float(value)
        return default if math.isnan(x) else x
    except (TypeError, ValueError):
        return default


def yes(value) -> bool:
    return str(value).strip().lower() in {"true", "yes", "1", "t", "y"}


def first_text(values) -> str:
    vals = [str(x).strip() for x in values if str(x).strip() not in {"", "NA", "nan", "None"}]
    return max(vals, key=len) if vals else ""


def parse_gtf(path: Path) -> pd.DataFrame:
    rows = []
    pat = re.compile(r'gene_id "([^"]+)"')
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) != 9 or parts[2] != "transcript":
                continue
            m = pat.search(parts[8])
            if not m:
                continue
            chrom = re.match(r"chr(\d+)", parts[0])
            rows.append((m.group(1), parts[0], f"chr{int(chrom.group(1)):02d}" if chrom else parts[0], int(parts[3]) - 1, int(parts[4])))
    out = pd.DataFrame(rows, columns=["gene_id", "reference_seqid", "Chromosome", "Start0", "End0"])
    return out.drop_duplicates("gene_id")


def normalize_genotype(a: str, b: str) -> str:
    if a not in ALLOWED or b not in ALLOWED:
        return ""
    codes = sorted((ALLOWED[a], ALLOWED[b]), key=lambda x: "DPM".index(x))
    return "/".join(codes)


def aggregate_ase(df: pd.DataFrame) -> pd.DataFrame:
    def reference_gene(pair: str) -> str:
        xs = pair.split("|")
        ref = next((x for x in xs if re.search(r"\.chr\d+B\.", x)), xs[0] if xs else "")
        return ref.replace("evm.model.", "evm.TU.")

    d = df.copy()
    d["gene_id"] = d["Allele_pair"].map(reference_gene)
    d["ASE_sample_count_num"] = pd.to_numeric(d["ASE_sample_count"], errors="coerce").fillna(0)
    records = []
    for gid, sub in d.groupby("gene_id", sort=False):
        rec = {"gene_id": gid}
        for individual in ("FL", "TN"):
            x = sub[sub["Individual"] == individual]
            rec[f"ase_{individual}_allele_rows"] = len(x)
            rec[f"ase_{individual}_max_sample_count"] = int(x["ASE_sample_count_num"].max()) if len(x) else 0
            rec[f"ase_{individual}_positive_alleles"] = int((x["ASE_in_at_least_one_sample"] == "YES").sum()) if len(x) else 0
            rec[f"ase_{individual}_ancestries"] = ";".join(sorted(set(x["Ancestry_call"]) - {"", "Unknown"})) if len(x) else ""
        records.append(rec)
    return pd.DataFrame(records)


def best_by_score(df: pd.DataFrame, key: str, score: str, ascending=False) -> pd.DataFrame:
    d = df.copy()
    d["_score"] = pd.to_numeric(d[score], errors="coerce")
    d = d.sort_values("_score", ascending=ascending, na_position="last")
    return d.drop_duplicates(key).drop(columns="_score")


def join_selected(base: pd.DataFrame, other: pd.DataFrame, key: str, cols: list[str], prefix: str) -> pd.DataFrame:
    existing = [x for x in cols if x in other.columns]
    d = other[[key] + existing].copy()
    d = d.rename(columns={key: "gene_id", **{c: f"{prefix}{c}" for c in existing}})
    return base.merge(d, on="gene_id", how="left", validate="one_to_one")


def main() -> None:
    a = parse_args()
    if a.output_dir.exists():
        raise SystemExit(f"[ERROR] output exists: {a.output_dir}")
    stage = a.output_dir.with_name(f"{a.output_dir.name}.building.{os.getpid()}")
    stage.mkdir(parents=True)
    a.checkpoint_dir.mkdir(parents=True, exist_ok=True)

    catalog = read_tsv(a.catalog)
    if len(catalog) != 25480 or catalog["gene_id"].duplicated().any():
        raise SystemExit(f"[ERROR] primary catalog expected 25480 unique genes; rows={len(catalog)} dup={catalog['gene_id'].duplicated().sum()}")

    coords = parse_gtf(a.gtf)
    mechanism = best_by_score(read_tsv(a.mechanism), "common_gene_id", "mechanism_rank", ascending=True)
    protein_raw = read_tsv(a.protein)
    protein = best_by_score(protein_raw, "reference_gene_alias", "proteomics_support_score", ascending=False)
    heterosis = read_tsv(a.heterosis).drop_duplicates("gene_africa")
    ase = aggregate_ase(read_tsv(a.ase)).drop_duplicates("gene_id")
    bins = read_tsv(a.bins)

    evidence = catalog.merge(coords, on="gene_id", how="left", validate="one_to_one")
    mechanism_cols = [
        "robust_link_n", "robust_global_link_n", "robust_compound_n", "named_robust_compound_n",
        "max_abs_partial_r", "min_permutation_global_fdr", "linked_axes", "linked_compounds",
        "cis_trans_robust_ASE_stage_n", "expression_above_better_parent_stage_n",
        "median_expression_BPH_log2", "mechanism_priority_score", "evidence_tier", "mechanism_rank",
    ]
    evidence = join_selected(evidence, mechanism, "common_gene_id", mechanism_cols, "mech_")
    protein_cols = [
        "allele_unit_id", "product", "corrected_candidate_tier", "corrected_comprehensive_score", "topic_memberships",
        "current_protein_sample_completeness", "current_protein_quantified_sample_n", "current_stagewise_min_FDR",
        "current_significant_TN_FL_stage_n", "current_strong_abslog2FC1_TN_FL_stage_n", "current_average_TN_FL_log2FC",
        "current_average_TN_FL_FDR", "current_genotype_stage_interaction_FDR", "current_WGCNA_max_abs_kME",
        "best_cross_omics_module_permutation_FDR", "protein_heterosis_stage_n", "above_better_parent_stage_n",
        "protein_median_log2_better_parent_effect", "actual_haplotype_specific_protein_claim_allowed",
        "proteomics_evidence_scope", "proteomics_support_score", "proteomics_followup_rank",
    ]
    evidence = join_selected(evidence, protein, "reference_gene_alias", protein_cols, "prot_")
    heterosis_cols = [
        "Stages_tested", "Above_better_parent_stages", "Median_log2_MPV_effect", "Max_log2_HPV_effect",
        "Median_log2_HPV_effect_among_ABPH_stages", "Recurrent_ABPH_ge3_stages", "family", "preferred_name",
        "description", "trait_module", "TN_homolog1_gene", "TN_homolog2_gene", "TN_homolog1_start0",
        "TN_homolog2_start0", "TN_homolog1_ancestry", "TN_homolog2_ancestry", "Diploid_state",
        "Ancestry_genotype", "Inference_level",
    ]
    evidence = join_selected(evidence, heterosis, "gene_africa", heterosis_cols, "het_")
    evidence = evidence.merge(ase, on="gene_id", how="left", validate="one_to_one")

    # Assign the FL diploid ancestry at the homologous chromosome fraction. chrB is the reference (FL_HapB).
    fl_bins = bins[bins["Individual"] == "FL"].copy()
    for c in ("Fraction_start", "Fraction_end"):
        fl_bins[c] = pd.to_numeric(fl_bins[c], errors="coerce")
    chrom_lengths = coords.groupby("Chromosome")["End0"].max().to_dict()
    fl_lookup = {(r.Chromosome, int(r.Bin_index)): r for r in fl_bins.itertuples()}
    fl_hap1, fl_hap2, fl_geno, frac = [], [], [], []
    for r in evidence.itertuples():
        if not r.Chromosome or pd.isna(r.Start0) or pd.isna(r.End0) or r.Chromosome not in chrom_lengths:
            fl_hap1.append(""); fl_hap2.append(""); fl_geno.append(""); frac.append(np.nan); continue
        f = min(max(((float(r.Start0) + float(r.End0)) / 2) / chrom_lengths[r.Chromosome], 0.0), 0.999999)
        b = fl_lookup.get((r.Chromosome, int(f * 100)))
        frac.append(f)
        if b is None:
            fl_hap1.append(""); fl_hap2.append(""); fl_geno.append("")
        else:
            fl_hap1.append(b.Hap1_ancestry); fl_hap2.append(b.Hap2_ancestry)
            fl_geno.append(normalize_genotype(b.Hap1_ancestry, b.Hap2_ancestry))
    evidence["chromosome_fraction"] = frac
    evidence["FL_HapA_ancestry"] = fl_hap1
    evidence["FL_HapB_ancestry"] = fl_hap2
    evidence["FL_resolved_genotype"] = fl_geno
    evidence["TN_resolved_genotype"] = [normalize_genotype(x, y) for x, y in zip(
        evidence["het_TN_homolog1_ancestry"].fillna(""), evidence["het_TN_homolog2_ancestry"].fillna("")
    )]

    # Versioned dictionary.
    dictionary = pd.DataFrame(RULES, columns=["Module_ID", "Module", "Family", "Regex", "Direction_prior"])
    dictionary.to_csv(stage / "functional_dictionary_v1.tsv", sep="\t", index=False)
    compiled = [(mid, module, family, re.compile(pattern, re.I), prior) for mid, module, family, pattern, prior in RULES]

    outputs = defaultdict(list)
    for _, row in evidence.iterrows():
        texts = [
            row.get("product", ""), row.get("prot_product", ""), row.get("trait_family", ""),
            row.get("trait_preferred_name", ""), row.get("trait_description", ""), row.get("trait_module", ""),
            row.get("BGC_descriptions", ""), row.get("BGC_preferred_names", ""), row.get("core_selection_reasons", ""),
            row.get("prot_topic_memberships", ""), row.get("het_family", ""), row.get("het_preferred_name", ""),
            row.get("het_description", ""), row.get("het_trait_module", ""),
        ]
        annotation = " | ".join(str(x) for x in texts if str(x).strip() not in {"", "NA", "nan"})
        matches = [(mid, module, family, prior) for mid, module, family, pattern, prior in compiled if pattern.search(annotation)]
        outputs["functional_annotation"].append(annotation)
        outputs["candidate_flag"].append("YES" if matches else "NO")
        outputs["module_ids"].append(";".join(dict.fromkeys(x[0] for x in matches)))
        outputs["modules"].append(";".join(dict.fromkeys(x[1] for x in matches)))
        outputs["families"].append(";".join(dict.fromkeys(x[2] for x in matches)))
        outputs["direction_priors"].append(";".join(dict.fromkeys(x[3] for x in matches)))
        outputs["primary_module"].append(matches[0][1] if matches else "")
        outputs["primary_family"].append(matches[0][2] if matches else "")

    for key, vals in outputs.items():
        evidence[key] = vals

    score, class_count, mol_count = [], [], []
    rna_s, prot_s, mech_s, ase_s, het_s = [], [], [], [], []
    actions, targets, secondary, tiers, reasons, limitations = [], [], [], [], [], []
    for _, row in evidence.iterrows():
        candidate = row["candidate_flag"] == "YES"
        rna = any([
            yes(row.get("expression_interaction_support", "")), yes(row.get("DE_interaction_development", "")),
            yes(row.get("DE_interaction_postharvest", "")), n(row.get("interaction_padj_development", 1), 1) <= 0.05,
            n(row.get("interaction_padj_postharvest", 1), 1) <= 0.05,
        ])
        prot = any([
            n(row.get("prot_current_significant_TN_FL_stage_n", 0)) >= 1,
            n(row.get("prot_current_strong_abslog2FC1_TN_FL_stage_n", 0)) >= 1,
            n(row.get("prot_best_cross_omics_module_permutation_FDR", 1), 1) <= 0.05,
        ])
        mech = any([
            n(row.get("mech_robust_link_n", 0)) >= 1,
            yes(row.get("strict_compound_module_support", "")),
            n(row.get("strict_RNA_compound_module_link_n", 0)) >= 1,
        ])
        ase_support = any([
            yes(row.get("robust_FL_ASE_support", "")), yes(row.get("robust_TN_ASE_support", "")),
            n(row.get("ase_FL_max_sample_count", 0)) >= 1, n(row.get("ase_TN_max_sample_count", 0)) >= 1,
        ])
        heter = yes(row.get("het_Recurrent_ABPH_ge3_stages", "")) or n(row.get("het_Above_better_parent_stages", 0)) >= 3
        tn = str(row.get("TN_resolved_genotype", "") or "")
        fl = str(row.get("FL_resolved_genotype", "") or "")
        classes = int(candidate) + int(rna) + int(prot) + int(mech) + int(ase_support) + int(heter) + int(bool(tn)) + int(bool(fl))
        molecular = sum((rna, prot, mech, ase_support, heter))
        sc = 2 * int(candidate) + int(rna) + 2 * int(prot) + 2 * int(mech) + int(ase_support) + int(heter) + 2 * int(bool(tn)) + int(bool(fl))

        modules = set(str(row.get("module_ids", "")).split(";"))
        avg = n(row.get("prot_current_average_TN_FL_log2FC", 0), 0)
        action, target, alt = "NOT_TRAIT_RELEVANT", "", ""
        why, limit = [], []
        if candidate:
            why.append("direct_functional_recall")
            if rna: why.append("replicated_RNA_timecourse")
            if prot: why.append("replicated_protein_timecourse")
            if mech: why.append("adjusted_multiomic_mechanism")
            if ase_support: why.append("ASE_or_cis_support")
            if heter: why.append("exploratory_recurrent_parent_F1_ABPH")
            if "M09" in modules:
                if fl and (avg > 0.25 or heter):
                    action, target = "EXCLUDE_TN_RISK", fl
                    why.append("TN_elevated_or_ABPH_rancidity_risk")
                else:
                    action = "VALIDATE_ONLY"
                    limit.append("risk_direction_not_resolved_for_action")
            elif "M01" in modules:
                if tn and (heter or avg > 0.25 or mech):
                    action = "TIMING_SCREEN"
                    if tn in {"D/D", "P/P"} and "M" in fl:
                        target = f"{tn[0]}/M"
                        alt = tn
                    else:
                        target = tn
                    why.append("candidate_for_earlier_or_longer_lipid_window")
                elif fl and avg < -0.25:
                    action, target = "RETAIN_FL", fl
                else:
                    action = "VALIDATE_ONLY"
                    limit.append("timing_effect_requires_stage_resolved_validation")
            elif "M08" in modules or "M04" in modules:
                if fl and (avg <= 0.0 or "FAD2" in row["families"]):
                    action, target = "RETAIN_FL", fl
                    why.append("FL_quality_or_stability_background")
                elif tn and avg > 0.25 and "FAD2" not in row["families"]:
                    action, target = "INTRODUCE_OR_TUNE", tn
                else:
                    action = "VALIDATE_ONLY"
                    limit.append("composition_or_antioxidant_direction_requires_validation")
            else:
                if tn and (heter or avg > 0.25):
                    action, target = "INTRODUCE_OR_TUNE", tn
                    why.append("TN_expression_or_protein_advantage")
                    if tn == "D/D" and "M" in fl: alt = "D/M"
                    elif tn == "P/P" and "M" in fl: alt = "P/M"
                    elif tn == "D/P": alt = "D/D;P/P"
                elif fl and avg < -0.25:
                    action, target = "RETAIN_FL", fl
                    why.append("FL_protein_abundance_advantage")
                else:
                    action = "VALIDATE_ONLY"
                    limit.append("functional_match_without_directional_donor_evidence")

        actionable = action in {"INTRODUCE_OR_TUNE", "RETAIN_FL", "EXCLUDE_TN_RISK", "TIMING_SCREEN"} and target in {"D/D", "P/P", "M/M", "D/P", "D/M", "P/M"}
        if actionable and sc >= 7 and classes >= 3 and molecular >= 2:
            tier = "Tier_A"
        elif actionable and sc >= 5 and classes >= 2 and molecular >= 1:
            tier = "Tier_B"
        else:
            tier = "Tier_C"
            if candidate and not actionable:
                limit.append("no_complete_actionable_DPM_target")
        if candidate and not tn: limit.append("TN_primitive_diplotype_unresolved")
        if candidate and not fl: limit.append("FL_primitive_diplotype_unresolved")
        limit.append("ancestry_is_not_functional_causality") if candidate else None
        if heter: limit.append("parent_F1_n1_exploratory")

        rna_s.append(rna); prot_s.append(prot); mech_s.append(mech); ase_s.append(ase_support); het_s.append(heter)
        score.append(sc); class_count.append(classes); mol_count.append(molecular)
        actions.append(action); targets.append(target); secondary.append(alt); tiers.append(tier)
        reasons.append(";".join(dict.fromkeys(why))); limitations.append(";".join(dict.fromkeys(limit)))

    evidence["RNA_support"] = rna_s
    evidence["protein_support"] = prot_s
    evidence["mechanism_support"] = mech_s
    evidence["ASE_support"] = ase_s
    evidence["exploratory_parent_F1_ABPH_support"] = het_s
    evidence["molecular_support_class_n"] = mol_count
    evidence["independent_evidence_class_n_recomputed"] = class_count
    evidence["breeding_evidence_score"] = score
    evidence["breeding_action"] = actions
    evidence["target_DPM_genotype"] = targets
    evidence["secondary_genotype_test"] = secondary
    evidence["breeding_tier"] = tiers
    evidence["inclusion_reasons"] = reasons
    evidence["limitations_and_validation"] = limitations

    evidence.to_csv(stage / "complete_25480_gene_evidence_matrix.tsv", sep="\t", index=False, na_rep="")
    candidates = evidence[evidence["candidate_flag"] == "YES"].copy()
    candidates = candidates.sort_values(["breeding_tier", "breeding_evidence_score", "Chromosome", "Start0"], ascending=[True, False, True, True])
    candidates.to_csv(stage / "all_trait_relevant_candidates.tsv", sep="\t", index=False, na_rep="")
    actionable = candidates[candidates["breeding_tier"].isin(["Tier_A", "Tier_B"])].copy()
    actionable.to_csv(stage / "tierA_tierB_actionable_breeding_targets.tsv", sep="\t", index=False, na_rep="")
    risk = candidates[candidates["module_ids"].str.contains(r"(^|;)M09(;|$)", regex=True, na=False)].copy()
    risk.to_csv(stage / "risk_alleles_and_linkage_drag.tsv", sep="\t", index=False, na_rep="")
    audit = candidates[(candidates["breeding_tier"] == "Tier_C") | (candidates["limitations_and_validation"].str.contains("unresolved|conflict", case=False, na=False))].copy()
    audit.to_csv(stage / "unresolved_or_conflicting_candidates_audit.tsv", sep="\t", index=False, na_rep="")

    module_rows = []
    for mid, module, *_ in RULES:
        if any(x["Module_ID"] == mid for x in module_rows):
            continue
        sub = candidates[candidates["module_ids"].str.contains(fr"(^|;){mid}(;|$)", regex=True, na=False)]
        module_rows.append({
            "Module_ID": mid, "Module": module, "Candidate_n": len(sub),
            "Tier_A_n": int((sub["breeding_tier"] == "Tier_A").sum()),
            "Tier_B_n": int((sub["breeding_tier"] == "Tier_B").sum()),
            "Risk_action_n": int((sub["breeding_action"] == "EXCLUDE_TN_RISK").sum()),
        })
    pd.DataFrame(module_rows).to_csv(stage / "module_summary.tsv", sep="\t", index=False)
    chrom_summary = actionable.groupby(["Chromosome", "breeding_action", "target_DPM_genotype"], dropna=False).size().reset_index(name="Candidate_n")
    chrom_summary.to_csv(stage / "chromosome_action_genotype_summary.tsv", sep="\t", index=False)

    qc = {
        "Primary_catalog_rows": len(catalog),
        "Primary_catalog_unique_gene_ids": catalog["gene_id"].nunique(),
        "Coordinate_mapped_n": int(evidence["Start0"].notna().sum()),
        "Coordinate_unmapped_n": int(evidence["Start0"].isna().sum()),
        "Protein_joined_n": int(evidence["prot_allele_unit_id"].fillna("").ne("").sum()),
        "Mechanism_joined_n": int(evidence["mech_mechanism_rank"].fillna("").ne("").sum()),
        "Heterosis_joined_n": int(evidence["het_Stages_tested"].fillna("").ne("").sum()),
        "ASE_joined_n": int(evidence["ase_FL_allele_rows"].fillna(0).astype(float).add(evidence["ase_TN_allele_rows"].fillna(0).astype(float)).gt(0).sum()),
        "Trait_relevant_candidate_n": len(candidates),
        "Tier_A_n": int((candidates["breeding_tier"] == "Tier_A").sum()),
        "Tier_B_n": int((candidates["breeding_tier"] == "Tier_B").sum()),
        "Tier_C_n": int((candidates["breeding_tier"] == "Tier_C").sum()),
        "Actionable_complete_DPM_n": len(actionable),
        "Actionable_unresolved_target_n": int((~actionable["target_DPM_genotype"].isin(["D/D", "P/P", "M/M", "D/P", "D/M", "P/M"])).sum()),
    }
    if qc["Coordinate_unmapped_n"] / len(catalog) > 0.01:
        raise SystemExit(f"[ERROR] >1% coordinate-unmapped: {qc['Coordinate_unmapped_n']}")
    if qc["Actionable_unresolved_target_n"]:
        raise SystemExit("[ERROR] actionable candidates contain unresolved target genotypes")
    pd.DataFrame(qc.items(), columns=["Metric", "Value"]).to_csv(stage / "QC_summary.tsv", sep="\t", index=False)

    action_counts = Counter(actionable["breeding_action"])
    genotype_counts = Counter(actionable["target_DPM_genotype"])
    report = [
        "# 全基因D/P/M育种设计阶段性结果", "",
        f"- 全量筛查基因：{len(catalog):,}",
        f"- 广泛功能召回候选：{len(candidates):,}",
        f"- Tier A：{qc['Tier_A_n']:,}",
        f"- Tier B：{qc['Tier_B_n']:,}",
        f"- Tier C/补证：{qc['Tier_C_n']:,}",
        f"- 可执行且D/P/M完整解析：{len(actionable):,}", "",
        "## 育种动作", "",
    ]
    report += [f"- {k}: {v:,}" for k, v in sorted(action_counts.items())]
    report += ["", "## 目标基因型", ""]
    report += [f"- {k}: {v:,}" for k, v in sorted(genotype_counts.items())]
    report += [
        "", "## 解释边界", "",
        "- Tier A/B表示优先验证和育种设计价值，不代表已证明因果。",
        "- TN亲本–F1表达为每阶段n=1，仅作为探索性加分。",
        "- 祖源只用于定义等位来源，不能代替功能验证。",
        "- FAD2、细胞壁酶和时序调控因子需按发育阶段验证剂量与方向。",
    ]
    (stage / "RESULTS_AND_BREEDING_DESIGN.md").write_text("\n".join(report) + "\n", encoding="utf-8")

    # Atomic publication and stage checkpoints.
    stage.rename(a.output_dir)
    for p in ("P1_ID_CROSSWALK", "P2_FULL_GENE_EVIDENCE_MATRIX", "P3_FUNCTIONAL_RECALL", "P4_DIRECTION_AND_RISK_REVIEW",
              "P5_ANCESTRY_DIPLOTYPE_ASSIGNMENT", "P6_SCORING_AND_TIERS", "P7_TARGET_GENOTYPE_DESIGN", "P8_QC_AND_SENSITIVITY"):
        (a.checkpoint_dir / f"{p}.PASS").write_text("PASS\n", encoding="utf-8")
    print(f"COMPREHENSIVE_DPM_ANALYSIS_PASS candidates={len(candidates)} actionable={len(actionable)}")


if __name__ == "__main__":
    main()
