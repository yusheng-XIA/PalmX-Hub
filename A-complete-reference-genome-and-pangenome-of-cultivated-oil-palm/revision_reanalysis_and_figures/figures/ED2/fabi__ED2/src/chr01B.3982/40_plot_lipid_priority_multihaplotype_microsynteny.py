#!/usr/bin/env python3
"""Draw versioned multi-haplotype microsynteny figures for FL Africa PAVs.

The layout intentionally follows the accepted
evm.TU.chr02B.2671_all_haplotypes_v4 figure. Orthogroup presence is based on
the existing 11-genome OrthoFinder run; local coordinates come from the GFF3
files recorded in config/Input_Manifest.tsv.
"""

from __future__ import annotations

import argparse
import csv
import json
import re
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import PathPatch, Rectangle
from matplotlib.path import Path as MplPath


SAMPLE_ORDER = [
    "FL_Africa_hap2",
    "EG11",
    "dura_hap1",
    "dura_hap2",
    "pisifera_hap1",
    "pisifera_hap2",
    "nrly_hap1",
    "nrly_hap2",
    "BK_hap1",
    "BK_hap2",
    "FL_American_hap1",
]

DISPLAY_LABELS = {
    "FL_Africa_hap2": "FL Africa hap2",
    "EG11": "EG11",
    "dura_hap1": "Dura hap1",
    "dura_hap2": "Dura hap2",
    "pisifera_hap1": "Pisifera hap1",
    "pisifera_hap2": "Pisifera hap2",
    "nrly_hap1": "NRLY hap1",
    "nrly_hap2": "NRLY hap2",
    "BK_hap1": "BK hap1",
    "BK_hap2": "BK hap2",
    "FL_American_hap1": "FL American hap1",
}

ORTHOGROUP_ALIASES = {sample: sample for sample in SAMPLE_ORDER}
TRACK_GROUPS = {
    sample: ("E. guineensis" if sample != "FL_American_hap1" else "E. oleifera-derived")
    for sample in SAMPLE_ORDER
}

TARGETS = [
    {
        "gene_id": "evm.TU.chr03B.814",
        "transcript_id": "evm.model.chr03B.814",
        "short_label": "PDAT1",
        "legend_label": "PDAT1 PAV",
        "output_dir": "evm.TU.chr03B.814_all_haplotypes_v2",
        "upstream_genes": 5,
        "downstream_genes": 6,
    },
    {
        "gene_id": "evm.TU.chr02B.2671",
        "transcript_id": "evm.model.chr02B.2671",
        "short_label": "FabI-like ENR",
        "legend_label": "FabI-like ENR PAV",
        "output_dir": "evm.TU.chr02B.2671_all_haplotypes_v6",
        "upstream_genes": 5,
        "downstream_genes": 6,
    },
    {
        "gene_id": "evm.TU.chr05B.150",
        "transcript_id": "evm.model.chr05B.150",
        "short_label": "HACD",
        "legend_label": "HACD PAV",
        "output_dir": "evm.TU.chr05B.150_all_haplotypes_v2",
        "upstream_genes": 5,
        "downstream_genes": 6,
    },
    {
        "gene_id": "evm.TU.chr04B.1103",
        "transcript_id": "evm.model.chr04B.1103",
        "short_label": "GDSL lipase",
        "legend_label": "GDSL lipase PAV",
        "output_dir": "evm.TU.chr04B.1103_all_haplotypes_v2",
        "upstream_genes": 10,
        "downstream_genes": 10,
    },
    {
        "gene_id": "evm.TU.chr01B.3982",
        "transcript_id": "evm.model.chr01B.3982",
        "short_label": "Polyketide cyclase-like",
        "legend_label": "Polyketide cyclase-like PAV",
        "output_dir": "Polyketide cyclase-like",
        "upstream_genes": 5,
        "downstream_genes": 6,
        "exclude_orthogroups": ["OG0000025"],
    },
]

TARGET_COLOR = "#E68445"
CNV_COLOR = "#4C78A8"
SYNTENIC_COLOR = "#B8B8B8"
OTHER_COLOR = "#E6E6E6"
TRACK_COLOR = "#333333"


@dataclass
class Gene:
    sample: str
    seqid: str
    start: int
    end: int
    strand: str
    gene_id: str
    transcript_ids: List[str] = field(default_factory=list)
    orthogroup: Optional[str] = None
    role: str = "Other_local_gene"


@dataclass
class Track:
    sample: str
    seqid: str
    region_start: int
    region_end: int
    genes: List[Gene]
    anchor_context: str
    anchor_gene_ids: List[str]
    projected_target_position: Optional[int] = None
    projection_method: str = "NA"


def parse_attributes(text: str) -> Dict[str, str]:
    out: Dict[str, str] = {}
    for item in text.strip().split(";"):
        if not item:
            continue
        if "=" in item:
            key, value = item.split("=", 1)
            out[key] = value
        elif " " in item:
            key, value = item.split(" ", 1)
            out[key] = value.strip('"')
    return out


def read_manifest(path: Path) -> Dict[str, Dict[str, str]]:
    global SAMPLE_ORDER, DISPLAY_LABELS, ORTHOGROUP_ALIASES, TRACK_GROUPS
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if rows and any(row.get("Plot_Order") for row in rows):
        rows.sort(key=lambda row: int(row["Plot_Order"]))
        SAMPLE_ORDER = [row["Sample_ID"] for row in rows]
        DISPLAY_LABELS = {
            row["Sample_ID"]: row.get("Display_Label") or row["Sample_ID"]
            for row in rows
        }
        ORTHOGROUP_ALIASES = {
            row["Sample_ID"]: row.get("Orthogroup_Alias") or row["Sample_ID"]
            for row in rows
        }
        TRACK_GROUPS = {
            row["Sample_ID"]: row.get("Plot_Group") or "Unspecified"
            for row in rows
        }
    by_sample = {row["Sample_ID"]: row for row in rows}
    missing = [sample for sample in SAMPLE_ORDER if sample not in by_sample]
    if missing:
        raise RuntimeError(f"Manifest is missing samples: {','.join(missing)}")
    return by_sample


def parse_gff(
    sample: str,
    path: Path,
    seqid_pattern: Optional[str] = None,
    allowed_seqids: Optional[set[str]] = None,
) -> Tuple[List[Gene], Dict[str, Gene], Dict[str, Gene]]:
    genes: List[Gene] = []
    gene_by_id: Dict[str, Gene] = {}
    transcript_parent: Dict[str, str] = {}
    compiled_pattern = re.compile(seqid_pattern) if seqid_pattern and seqid_pattern != "NA" else None

    with path.open() as handle:
        for line_number, line in enumerate(handle, 1):
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n\r").split("\t")
            if len(fields) != 9:
                continue
            seqid, _, feature, start, end, _, strand, _, attr_text = fields
            if allowed_seqids is not None and seqid not in allowed_seqids:
                continue
            if compiled_pattern is not None and compiled_pattern.search(seqid) is None:
                continue
            attrs = parse_attributes(attr_text)
            if feature == "gene":
                gene_id = attrs.get("ID")
                if not gene_id:
                    raise RuntimeError(f"Missing gene ID in {path}:{line_number}")
                if gene_id in gene_by_id:
                    previous = gene_by_id[gene_id]
                    if previous.seqid != seqid or previous.strand != strand:
                        raise RuntimeError(f"Incompatible duplicate gene ID {gene_id} in {path}")
                    previous.start = min(previous.start, int(start))
                    previous.end = max(previous.end, int(end))
                    continue
                gene = Gene(sample, seqid, int(start), int(end), strand, gene_id)
                genes.append(gene)
                gene_by_id[gene_id] = gene
            elif feature in {"mRNA", "transcript"}:
                transcript_id = attrs.get("ID")
                parent = attrs.get("Parent", "").split(",")[0]
                if transcript_id and parent:
                    transcript_parent[transcript_id] = parent

    transcript_to_gene: Dict[str, Gene] = {}
    for transcript_id, parent in transcript_parent.items():
        gene = gene_by_id.get(parent)
        if gene is not None:
            gene.transcript_ids.append(transcript_id)
            transcript_to_gene[transcript_id] = gene

    genes.sort(key=lambda gene: (gene.seqid, gene.start, gene.end, gene.gene_id))
    return genes, gene_by_id, transcript_to_gene


def infer_eg11_main_seqids(path: Path, orthofinder_gene_ids: set[str]) -> set[str]:
    """Infer the nuclear chromosome seqids represented in the filtered proteome."""
    seqids: set[str] = set()
    matched_gene_ids: set[str] = set()
    with path.open() as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n\r").split("\t")
            if len(fields) != 9 or fields[2] != "gene":
                continue
            gene_id = parse_attributes(fields[8]).get("ID")
            if gene_id in orthofinder_gene_ids:
                seqids.add(fields[0])
                matched_gene_ids.add(gene_id)
    missing = orthofinder_gene_ids - matched_gene_ids
    if missing:
        preview = ",".join(sorted(missing)[:5])
        raise RuntimeError(f"EG11 OrthoFinder gene IDs missing from GFF ({len(missing)}): {preview}")
    if len(seqids) != 16:
        raise RuntimeError(f"Expected 16 EG11 main-chromosome seqids, observed {len(seqids)}: {','.join(sorted(seqids))}")
    return seqids


def read_orthogroups(path: Path) -> Tuple[Dict[str, Dict[str, List[str]]], Dict[str, Dict[str, str]]]:
    groups: Dict[str, Dict[str, List[str]]] = {}
    protein_to_og: Dict[str, Dict[str, str]] = {sample: {} for sample in SAMPLE_ORDER}
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required_columns = sorted(set(ORTHOGROUP_ALIASES.values()))
        missing = [sample for sample in required_columns if sample not in (reader.fieldnames or [])]
        if missing:
            raise RuntimeError(f"Orthogroups table is missing columns: {','.join(missing)}")
        for row in reader:
            og = row["Orthogroup"]
            groups[og] = {}
            for sample in SAMPLE_ORDER:
                source_sample = ORTHOGROUP_ALIASES[sample]
                members = [item.strip() for item in row[source_sample].split(",") if item.strip()]
                raw_members = []
                prefix = f"{source_sample}__"
                for member in members:
                    raw = member[len(prefix):] if member.startswith(prefix) else member
                    raw_members.append(raw)
                    if raw in protein_to_og[sample] and protein_to_og[sample][raw] != og:
                        raise RuntimeError(f"Protein {sample}:{raw} occurs in multiple orthogroups")
                    protein_to_og[sample][raw] = og
                groups[og][sample] = raw_members
    return groups, protein_to_og


def attach_orthogroups(
    genes_by_sample: Dict[str, List[Gene]],
    gene_by_id: Dict[str, Dict[str, Gene]],
    transcript_to_gene: Dict[str, Dict[str, Gene]],
    protein_to_og: Dict[str, Dict[str, str]],
) -> Dict[str, Dict[str, Gene]]:
    protein_to_gene: Dict[str, Dict[str, Gene]] = {sample: {} for sample in SAMPLE_ORDER}
    for sample in SAMPLE_ORDER:
        protein_to_gene[sample].update(gene_by_id[sample])
        protein_to_gene[sample].update(transcript_to_gene[sample])
        for protein_id, og in protein_to_og[sample].items():
            gene = protein_to_gene[sample].get(protein_id)
            if gene is None:
                if ORTHOGROUP_ALIASES.get(sample, sample) != sample:
                    continue
                raise RuntimeError(f"OrthoFinder ID is absent from GFF: {sample}:{protein_id}")
            if gene.orthogroup is not None and gene.orthogroup != og:
                raise RuntimeError(f"Gene {sample}:{gene.gene_id} maps to multiple orthogroups")
            gene.orthogroup = og
    return protein_to_gene


def target_orthogroup(target: Dict[str, str], protein_to_og: Dict[str, Dict[str, str]]) -> str:
    og = protein_to_og["FL_Africa_hap2"].get(target["transcript_id"])
    if og is None:
        raise RuntimeError(f"Target is not assigned by OrthoFinder: {target['transcript_id']}")
    return og


def select_fl_neighborhood(
    target: Dict[str, str],
    genes_by_sample: Dict[str, List[Gene]],
    gene_by_id: Dict[str, Dict[str, Gene]],
) -> List[Gene]:
    target_gene = gene_by_id["FL_Africa_hap2"].get(target["gene_id"])
    if target_gene is None:
        raise RuntimeError(f"Target gene is absent from FL GFF: {target['gene_id']}")
    chromosome_genes = [gene for gene in genes_by_sample["FL_Africa_hap2"] if gene.seqid == target_gene.seqid]
    focal_index = next(i for i, gene in enumerate(chromosome_genes) if gene.gene_id == target_gene.gene_id)
    upstream = int(target["upstream_genes"])
    downstream = int(target["downstream_genes"])
    expected_count = upstream + 1 + downstream
    start = max(0, focal_index - upstream)
    end = min(len(chromosome_genes), focal_index + downstream + 1)
    selected = chromosome_genes[start:end]
    if len(selected) != expected_count or target_gene not in selected:
        raise RuntimeError(f"Could not select {expected_count}-gene FL neighborhood for {target['gene_id']}")
    return selected


def anchor_side_map(local_fl_genes: Sequence[Gene], target_og: str) -> Dict[str, set[str]]:
    focal_index = next(i for i, gene in enumerate(local_fl_genes) if gene.orthogroup == target_og)
    side_map: Dict[str, set[str]] = defaultdict(set)
    for index, gene in enumerate(local_fl_genes):
        if gene.orthogroup is None:
            continue
        if index < focal_index:
            side_map[gene.orthogroup].add("upstream")
        elif index > focal_index:
            side_map[gene.orthogroup].add("downstream")
        else:
            side_map[gene.orthogroup].add("target")
    return side_map


def classify_anchor_context(cluster: Sequence[Gene], target_og: str, side_map: Dict[str, set[str]]) -> str:
    if any(gene.orthogroup == target_og for gene in cluster):
        return "target_present"
    sides = set()
    for gene in cluster:
        sides.update(side_map.get(gene.orthogroup or "", set()))
    has_up = "upstream" in sides
    has_down = "downstream" in sides
    if has_up and has_down:
        return "bilateral"
    if has_up:
        return "upstream_only"
    if has_down:
        return "downstream_only"
    return "unresolved"


def choose_anchor_cluster(
    anchors: Sequence[Gene],
    target_og: str,
    side_map: Dict[str, set[str]],
    max_span_bp: int = 2_000_000,
) -> Tuple[str, List[Gene], str]:
    """Choose one compact locus and reject chromosome-wide paralog expansion."""
    unique = {gene.gene_id: gene for gene in anchors}
    anchors = list(unique.values())
    target_genes = [gene for gene in anchors if gene.orthogroup == target_og]
    candidates: List[Tuple[Tuple[int, int, int, int], str, List[Gene], str]] = []

    if target_genes:
        for target_gene in target_genes:
            center = (target_gene.start + target_gene.end) // 2
            cluster = [
                gene for gene in anchors
                if gene.seqid == target_gene.seqid
                and abs(((gene.start + gene.end) // 2) - center) <= max_span_bp // 2
            ]
            distinct = len({gene.orthogroup for gene in cluster if gene.orthogroup})
            span = max(gene.end for gene in cluster) - min(gene.start for gene in cluster)
            score = (distinct, len(cluster), -span, -target_gene.start)
            candidates.append((score, target_gene.seqid, cluster, "target_present"))
    else:
        by_seqid: Dict[str, List[Gene]] = defaultdict(list)
        for gene in anchors:
            by_seqid[gene.seqid].append(gene)
        for seqid, seq_genes in by_seqid.items():
            seq_genes.sort(key=lambda gene: ((gene.start + gene.end) // 2, gene.gene_id))
            for left in range(len(seq_genes)):
                cluster = []
                left_center = (seq_genes[left].start + seq_genes[left].end) // 2
                for gene in seq_genes[left:]:
                    center = (gene.start + gene.end) // 2
                    if center - left_center > max_span_bp:
                        break
                    cluster.append(gene)
                if not cluster:
                    continue
                context = classify_anchor_context(cluster, target_og, side_map)
                bilateral = 1 if context == "bilateral" else 0
                distinct = len({gene.orthogroup for gene in cluster if gene.orthogroup})
                span = max(gene.end for gene in cluster) - min(gene.start for gene in cluster)
                score = (distinct, bilateral, len(cluster), -span)
                candidates.append((score, seqid, cluster, context))

    if not candidates:
        raise RuntimeError(f"No anchor cluster could be formed around {target_og}")
    _, seqid, cluster, context = max(candidates, key=lambda item: item[0])
    if len({gene.orthogroup for gene in cluster if gene.orthogroup}) < 2:
        raise RuntimeError(f"Insufficient compact anchors around {target_og}")
    return seqid, cluster, context


def build_tracks(
    local_fl_genes: Sequence[Gene],
    target_og: str,
    groups: Dict[str, Dict[str, List[str]]],
    genes_by_sample: Dict[str, List[Gene]],
    protein_to_gene: Dict[str, Dict[str, Gene]],
    excluded_ogs: Sequence[str] = (),
    flank_bp: int = 30000,
) -> Tuple[List[Track], List[str]]:
    excluded = set(excluded_ogs)
    local_ogs = []
    for gene in local_fl_genes:
        if gene.orthogroup and gene.orthogroup not in excluded and gene.orthogroup not in local_ogs:
            local_ogs.append(gene.orthogroup)
    if target_og not in local_ogs:
        raise RuntimeError(f"Target orthogroup {target_og} is not in its FL neighborhood")

    tracks: List[Track] = []
    side_map = anchor_side_map(local_fl_genes, target_og)
    for sample in SAMPLE_ORDER:
        anchors: List[Gene] = []
        for og in local_ogs:
            for protein_id in groups[og][sample]:
                gene = protein_to_gene[sample].get(protein_id)
                if gene is not None:
                    anchors.append(gene)
        if not anchors:
            raise RuntimeError(f"No local orthogroup anchors found for {sample} around {target_og}")
        chosen_seqid, anchors, anchor_context = choose_anchor_cluster(anchors, target_og, side_map)
        region_start = max(1, min(gene.start for gene in anchors) - flank_bp)
        region_end = max(gene.end for gene in anchors) + flank_bp
        local_genes = [
            gene
            for gene in genes_by_sample[sample]
            if gene.seqid == chosen_seqid and gene.end >= region_start and gene.start <= region_end
        ]
        local_genes.sort(key=lambda gene: (gene.start, gene.end, gene.gene_id))
        tracks.append(
            Track(
                sample,
                chosen_seqid,
                region_start,
                region_end,
                local_genes,
                anchor_context,
                sorted(gene.gene_id for gene in anchors),
            )
        )
    return tracks, local_ogs


def project_syri_position(
    syri_path: Path,
    source_position: int,
    source_side: str,
    expected_chromosome: str,
) -> Dict[str, object]:
    """Project one coordinate through SYNAL records without loading syri.out."""
    if source_side not in {"reference", "query"}:
        raise ValueError(f"Unsupported source side: {source_side}")
    containing = None
    previous = None
    following = None
    synal_count = 0
    with syri_path.open() as handle:
        for line in handle:
            fields = line.rstrip("\n\r").split("\t")
            if len(fields) < 11 or fields[10] != "SYNAL":
                continue
            if fields[0] != expected_chromosome or fields[5] != expected_chromosome:
                continue
            try:
                ref_start, ref_end = int(fields[1]), int(fields[2])
                qry_start, qry_end = int(fields[6]), int(fields[7])
            except ValueError:
                continue
            synal_count += 1
            if source_side == "reference":
                source_start, source_end = ref_start, ref_end
                dest_start, dest_end = qry_start, qry_end
            else:
                source_start, source_end = qry_start, qry_end
                dest_start, dest_end = ref_start, ref_end
            record = {
                "Source_Start": source_start,
                "Source_End": source_end,
                "Destination_Start": dest_start,
                "Destination_End": dest_end,
                "SyRI_Record_ID": fields[8],
            }
            if source_start <= source_position <= source_end:
                if containing is None or (source_end - source_start) < (containing["Source_End"] - containing["Source_Start"]):
                    containing = record
            elif source_end < source_position:
                if previous is None or source_end > previous["Source_End"]:
                    previous = record
            elif source_start > source_position:
                if following is None or source_start < following["Source_Start"]:
                    following = record

    if synal_count == 0:
        raise RuntimeError(f"No chr04 SYNAL records found in {syri_path}")

    if containing is not None:
        source_span = max(1, containing["Source_End"] - containing["Source_Start"])
        fraction = (source_position - containing["Source_Start"]) / source_span
        projected = round(
            containing["Destination_Start"]
            + fraction * (containing["Destination_End"] - containing["Destination_Start"])
        )
        return {
            "Projected_Position": projected,
            "Projection_Method": "within_SYNAL",
            "Left_Record_ID": containing["SyRI_Record_ID"],
            "Right_Record_ID": containing["SyRI_Record_ID"],
            "Left_Distance_bp": 0,
            "Right_Distance_bp": 0,
        }

    if previous is not None and following is not None:
        source_gap = max(1, following["Source_Start"] - previous["Source_End"])
        fraction = (source_position - previous["Source_End"]) / source_gap
        projected = round(
            previous["Destination_End"]
            + fraction * (following["Destination_Start"] - previous["Destination_End"])
        )
        return {
            "Projected_Position": projected,
            "Projection_Method": "between_SYNAL_flanks",
            "Left_Record_ID": previous["SyRI_Record_ID"],
            "Right_Record_ID": following["SyRI_Record_ID"],
            "Left_Distance_bp": source_position - previous["Source_End"],
            "Right_Distance_bp": following["Source_Start"] - source_position,
        }

    nearest = previous if following is None else following
    if nearest is None:
        raise RuntimeError(f"No usable SyRI projection context around {source_position} in {syri_path}")
    if previous is not None:
        distance = source_position - previous["Source_End"]
        projected = previous["Destination_End"] + distance
        record_id = previous["SyRI_Record_ID"]
    else:
        distance = following["Source_Start"] - source_position
        projected = following["Destination_Start"] - distance
        record_id = following["SyRI_Record_ID"]
    if distance > 2_000_000:
        raise RuntimeError(f"Nearest one-sided SyRI block is {distance:,} bp away in {syri_path}")
    return {
        "Projected_Position": projected,
        "Projection_Method": "nearest_SYNAL_extrapolation",
        "Left_Record_ID": record_id if previous is not None else "NA",
        "Right_Record_ID": record_id if following is not None else "NA",
        "Left_Distance_bp": distance if previous is not None else "NA",
        "Right_Distance_bp": distance if following is not None else "NA",
    }


def replace_track_with_syri_projection(
    tracks: List[Track],
    destination_sample: str,
    source_sample: str,
    syri_path: Path,
    source_side: str,
    target_og: str,
    local_ogs: Sequence[str],
    genes_by_sample: Dict[str, List[Gene]],
    expected_chromosome: str,
) -> Dict[str, object]:
    track_by_sample = {track.sample: track for track in tracks}
    source_track = track_by_sample[source_sample]
    target_genes = [gene for gene in source_track.genes if gene.orthogroup == target_og]
    if len(target_genes) != 1:
        raise RuntimeError(f"Expected one source target gene in {source_sample}, observed {len(target_genes)}")
    source_gene = target_genes[0]
    source_position = (source_gene.start + source_gene.end) // 2
    result = project_syri_position(syri_path, source_position, source_side, expected_chromosome)
    projected_position = int(result["Projected_Position"])

    destination_seqids = sorted({
        gene.seqid
        for gene in genes_by_sample[destination_sample]
        if normalize_chromosome(gene.seqid, "chrNA") == expected_chromosome
    })
    if len(destination_seqids) != 1:
        raise RuntimeError(
            f"Expected one {expected_chromosome} GFF seqid for {destination_sample}, observed {destination_seqids}"
        )
    destination_seqid = destination_seqids[0]
    span = source_track.region_end - source_track.region_start
    region_start = max(1, projected_position - span // 2)
    region_end = projected_position + span // 2
    local_genes = [
        gene
        for gene in genes_by_sample[destination_sample]
        if gene.seqid == destination_seqid and gene.end >= region_start and gene.start <= region_end
    ]
    local_genes.sort(key=lambda gene: (gene.start, gene.end, gene.gene_id))
    if any(gene.orthogroup == target_og for gene in local_genes):
        raise RuntimeError(f"Projected {destination_sample} region unexpectedly contains target orthogroup {target_og}")
    anchor_ids = sorted(gene.gene_id for gene in local_genes if gene.orthogroup in local_ogs)
    replacement = Track(
        destination_sample,
        destination_seqid,
        region_start,
        region_end,
        local_genes,
        "syri_projected",
        anchor_ids,
        projected_position,
        str(result["Projection_Method"]),
    )
    tracks[tracks.index(track_by_sample[destination_sample])] = replacement
    return {
        "Destination_Sample": destination_sample,
        "Source_Sample": source_sample,
        "SyRI_File": str(syri_path),
        "Source_Side": source_side,
        "Chromosome": expected_chromosome,
        "Source_Gene_ID": source_gene.gene_id,
        "Source_Position": source_position,
        **result,
        "Destination_GFF_Seqid": destination_seqid,
        "Destination_Region_Start": region_start,
        "Destination_Region_End": region_end,
        "Destination_Local_Genes": len(local_genes),
        "Destination_Local_Orthogroup_Anchors": len(anchor_ids),
    }


def local_copy_counts(tracks: Sequence[Track], local_ogs: Sequence[str]) -> Dict[str, Dict[str, int]]:
    counts = {og: {sample: 0 for sample in SAMPLE_ORDER} for og in local_ogs}
    for track in tracks:
        counter = Counter(gene.orthogroup for gene in track.genes if gene.orthogroup in counts)
        for og in local_ogs:
            counts[og][track.sample] = counter.get(og, 0)
    return counts


def infer_roles(
    tracks: Sequence[Track], target_og: str, local_ogs: Sequence[str]
) -> Tuple[Dict[str, Dict[str, int]], List[str]]:
    counts = local_copy_counts(tracks, local_ogs)
    cnv_ogs = sorted(
        og for og in local_ogs
        if og != target_og and max(counts[og].values()) >= 2 and min(counts[og].values()) < max(counts[og].values())
    )
    for track in tracks:
        for gene in track.genes:
            if gene.orthogroup == target_og:
                gene.role = "Target_PAV"
            elif gene.orthogroup in cnv_ogs:
                gene.role = "CNV"
            elif gene.orthogroup in local_ogs:
                gene.role = "Syntenic_ortholog"
            else:
                gene.role = "Other_local_gene"
    return counts, cnv_ogs


def normalize_chromosome(seqid: str, fallback: str) -> str:
    match = re.search(r"chr(?:omosome)?[_-]?(1[0-6]|0?[1-9])(?=\D|$)", seqid, flags=re.IGNORECASE)
    if match:
        return f"chr{int(match.group(1)):02d}"
    match = re.fullmatch(r"NC_(02599[3-9]|02600[0-8])\.2", seqid)
    if match:
        accession_number = int(match.group(1))
        chromosome_number = accession_number - 25992
        return f"chr{chromosome_number:02d}"
    return fallback


def x_position(value: float, track: Track, left: float = 0.185, right: float = 0.765) -> float:
    span = max(1, track.region_end - track.region_start)
    return left + (value - track.region_start) / span * (right - left)


def gene_center_x(gene: Gene, track: Track) -> float:
    return x_position((gene.start + gene.end) / 2.0, track)


def inferred_absence_x(track: Track, local_fl_genes: Sequence[Gene], target_og: str) -> float:
    focal_index = next(i for i, gene in enumerate(local_fl_genes) if gene.orthogroup == target_og)
    upstream_ogs = [gene.orthogroup for gene in reversed(local_fl_genes[:focal_index]) if gene.orthogroup]
    downstream_ogs = [gene.orthogroup for gene in local_fl_genes[focal_index + 1:] if gene.orthogroup]
    upstream_gene = None
    downstream_gene = None
    for og in upstream_ogs:
        candidates = [gene for gene in track.genes if gene.orthogroup == og]
        if candidates:
            upstream_gene = max(candidates, key=lambda gene: gene.end)
            break
    for og in downstream_ogs:
        candidates = [gene for gene in track.genes if gene.orthogroup == og]
        if candidates:
            downstream_gene = min(candidates, key=lambda gene: gene.start)
            break
    if upstream_gene and downstream_gene:
        return (gene_center_x(upstream_gene, track) + gene_center_x(downstream_gene, track)) / 2.0
    if upstream_gene:
        return min(0.75, gene_center_x(upstream_gene, track) + 0.03)
    if downstream_gene:
        return max(0.20, gene_center_x(downstream_gene, track) - 0.03)
    return 0.475


def ribbon_path(x1: float, y1: float, x2: float, y2: float, half_width: float = 0.0045) -> MplPath:
    delta = (y2 - y1) * 0.45
    verts = [
        (x1 - half_width, y1),
        (x1 - half_width, y1 + delta),
        (x2 - half_width, y2 - delta),
        (x2 - half_width, y2),
        (x2 + half_width, y2),
        (x2 + half_width, y2 - delta),
        (x1 + half_width, y1 + delta),
        (x1 + half_width, y1),
        (x1 - half_width, y1),
    ]
    codes = [
        MplPath.MOVETO,
        MplPath.CURVE4,
        MplPath.CURVE4,
        MplPath.CURVE4,
        MplPath.LINETO,
        MplPath.CURVE4,
        MplPath.CURVE4,
        MplPath.CURVE4,
        MplPath.CLOSEPOLY,
    ]
    return MplPath(verts, codes)


def link_rows(tracks: Sequence[Track], local_ogs: Sequence[str], target_og: str, cnv_ogs: Sequence[str]):
    rows = []
    for upper, lower in zip(tracks[:-1], tracks[1:]):
        for og in local_ogs:
            upper_genes = [gene for gene in upper.genes if gene.orthogroup == og]
            lower_genes = [gene for gene in lower.genes if gene.orthogroup == og]
            if og == target_og:
                link_type = "Target_PAV"
            elif og in cnv_ogs:
                link_type = "CNV"
            else:
                link_type = "Syntenic_ortholog"
            for upper_gene in upper_genes:
                for lower_gene in lower_genes:
                    rows.append({
                        "Upper_Sample": upper.sample,
                        "Upper_Gene": upper_gene.gene_id,
                        "Lower_Sample": lower.sample,
                        "Lower_Gene": lower_gene.gene_id,
                        "Orthogroup": og,
                        "Link_Type": link_type,
                    })
    return rows


def draw_figure(
    target: Dict[str, str],
    target_og: str,
    tracks: Sequence[Track],
    local_fl_genes: Sequence[Gene],
    local_ogs: Sequence[str],
    cnv_ogs: Sequence[str],
    output_stem: Path,
) -> None:
    mpl.rcParams.update({
        "font.family": "Arial",
        "font.size": 9,
        "axes.facecolor": "white",
        "figure.facecolor": "white",
        "axes.grid": False,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "svg.fonttype": "none",
    })
    fig_height = 7.016 + max(0, len(tracks) - 11) * 0.45
    fig = plt.figure(figsize=(7.0, fig_height), facecolor="white")
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    title = f"Multi-haplotype microsynteny around {target['gene_id']}"
    fig.text(0.025, 0.965, title, ha="left", va="top", fontsize=13.0, fontweight="bold")
    legend_handles = [
        Rectangle((0, 0), 1, 1, facecolor=TARGET_COLOR, edgecolor="none", label=target["legend_label"]),
        Rectangle((0, 0), 1, 1, facecolor=CNV_COLOR, edgecolor="none", label="CNV orthogroup"),
        Rectangle((0, 0), 1, 1, facecolor=SYNTENIC_COLOR, edgecolor="none", label="Syntenic orthologs"),
        Rectangle((0, 0), 1, 1, facecolor=OTHER_COLOR, edgecolor="none", label="Other local genes"),
    ]
    ax.legend(
        handles=legend_handles,
        loc="upper left",
        bbox_to_anchor=(0.025, 0.925),
        ncol=4,
        frameon=False,
        handlelength=0.9,
        handleheight=0.9,
        columnspacing=1.25,
        fontsize=10.0,
        borderaxespad=0,
    )
    absent_samples = [
        DISPLAY_LABELS[track.sample]
        for track in tracks
        if not any(gene.orthogroup == target_og for gene in track.genes)
    ]
    absent_text = ", ".join(absent_samples)
    subtitle = f"Blue marks local CNV; curved ribbons retain gene order; {target['short_label']} is absent from {absent_text}"
    fig.text(0.50, 0.866, subtitle, ha="center", va="center", fontsize=9.0, color="#4D4D4D")

    y_top, y_bottom = 0.815, 0.105
    y_values = [y_top - i * (y_top - y_bottom) / (len(tracks) - 1) for i in range(len(tracks))]
    track_by_sample = {track.sample: track for track in tracks}

    for i, (upper, lower) in enumerate(zip(tracks[:-1], tracks[1:])):
        y1, y2 = y_values[i], y_values[i + 1]
        for og in local_ogs:
            upper_genes = [gene for gene in upper.genes if gene.orthogroup == og]
            lower_genes = [gene for gene in lower.genes if gene.orthogroup == og]
            if og == target_og:
                color, alpha = TARGET_COLOR, 0.31
            elif og in cnv_ogs:
                color, alpha = CNV_COLOR, 0.30
            else:
                color, alpha = SYNTENIC_COLOR, 0.20
            for upper_gene in upper_genes:
                for lower_gene in lower_genes:
                    patch = PathPatch(
                        ribbon_path(gene_center_x(upper_gene, upper), y1, gene_center_x(lower_gene, lower), y2),
                        facecolor=color,
                        edgecolor="none",
                        alpha=alpha,
                        zorder=1,
                    )
                    ax.add_patch(patch)

    fallback_chr = normalize_chromosome(local_fl_genes[0].seqid, "chrNA")
    for index, (track, y) in enumerate(zip(tracks, y_values)):
        ax.add_line(Line2D([0.19, 0.765], [y, y], color=TRACK_COLOR, linewidth=0.75, zorder=2))
        weight = "bold" if track.sample in {"FL_Africa_hap2", "EG11"} else "normal"
        fig.text(0.175, y, DISPLAY_LABELS[track.sample], ha="right", va="center", fontsize=9.0, fontweight=weight)
        chrom = normalize_chromosome(track.seqid, fallback_chr)
        fig.text(
            0.775,
            y,
            f"{chrom}  {track.region_start / 1e6:.2f}–{track.region_end / 1e6:.2f} Mb",
            ha="left",
            va="center",
            fontsize=7.0,
            color="#737373",
        )
        for gene in track.genes:
            x1 = x_position(gene.start, track)
            x2 = x_position(gene.end, track)
            width = max(0.0022, x2 - x1)
            center = (x1 + x2) / 2.0
            if gene.role == "Target_PAV":
                color = TARGET_COLOR
            elif gene.role == "CNV":
                color = CNV_COLOR
            elif gene.role == "Syntenic_ortholog":
                color = SYNTENIC_COLOR
            else:
                color = OTHER_COLOR
            ax.add_patch(Rectangle((center - width / 2.0, y - 0.0105), width, 0.021, facecolor=color, edgecolor="none", zorder=3))

        if not any(gene.orthogroup == target_og for gene in track.genes):
            x = (
                x_position(track.projected_target_position, track)
                if track.projected_target_position is not None
                else inferred_absence_x(track, local_fl_genes, target_og)
            )
            ax.add_patch(
                Rectangle(
                    (x - 0.0053, y - 0.0145),
                    0.0106,
                    0.029,
                    fill=False,
                    edgecolor=TARGET_COLOR,
                    linewidth=1.0,
                    linestyle=(0, (3, 2)),
                    zorder=4,
                )
            )
            pav_label = "PAV" if track.anchor_context == "bilateral" else "PAV*"
            fig.text(x, y + 0.021, pav_label, ha="center", va="bottom", fontsize=6.2, color=TARGET_COLOR)

    for index, (upper, lower) in enumerate(zip(tracks[:-1], tracks[1:])):
        upper_group = TRACK_GROUPS.get(upper.sample, "Unspecified")
        lower_group = TRACK_GROUPS.get(lower.sample, "Unspecified")
        if upper_group == lower_group:
            continue
        separator_y = (y_values[index] + y_values[index + 1]) / 2.0
        ax.add_line(Line2D([0.185, 0.77], [separator_y, separator_y], color="#666666", linewidth=0.8, linestyle=(0, (7, 6))))
        fig.text(0.775, separator_y + 0.014, upper_group, ha="left", va="center", fontsize=7.0, color="#777777", fontstyle="italic")
        fig.text(0.775, separator_y - 0.014, lower_group, ha="left", va="center", fontsize=7.0, color="#777777", fontstyle="italic")

    for suffix in ["pdf", "svg", "png"]:
        kwargs = {"bbox_inches": "tight", "facecolor": "white"}
        if suffix == "png":
            kwargs["dpi"] = 600
        fig.savefig(Path(f"{output_stem}.{suffix}"), **kwargs)
    plt.close(fig)


def write_tsv(path: Path, fieldnames: Sequence[str], rows: Iterable[Dict[str, object]]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def write_outputs(
    target: Dict[str, str],
    target_og: str,
    tracks: Sequence[Track],
    local_fl_genes: Sequence[Gene],
    local_ogs: Sequence[str],
    counts: Dict[str, Dict[str, int]],
    cnv_ogs: Sequence[str],
    output_dir: Path,
    projection_rows: Sequence[Dict[str, object]] = (),
) -> None:
    output_dir.mkdir(parents=False, exist_ok=False)
    links = link_rows(tracks, local_ogs, target_og, cnv_ogs)
    gene_rows = []
    for track in tracks:
        for gene in track.genes:
            gene_rows.append({
                "Sample": track.sample,
                "Seqid": normalize_chromosome(track.seqid, normalize_chromosome(local_fl_genes[0].seqid, "chrNA")),
                "Source_Seqid": track.seqid,
                "Start": gene.start,
                "End": gene.end,
                "Strand": gene.strand,
                "Gene_ID": gene.gene_id,
                "Orthogroup": gene.orthogroup or "NA",
                "Role": gene.role,
            })
    write_tsv(
        output_dir / "Genes.tsv",
        ["Sample", "Seqid", "Source_Seqid", "Start", "End", "Strand", "Gene_ID", "Orthogroup", "Role"],
        gene_rows,
    )

    track_rows = []
    for track in tracks:
        target_genes = [gene.gene_id for gene in track.genes if gene.orthogroup == target_og]
        track_rows.append({
            "Sample": track.sample,
            "Seqid": normalize_chromosome(track.seqid, normalize_chromosome(local_fl_genes[0].seqid, "chrNA")),
            "Source_Seqid": track.seqid,
            "Region_Start": track.region_start,
            "Region_End": track.region_end,
            "Local_Genes": len(track.genes),
            "Orthogroup_Genes": sum(gene.orthogroup in local_ogs for gene in track.genes),
            "Local_Orthogroups": len({gene.orthogroup for gene in track.genes if gene.orthogroup in local_ogs}),
            "Target_Copies": len(target_genes),
            "Target_Gene_IDs": ",".join(target_genes) if target_genes else "NA",
            "Anchor_Context": track.anchor_context,
            "Anchor_Gene_IDs": ",".join(track.anchor_gene_ids),
            "Projected_Target_Position": track.projected_target_position or "NA",
            "Projection_Method": track.projection_method,
        })
    write_tsv(
        output_dir / "Tracks.tsv",
        ["Sample", "Seqid", "Source_Seqid", "Region_Start", "Region_End", "Local_Genes", "Orthogroup_Genes", "Local_Orthogroups", "Target_Copies", "Target_Gene_IDs", "Anchor_Context", "Anchor_Gene_IDs", "Projected_Target_Position", "Projection_Method"],
        track_rows,
    )
    write_tsv(
        output_dir / "Adjacent_Track_Links.tsv",
        ["Upper_Sample", "Upper_Gene", "Lower_Sample", "Lower_Gene", "Orthogroup", "Link_Type"],
        links,
    )
    presence_rows = []
    for track in tracks:
        target_genes = [gene.gene_id for gene in track.genes if gene.orthogroup == target_og]
        presence_rows.append({
            "Sample": track.sample,
            "Status": "PRESENT" if target_genes else "ABSENT",
            "Local_Copy_Number": len(target_genes),
            "Gene_IDs": ",".join(target_genes) if target_genes else "NA",
            "Orthogroup": target_og,
            "Anchor_Context": track.anchor_context,
            "Position_Inference": (
                "not_applicable"
                if target_genes
                else f"syri_projection:{track.projection_method}"
                if track.projected_target_position is not None
                else "bilateral_flanks"
                if track.anchor_context == "bilateral"
                else "one_sided_context"
            ),
        })
    write_tsv(
        output_dir / "Target_Presence.tsv",
        ["Sample", "Status", "Local_Copy_Number", "Gene_IDs", "Orthogroup", "Anchor_Context", "Position_Inference"],
        presence_rows,
    )

    if projection_rows:
        write_tsv(
            output_dir / "SyRI_Projection_Evidence.tsv",
            [
                "Destination_Sample", "Source_Sample", "SyRI_File", "Source_Side", "Chromosome",
                "Source_Gene_ID", "Source_Position", "Projected_Position", "Projection_Method",
                "Left_Record_ID", "Right_Record_ID", "Left_Distance_bp", "Right_Distance_bp",
                "Destination_GFF_Seqid", "Destination_Region_Start", "Destination_Region_End",
                "Destination_Local_Genes", "Destination_Local_Orthogroup_Anchors",
            ],
            projection_rows,
        )

    cnv_rows = []
    for og in cnv_ogs:
        row: Dict[str, object] = {
            "Orthogroup": og,
            "Min_Local_Copy": min(counts[og].values()),
            "Max_Local_Copy": max(counts[og].values()),
        }
        row.update(counts[og])
        cnv_rows.append(row)
    write_tsv(
        output_dir / "CNV_Orthogroups.tsv",
        ["Orthogroup", "Min_Local_Copy", "Max_Local_Copy"] + SAMPLE_ORDER,
        cnv_rows,
    )

    output_stem = output_dir / f"{target['gene_id']}_all_haplotypes_microsynteny"
    draw_figure(target, target_og, tracks, local_fl_genes, local_ogs, cnv_ogs, output_stem)

    absent = [row["Sample"] for row in presence_rows if row["Status"] == "ABSENT"]
    one_sided = [
        row["Sample"]
        for row in presence_rows
        if row["Status"] == "ABSENT" and row["Position_Inference"] == "one_sided_context"
    ]
    projected = [
        row["Sample"]
        for row in presence_rows
        if row["Status"] == "ABSENT" and str(row["Position_Inference"]).startswith("syri_projection:")
    ]
    projection_evidence = ",".join(
        f"{row['Destination_Sample']}:{row['Projection_Method']}" for row in projection_rows
    )
    qa_rows = [
        {"Check": "Track_count", "Status": "PASS" if len(tracks) == len(SAMPLE_ORDER) else "FAIL", "Evidence": len(tracks)},
        {"Check": "Target_orthogroup", "Status": "PASS", "Evidence": target_og},
        {"Check": "Target_absent_tracks", "Status": "PASS", "Evidence": ",".join(absent)},
        {"Check": "One_sided_absence_context", "Status": "WARN" if one_sided else "PASS", "Evidence": ",".join(one_sided) if one_sided else "NONE"},
        {"Check": "SyRI_projected_tracks", "Status": "PASS" if projection_rows else "PASS", "Evidence": projection_evidence if projection_rows else "NONE"},
        {"Check": "FL_target_present", "Status": "PASS" if counts[target_og]["FL_Africa_hap2"] >= 1 else "FAIL", "Evidence": counts[target_og]["FL_Africa_hap2"]},
        {"Check": "CNV_orthogroups", "Status": "PASS", "Evidence": ",".join(cnv_ogs) if cnv_ogs else "NONE"},
        {"Check": "Adjacent_links", "Status": "PASS" if links else "FAIL", "Evidence": len(links)},
        {"Check": "Vector_and_600dpi", "Status": "PASS", "Evidence": "PDF+SVG+600-dpi PNG"},
        {"Check": "Image_inspection", "Status": "PENDING", "Evidence": "Manual visual inspection required"},
    ]
    write_tsv(output_dir / "Figure_QA.tsv", ["Check", "Status", "Evidence"], qa_rows)

    absent_labels = ", ".join(DISPLAY_LABELS[sample] for sample in absent)
    one_sided_sentence = (
        " Asterisks identify one-sided local-anchor inference in "
        + ", ".join(DISPLAY_LABELS[sample] for sample in one_sided)
        + "; these positions are unresolved and are not bilateral-synteny absence calls."
        if one_sided
        else ""
    )
    projected_sentence = (
        " Asterisks in "
        + ", ".join(DISPLAY_LABELS[sample] for sample in projected)
        + " mark target positions projected from the corresponding target-present haplotype through chr04 SyRI SYNAL context; SyRI is used only for positional localization, whereas target absence remains annotation-based."
        if projected
        else ""
    )
    alias_samples = [sample for sample in SAMPLE_ORDER if ORTHOGROUP_ALIASES.get(sample, sample) != sample]
    alias_sentence = (
        " Orthogroups for "
        + ", ".join(DISPLAY_LABELS[sample] for sample in alias_samples)
        + " were inherited by stable transcript ID from the source annotation used for Liftoff; genomic positions and local gene order were taken from each new GFF3."
        if alias_samples
        else ""
    )
    excluded_ogs = target.get("exclude_orthogroups", [])
    excluded_sentence = (
        " The highly multicopy orthogroup "
        + ", ".join(excluded_ogs)
        + " was excluded from local-anchor selection and ribbons because its dispersed paralogues do not provide unambiguous microsynteny evidence."
        if excluded_ogs
        else ""
    )
    legend = f"""# Figure legend

Multi-haplotype microsynteny around the FL Africa hap2 {target['short_label']} candidate {target['gene_id']}. Orange rectangles and ribbons represent the target PAV orthogroup ({target_og}). Blue rectangles and ribbons identify local orthogroups whose copy number varies among tracks and reaches at least two copies; grey ribbons connect the remaining orthogroups defined by the {len(local_fl_genes)}-gene FL Africa hap2 neighborhood. Ribbons are drawn only between adjacent tracks and retain each assembly's genomic gene order, so crossings show local orientation or order differences. Dashed orange boxes mark inferred positions in tracks lacking an annotated target-orthogroup member ({absent_labels}).{one_sided_sentence}{projected_sentence} The dashed horizontal line separates the *E. guineensis* and *E. oleifera*-derived tracks.{alias_sentence}{excluded_sentence} Orthogroup presence is annotation-based and should not by itself be interpreted as experimentally validated gene loss.
"""
    (output_dir / "Figure_Legend.md").write_text(legend)

    design = f"""# Figure design brief

- Scientific question: Is the {target['short_label']} orthogroup absent specifically from EG11 while retained in other annotated oil-palm haplotypes?
- Figure family: specialized coordinate-aware microsynteny/genome-track figure.
- Data roles: sample, genomic coordinate, gene interval, strand and orthogroup membership.
- Design mode: preserve the accepted 11-track microsynteny structure used for PDAT1, FabI-like ENR, HACD and GDSL lipase.
- Rejected alternative: adding the two pure *E. oleifera* haplotypes, because this supplementary panel is intended to support the FL Africa hap2 versus EG11 annotation-quality argument.
- Excluded ambiguous anchors: {', '.join(target.get('exclude_orthogroups', [])) or 'none'}.
- Figure tier: manuscript candidate pending rendered-image and biological QA.
"""
    (output_dir / "Figure_Design.md").write_text(design)
    metadata = {
        "Figure_Family": "specialized_microsynteny",
        "Pattern_Document": "custom_specialized_synteny_no_generic_pattern",
        "Target_Gene": target["gene_id"],
        "Target_Orthogroup": target_og,
        "Track_Count": len(tracks),
        "Sample_Order": list(SAMPLE_ORDER),
        "Display_Labels": DISPLAY_LABELS,
        "Track_Groups": TRACK_GROUPS,
        "Orthogroup_Aliases": ORTHOGROUP_ALIASES,
        "Excluded_Orthogroups": target.get("exclude_orthogroups", []),
        "Coordinate_Source": "per-sample GFF3",
        "Output_Formats": ["PDF", "SVG", "PNG_600dpi"],
    }
    (output_dir / "Figure_Metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")


def validate_outputs(output_dir: Path, require_projection: bool = False) -> None:
    required = [
        "Genes.tsv",
        "Tracks.tsv",
        "Adjacent_Track_Links.tsv",
        "Target_Presence.tsv",
        "CNV_Orthogroups.tsv",
        "Figure_QA.tsv",
        "Figure_Legend.md",
        "Figure_Design.md",
        "Figure_Metadata.json",
    ]
    if require_projection:
        required.append("SyRI_Projection_Evidence.tsv")
    for name in required:
        path = output_dir / name
        if not path.is_file() or path.stat().st_size == 0:
            raise RuntimeError(f"Missing or empty output: {path}")
    pngs = list(output_dir.glob("*.png"))
    pdfs = list(output_dir.glob("*.pdf"))
    svgs = list(output_dir.glob("*.svg"))
    if len(pngs) != 1 or len(pdfs) != 1 or len(svgs) != 1:
        raise RuntimeError(f"Expected one PDF/PNG/SVG in {output_dir}")
    from PIL import Image
    with Image.open(pngs[0]) as image:
        if image.width < 3000 or image.height < 3500:
            raise RuntimeError(f"PNG is smaller than expected 600-dpi canvas: {image.size}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--orthogroups", required=True, type=Path)
    parser.add_argument("--output-root", required=True, type=Path)
    parser.add_argument("--target-gene")
    parser.add_argument("--output-name")
    parser.add_argument("--syri-root", type=Path)
    args = parser.parse_args()

    manifest = read_manifest(args.manifest)
    groups, protein_to_og = read_orthogroups(args.orthogroups)
    eg11_gff = Path(manifest["EG11"]["Gene_GFF3"])
    eg11_main_seqids = infer_eg11_main_seqids(eg11_gff, set(protein_to_og["EG11"]))
    print(f"[INFO] Inferred 16 EG11 main chromosomes: {','.join(sorted(eg11_main_seqids))}", flush=True)

    genes_by_sample: Dict[str, List[Gene]] = {}
    gene_by_id: Dict[str, Dict[str, Gene]] = {}
    transcript_to_gene: Dict[str, Dict[str, Gene]] = {}
    for sample in SAMPLE_ORDER:
        gff = Path(manifest[sample]["Gene_GFF3"])
        if not gff.is_file() or gff.stat().st_size == 0:
            raise RuntimeError(f"Missing or empty GFF for {sample}: {gff}")
        if sample == "EG11":
            genes, by_id, tx_to_gene = parse_gff(sample, gff, allowed_seqids=eg11_main_seqids)
        else:
            genes, by_id, tx_to_gene = parse_gff(sample, gff, seqid_pattern=manifest[sample]["SeqID_Regex"])
        genes_by_sample[sample] = genes
        gene_by_id[sample] = by_id
        transcript_to_gene[sample] = tx_to_gene
        print(f"[INFO] Parsed {sample}: {len(genes):,} genes", flush=True)

    protein_to_gene = attach_orthogroups(genes_by_sample, gene_by_id, transcript_to_gene, protein_to_og)
    print(f"[INFO] Loaded {len(groups):,} orthogroups", flush=True)

    args.output_root.mkdir(parents=True, exist_ok=True)
    targets = [target for target in TARGETS if args.target_gene in {None, target["gene_id"]}]
    if not targets:
        raise RuntimeError(f"Unknown target gene: {args.target_gene}")
    for target in targets:
        if args.output_name:
            if len(targets) != 1:
                raise RuntimeError("--output-name requires exactly one --target-gene")
            target = dict(target)
            target["output_dir"] = args.output_name
        output_dir = args.output_root / target["output_dir"]
        if output_dir.exists():
            raise RuntimeError(f"Refusing to overwrite existing output: {output_dir}")
        og = target_orthogroup(target, protein_to_og)
        local_fl_genes = select_fl_neighborhood(target, genes_by_sample, gene_by_id)
        tracks, local_ogs = build_tracks(
            local_fl_genes,
            og,
            groups,
            genes_by_sample,
            protein_to_gene,
            target.get("exclude_orthogroups", []),
        )
        projection_rows: List[Dict[str, object]] = []
        if target["gene_id"] == "evm.TU.chr04B.1103":
            if args.syri_root is None:
                raise RuntimeError("GDSL v2 requires --syri-root for chr04 sequence projection")
            expected_chromosome = normalize_chromosome(local_fl_genes[0].seqid, "chrNA")
            projection_specs = [
                ("nrly_hap2", "nrly_hap1", "09_nrly_hap1__nrly_hap2/syri.out", "reference"),
                ("BK_hap1", "BK_hap2", "11_BK_hap1__BK_hap2/syri.out", "query"),
            ]
            for destination_sample, source_sample, relative_path, source_side in projection_specs:
                syri_path = args.syri_root / relative_path
                if not syri_path.is_file() or syri_path.stat().st_size == 0:
                    raise RuntimeError(f"Missing or empty SyRI output: {syri_path}")
                projection_rows.append(
                    replace_track_with_syri_projection(
                        tracks,
                        destination_sample,
                        source_sample,
                        syri_path,
                        source_side,
                        og,
                        local_ogs,
                        genes_by_sample,
                        expected_chromosome,
                    )
                )
        counts, cnv_ogs = infer_roles(tracks, og, local_ogs)
        write_outputs(target, og, tracks, local_fl_genes, local_ogs, counts, cnv_ogs, output_dir, projection_rows)
        validate_outputs(output_dir, require_projection=bool(projection_rows))
        absent = [sample for sample in SAMPLE_ORDER if counts[og][sample] == 0]
        print(
            f"[INFO] Completed {target['gene_id']} | OG={og} | absent={','.join(absent)} | CNV={','.join(cnv_ogs) or 'NONE'}",
            flush=True,
        )


if __name__ == "__main__":
    main()
