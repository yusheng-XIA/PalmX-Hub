#!/usr/bin/env python3
"""Build an exploratory 32-chromosome diploid ancestry/breeding atlas."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import platform
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Patch
import numpy as np
import pandas as pd


PAIRING = {
    "TN": ("TN_h1", "TN_h2"),
    "FL": ("FL_HapA", "FL_HapB"),
}
CALL_MAP = {
    "Dura_like": "Dura",
    "Pisifera_like": "Pisifera",
    "Meizhou4_like": "Meizhou4",
    "TN_derived_unresolved": "TN-derived unresolved",
    "Mixed": "Mixed",
    "Unknown_low_support": "Unknown",
}
PRIMITIVE = {"Dura", "Pisifera", "Meizhou4"}
COLORS = {
    "Dura": "#0072B2",
    "Pisifera": "#E69F00",
    "Meizhou4": "#009E73",
    "TN-derived unresolved": "#8E6C8A",
    "Mixed": "#CC79A7",
    "Unknown": "#D9D9D9",
}
STATE_COLORS = {
    "Ancestry-homozygous": "#4C956C",
    "Ancestry-heterozygous": "#6C5CE7",
    "Unresolved": "#BDBDBD",
}
TRAIT_COLORS = {
    "Oil amount/composition": "#D55E00",
    "Fruit/mesocarp": "#56B4E9",
    "Rancidity risk": "#B2182B",
}


def args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--ancestry-windows", required=True, type=Path)
    p.add_argument("--ase-overlap", required=True, type=Path)
    p.add_argument("--lipid-candidates", required=True, type=Path)
    p.add_argument("--mesocarp-candidates", required=True, type=Path)
    p.add_argument("--rancidity-candidates", required=True, type=Path)
    p.add_argument("--heterosis-gene-stage", required=True, type=Path)
    p.add_argument("--trait-gene-catalog", required=True, type=Path)
    p.add_argument("--run-root", required=True, type=Path)
    p.add_argument("--bins", type=int, default=100)
    p.add_argument("--long-block-percent", type=float, default=5.0)
    return p.parse_args()


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def norm_gene_id(value: object) -> str:
    if pd.isna(value):
        return ""
    return str(value).replace("evm.TU.", "evm.model.")


def ancestry_state(a: str, b: str) -> tuple[str, str]:
    if a in PRIMITIVE and b in PRIMITIVE:
        if a == b:
            return "Ancestry-homozygous", f"{a}/{b}"
        return "Ancestry-heterozygous", "/".join(sorted((a, b)))
    return "Unresolved", f"{a}/{b}"


def weighted_call(sub: pd.DataFrame, lo: float, hi: float, chr_len: int) -> str:
    starts = sub["Start0"].to_numpy(float) / chr_len
    ends = sub["End0"].to_numpy(float) / chr_len
    overlap = np.maximum(0.0, np.minimum(ends, hi) - np.maximum(starts, lo))
    if overlap.sum() == 0:
        center = (lo + hi) / 2
        mids = (starts + ends) / 2
        return sub.iloc[int(np.argmin(np.abs(mids - center)))]["Ancestry"]
    scores: dict[str, float] = {}
    for call, weight in zip(sub["Ancestry"], overlap):
        scores[call] = scores.get(call, 0.0) + float(weight)
    return max(scores, key=scores.get)


def build_bins(windows: pd.DataFrame, n_bins: int) -> tuple[pd.DataFrame, pd.DataFrame]:
    windows = windows.copy()
    windows["Ancestry"] = windows["Final_ancestry_call"].map(CALL_MAP).fillna("Unknown")
    lengths = (windows.groupby(["Target_ID", "Chromosome"])["End0"].max()
               .rename("Chromosome_length_bp").reset_index())
    lookup = {}
    for (target, chrom), sub in windows.groupby(["Target_ID", "Chromosome"], sort=False):
        lookup[(target, chrom)] = sub.sort_values("Start0")

    rows = []
    for individual, (hap1, hap2) in PAIRING.items():
        chr1 = set(windows.loc[windows.Target_ID == hap1, "Chromosome"])
        chr2 = set(windows.loc[windows.Target_ID == hap2, "Chromosome"])
        chroms = sorted(chr1 & chr2, key=lambda x: int(str(x).replace("chr", "")))
        if len(chroms) != 16:
            raise ValueError(f"{individual}: expected 16 paired chromosomes, found {len(chroms)}")
        for chrom in chroms:
            s1, s2 = lookup[(hap1, chrom)], lookup[(hap2, chrom)]
            l1, l2 = int(s1.End0.max()), int(s2.End0.max())
            for i in range(n_bins):
                lo, hi = i / n_bins, (i + 1) / n_bins
                a1 = weighted_call(s1, lo, hi, l1)
                a2 = weighted_call(s2, lo, hi, l2)
                state, genotype = ancestry_state(a1, a2)
                rows.append({
                    "Individual": individual, "Chromosome": chrom,
                    "Bin_index": i, "Fraction_start": lo, "Fraction_end": hi,
                    "Hap1_ID": hap1, "Hap2_ID": hap2,
                    "Hap1_start0": round(lo*l1), "Hap1_end0": round(hi*l1),
                    "Hap2_start0": round(lo*l2), "Hap2_end0": round(hi*l2),
                    "Hap1_ancestry": a1, "Hap2_ancestry": a2,
                    "Diploid_state": state, "Ancestry_genotype": genotype,
                    "Coordinate_method": "chromosome_fraction_1pct_exploratory",
                })
    out = pd.DataFrame(rows)
    return out, lengths


def merge_blocks(bins: pd.DataFrame) -> pd.DataFrame:
    rows = []
    keys = ["Individual", "Chromosome"]
    for (individual, chrom), sub in bins.groupby(keys, sort=False):
        sub = sub.sort_values("Bin_index").reset_index(drop=True)
        group = (sub[["Hap1_ancestry", "Hap2_ancestry", "Diploid_state"]]
                 .ne(sub[["Hap1_ancestry", "Hap2_ancestry", "Diploid_state"]].shift())
                 .any(axis=1).cumsum())
        for _, b in sub.groupby(group, sort=False):
            first, last = b.iloc[0], b.iloc[-1]
            rows.append({
                "Individual": individual, "Chromosome": chrom,
                "Fraction_start": first.Fraction_start,
                "Fraction_end": last.Fraction_end,
                "Length_percent_chromosome": 100*(last.Fraction_end-first.Fraction_start),
                "Hap1_ID": first.Hap1_ID, "Hap2_ID": first.Hap2_ID,
                "Hap1_start0": int(first.Hap1_start0), "Hap1_end0": int(last.Hap1_end0),
                "Hap2_start0": int(first.Hap2_start0), "Hap2_end0": int(last.Hap2_end0),
                "Hap1_ancestry": first.Hap1_ancestry,
                "Hap2_ancestry": first.Hap2_ancestry,
                "Diploid_state": first.Diploid_state,
                "Ancestry_genotype": first.Ancestry_genotype,
                "Coordinate_method": first.Coordinate_method,
            })
    return pd.DataFrame(rows)


def candidate_queries(lipid: pd.DataFrame, meso: pd.DataFrame,
                      rancid: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for _, r in lipid.iterrows():
        for individual, col in (("FL", "FL_gene_id"), ("TN", "TN_gene_id")):
            gene = norm_gene_id(r.get(col, ""))
            if gene:
                rows.append({
                    "Candidate_source": "Frozen_Fig2", "Individual": individual,
                    "Query_gene_ID": gene, "Preferred_name": "",
                    "Trait_module": "Oil amount/composition",
                    "Product": r.get("product", ""),
                    "Existing_priority": r.get("figure2_evidence_count", np.nan),
                    "Design_test": "Compare source/source, source/alternative, and alternative/alternative diplotypes",
                })
    top = meso.sort_values("priority_score", ascending=False).head(15)
    for _, r in top.iterrows():
        for individual, col in (("FL", "gene_africa"), ("TN", "gene_dura")):
            gene = norm_gene_id(r.get(col, ""))
            if gene:
                rows.append({
                    "Candidate_source": "Top15_mesocarp", "Individual": individual,
                    "Query_gene_ID": gene,
                    "Preferred_name": r.get("preferred_name", ""),
                    "Trait_module": "Fruit/mesocarp", "Product": r.get("description", ""),
                    "Existing_priority": r.get("priority_score", np.nan),
                    "Design_test": "Test homozygous and heterozygous diplotypes for fruit/yield phenotype",
                })
    for _, r in rancid.iterrows():
        rows.append({
            "Candidate_source": "Direct_rancidity", "Individual": "TN",
            "Query_gene_ID": norm_gene_id(r.get("reference_gene_alias", "")),
            "Preferred_name": r.get("trait_preferred_name", ""),
            "Trait_module": "Rancidity risk", "Product": r.get("product", ""),
            "Existing_priority": r.get("robust_B_stage_n_TN", np.nan),
            "Design_test": "Avoid/suppress risk allele only after replicated FFA, acid-value, peroxide and enzyme validation",
        })
    return pd.DataFrame(rows).drop_duplicates(["Individual", "Query_gene_ID", "Trait_module"])


def map_candidates(queries: pd.DataFrame, ase: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    ase = ase.copy()
    ase["ASE_gene_norm"] = ase["ASE_gene_ID"].map(norm_gene_id)
    ase["Genome_gene_norm"] = ase["Genome_mRNA_ID"].map(norm_gene_id)
    mapped, missing = [], []
    for _, q in queries.iterrows():
        sub = ase[(ase.Individual == q.Individual) &
                  ((ase.ASE_gene_norm == q.Query_gene_ID) |
                   (ase.Genome_gene_norm == q.Query_gene_ID))]
        if sub.empty:
            missing.append(q.to_dict())
            continue
        pair_id = sub.iloc[0].Allele_pair
        pair = ase[(ase.Individual == q.Individual) & (ase.Allele_pair == pair_id)].copy()
        # Prefer one representative record per phased target.
        pair = pair.sort_values(["Target_ID", "ASE_sample_count"], ascending=[True, False])
        pair = pair.drop_duplicates("Target_ID")
        rec = q.to_dict()
        rec["Allele_pair"] = pair_id
        rec["Chromosome"] = str(pair.iloc[0].Chromosome)
        rec["Mapped_homolog_count"] = len(pair)
        for idx, (_, a) in enumerate(pair.iterrows(), 1):
            rec[f"Homolog{idx}_target"] = a.Target_ID
            rec[f"Homolog{idx}_gene"] = a.ASE_gene_ID
            rec[f"Homolog{idx}_start0"] = int(a.Start0)
            rec[f"Homolog{idx}_end0"] = int(a.End0)
            rec[f"Homolog{idx}_ancestry"] = CALL_MAP.get(a.Ancestry_call, "Unknown")
            rec[f"Homolog{idx}_ASE"] = a.ASE_in_at_least_one_sample
            rec[f"Homolog{idx}_ASE_sample_count"] = int(a.ASE_sample_count)
            rec[f"Homolog{idx}_ASE_bias_count"] = int(a.ASE_bias_toward_this_haplotype_count)
        if len(pair) == 2:
            state, genotype = ancestry_state(rec["Homolog1_ancestry"], rec["Homolog2_ancestry"])
        else:
            state, genotype = "Unresolved", "Unpaired_gene_model"
        rec["Candidate_diploid_state"] = state
        rec["Candidate_ancestry_genotype"] = genotype
        rec["Evidence_boundary"] = "ASE+local ancestry; no direct proof of heterosis or optimum diplotype"
        mapped.append(rec)
    return pd.DataFrame(mapped), pd.DataFrame(missing)


def heterosis_ancestry_context(heterosis: pd.DataFrame, trait_catalog: pd.DataFrame,
                               ase: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    gene_col = "gene_id_africa_hap2"
    if gene_col not in heterosis or "effect_class" not in heterosis:
        raise ValueError("heterosis gene-stage input lacks required columns")
    h = heterosis.copy()
    h["above_better"] = h.effect_class.eq("above_better_parent")
    summary = (h.groupby(gene_col).agg(
        Stages_tested=("stage", "nunique"),
        Above_better_parent_stages=("above_better", "sum"),
        Median_log2_MPV_effect=("log2_mid_parent_effect", "median"),
        Max_log2_HPV_effect=("log2_better_parent_effect", "max"),
    ).reset_index().rename(columns={gene_col: "gene_africa"}))
    positive = h[h.above_better].groupby(gene_col)["log2_better_parent_effect"].median()
    summary["Median_log2_HPV_effect_among_ABPH_stages"] = summary.gene_africa.map(positive)
    summary["Recurrent_ABPH_ge3_stages"] = summary.Above_better_parent_stages.ge(3)

    tc = trait_catalog.drop_duplicates("gene_africa")
    keep = [c for c in ["gene_africa", "family", "preferred_name", "description", "trait_module"] if c in tc]
    summary = summary.merge(tc[keep], on="gene_africa", how="left")
    summary["gene_norm"] = summary.gene_africa.map(norm_gene_id)

    tn = ase[ase.Individual.eq("TN")].copy()
    tn["ASE_gene_norm"] = tn.ASE_gene_ID.map(norm_gene_id)
    gene_to_pair = (tn.sort_values("ASE_sample_count", ascending=False)
                    .drop_duplicates("ASE_gene_norm").set_index("ASE_gene_norm")["Allele_pair"].to_dict())
    pair_records = {}
    for pair_id, pair in tn.groupby("Allele_pair"):
        pair = pair.sort_values("Target_ID").drop_duplicates("Target_ID")
        if len(pair) == 2:
            calls = [CALL_MAP.get(x, "Unknown") for x in pair.Ancestry_call]
            state, genotype = ancestry_state(calls[0], calls[1])
            pair_records[pair_id] = {
                "Chromosome": pair.iloc[0].Chromosome,
                "TN_homolog1_gene": pair.iloc[0].Genome_mRNA_ID,
                "TN_homolog2_gene": pair.iloc[1].Genome_mRNA_ID,
                "TN_homolog1_start0": int(pair.iloc[0].Start0),
                "TN_homolog2_start0": int(pair.iloc[1].Start0),
                "TN_homolog1_ancestry": calls[0], "TN_homolog2_ancestry": calls[1],
                "Diploid_state": state, "Ancestry_genotype": genotype,
            }
    summary["Allele_pair"] = summary.gene_norm.map(gene_to_pair)
    context = pd.DataFrame.from_dict(pair_records, orient="index")
    context.index.name = "Allele_pair"; context = context.reset_index()
    summary = summary.merge(context, on="Allele_pair", how="left")
    summary["Diploid_state"] = summary.Diploid_state.fillna("Unresolved")
    summary["Inference_level"] = "exploratory_parent_n1_not_formal_heterosis"
    recurrent = summary[summary.Recurrent_ABPH_ge3_stages].copy()
    recurrent["Trait_relevant"] = recurrent.trait_module.notna()
    state_summary = (summary.groupby(["Diploid_state", "Recurrent_ABPH_ge3_stages"])
                     .size().rename("Gene_count").reset_index())
    totals = summary.groupby("Diploid_state").size().rename("State_total_genes").reset_index()
    state_summary = state_summary.merge(totals, on="Diploid_state")
    state_summary["Percent_within_state"] = 100*state_summary.Gene_count/state_summary.State_total_genes
    state_summary["Inference_level"] = "descriptive_exploratory_parent_n1"
    return summary, recurrent, state_summary


def composition(bins: pd.DataFrame) -> pd.DataFrame:
    c = (bins.groupby(["Individual", "Diploid_state"]).size().rename("Bin_count").reset_index())
    c["Percent"] = c.groupby("Individual")["Bin_count"].transform(lambda x: 100*x/x.sum())
    return c


def rounded_track(ax, y: float, height: float = 0.26) -> FancyBboxPatch:
    patch = FancyBboxPatch((0, y-height/2), 1, height,
                           boxstyle=f"round,pad=0,rounding_size={height/2}",
                           facecolor="#F3F3F3", edgecolor="#7A7A7A", linewidth=0.45)
    ax.add_patch(patch)
    return patch


def draw_homologs(ax, bins: pd.DataFrame, individual: str, title: str) -> None:
    sub = bins[bins.Individual == individual]
    chroms = sorted(sub.Chromosome.unique(), key=lambda x: int(str(x).replace("chr", "")))
    for row, chrom in enumerate(chroms):
        y0 = len(chroms)-1-row
        c = sub[sub.Chromosome == chrom].sort_values("Bin_index")
        for offset, col in ((0.18, "Hap1_ancestry"), (-0.18, "Hap2_ancestry")):
            capsule = rounded_track(ax, y0+offset)
            for _, b in c.iterrows():
                rect = plt.Rectangle((b.Fraction_start, y0+offset-0.12),
                                     b.Fraction_end-b.Fraction_start, 0.24,
                                     facecolor=COLORS[b[col]], edgecolor="none")
                rect.set_clip_path(capsule)
                ax.add_patch(rect)
        ax.text(-0.035, y0, str(chrom).replace("chr", "Chr"), ha="right", va="center", fontsize=6.5)
    ax.set_xlim(-0.08, 1.01); ax.set_ylim(-0.7, len(chroms)-0.3)
    ax.set_xticks([0, .25, .5, .75, 1]); ax.set_xticklabels(["0", "25", "50", "75", "100"])
    ax.set_yticks([]); ax.set_xlabel("Relative chromosome position (%)", fontsize=7)
    ax.set_title(title, loc="left", fontsize=10, fontweight="bold")
    for s in ax.spines.values(): s.set_visible(False)
    ax.tick_params(axis="x", labelsize=6, length=2)


def draw_target_framework(ax, candidates: pd.DataFrame) -> None:
    chroms = [f"chr{i:02d}" for i in range(1, 17)]
    for row, chrom in enumerate(chroms):
        y0 = len(chroms)-1-row
        rounded_track(ax, y0+0.18); rounded_track(ax, y0-0.18)
        ax.text(-0.035, y0, f"Chr{row+1:02d}", ha="right", va="center", fontsize=6.5)
    if not candidates.empty:
        # Candidate positions use their own chromosome assembly length proxy.
        for _, r in candidates.iterrows():
            chrom = str(r.Chromosome).replace("A", "").replace("B", "")
            if chrom not in chroms: continue
            y0 = len(chroms)-1-chroms.index(chrom)
            starts = [r.get("Homolog1_start0", np.nan), r.get("Homolog2_start0", np.nan)]
            vals = [float(x) for x in starts if pd.notna(x)]
            if not vals: continue
            # Use a conservative display denominator; exact homolog coordinates remain in TSV.
            pos = min(0.99, max(0.01, np.mean(vals)/180_000_000))
            marker = "x" if r.Trait_module == "Rancidity risk" else "D"
            ax.scatter([pos], [y0], s=14, marker=marker, color=TRAIT_COLORS[r.Trait_module],
                       linewidths=0.7, zorder=5)
    ax.set_xlim(-0.08, 1.01); ax.set_ylim(-0.7, len(chroms)-0.3)
    ax.set_yticks([]); ax.set_xticks([0, .25, .5, .75, 1]); ax.set_xticklabels(["0", "25", "50", "75", "100"])
    ax.set_xlabel("Schematic relative position (%)", fontsize=7)
    ax.set_title("C  Target 32-chromosome framework\nancestry state must be chosen by diplotype phenotype", loc="left", fontsize=10, fontweight="bold")
    for s in ax.spines.values(): s.set_visible(False)
    ax.tick_params(axis="x", labelsize=6, length=2)


def render_figures(bins: pd.DataFrame, comp: pd.DataFrame, candidates: pd.DataFrame,
                   figure_dir: Path) -> None:
    plt.rcParams.update({"font.family": "DejaVu Sans", "pdf.fonttype": 42, "svg.fonttype": "none"})
    fig = plt.figure(figsize=(13.5, 10.5), constrained_layout=True)
    gs = fig.add_gridspec(2, 3, height_ratios=[4.5, 1.3])
    ax1 = fig.add_subplot(gs[0, 0]); draw_homologs(ax1, bins, "TN", "A  TN phased diploid ancestry (32 chromosomes)")
    ax2 = fig.add_subplot(gs[0, 1]); draw_homologs(ax2, bins, "FL", "B  FL phased diploid ancestry (32 chromosomes)")
    ax3 = fig.add_subplot(gs[0, 2]); draw_target_framework(ax3, candidates)
    ax4 = fig.add_subplot(gs[1, :2])
    pivot = comp.pivot(index="Individual", columns="Diploid_state", values="Percent").fillna(0)
    order = ["Ancestry-homozygous", "Ancestry-heterozygous", "Unresolved"]
    left = np.zeros(len(pivot))
    for state in order:
        vals = pivot.get(state, pd.Series(0, index=pivot.index)).to_numpy()
        ax4.barh(pivot.index, vals, left=left, color=STATE_COLORS[state], height=.55, label=state)
        for y, (l, v) in enumerate(zip(left, vals)):
            if v >= 5: ax4.text(l+v/2, y, f"{v:.1f}%", ha="center", va="center", fontsize=7, color="white")
        left += vals
    ax4.set_xlim(0,100); ax4.set_xlabel("Paired chromosome fraction (%)", fontsize=8)
    ax4.set_title("D  Ancestry state of homologous chromosome fractions", loc="left", fontsize=10, fontweight="bold")
    ax4.legend(frameon=False, ncol=3, fontsize=7, loc="lower center", bbox_to_anchor=(.5,-.48))
    for s in ("top","right","left"): ax4.spines[s].set_visible(False)
    ax4.tick_params(labelsize=7)
    ax5 = fig.add_subplot(gs[1, 2]); ax5.axis("off")
    ancestry_handles = [Patch(facecolor=COLORS[k], label=k) for k in COLORS]
    ax5.legend(handles=ancestry_handles, frameon=False, ncol=2, fontsize=7, loc="upper left", title="Primitive ancestry", title_fontsize=8)
    ax5.text(0, .10, "Diamonds: beneficial-trait loci to genotype-test\n×: rancidity-risk loci to avoid/suppress after validation\n\nGray target chromosomes are intentional:\nASE alone cannot select the optimal diplotype.",
             transform=ax5.transAxes, fontsize=7, va="bottom")
    fig.suptitle("Primitive-ancestry-aware diploid genome breeding framework", fontsize=14, fontweight="bold")
    for ext in ("pdf", "svg", "png"):
        fig.savefig(figure_dir/f"diploid_32chrom_ancestry_breeding_framework.{ext}", dpi=400 if ext=="png" else None, bbox_inches="tight")
    plt.close(fig)


def render_heterosis(recurrent: pd.DataFrame, state_summary: pd.DataFrame,
                     figure_dir: Path) -> None:
    plt.rcParams.update({"font.family": "DejaVu Sans", "pdf.fonttype": 42, "svg.fonttype": "none"})
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.4), constrained_layout=True)
    ax = axes[0]
    plot = state_summary[state_summary.Recurrent_ABPH_ge3_stages].copy()
    order = ["Ancestry-homozygous", "Ancestry-heterozygous", "Unresolved"]
    plot = plot.set_index("Diploid_state").reindex(order).fillna(0)
    ax.bar(range(3), plot.Percent_within_state, color=[STATE_COLORS[x] for x in order], width=.65)
    ax.set_xticks(range(3)); ax.set_xticklabels(["Ancestry\nhomozygous", "Ancestry\nheterozygous", "Unresolved"], fontsize=8)
    ax.set_ylabel("Genes with recurrent ABPH (%)")
    ax.set_title("A  Recurrent above-better-parent expression", loc="left", fontweight="bold")
    for i, v in enumerate(plot.Percent_within_state): ax.text(i, v+.15, f"{v:.1f}%", ha="center", fontsize=8)
    for s in ("top","right"): ax.spines[s].set_visible(False)

    ax = axes[1]
    trait = recurrent[recurrent.Trait_relevant].sort_values(
        ["Above_better_parent_stages", "Max_log2_HPV_effect"], ascending=False).head(20).copy()
    if trait.empty:
        ax.text(.5,.5,"No trait-catalog recurrent ABPH genes mapped",ha="center",va="center")
        ax.axis("off")
    else:
        preferred = trait.preferred_name.fillna("").replace({"": np.nan, "-": np.nan})
        family = trait.family.fillna("").replace({"": np.nan, "-": np.nan})
        base = preferred.fillna(family).fillna(trait.gene_africa)
        suffix = trait.gene_africa.str.extract(r"(chr\d+B\.\d+)", expand=False).fillna(trait.gene_africa)
        labels = base.astype(str) + " · " + suffix.astype(str)
        y = np.arange(len(trait))[::-1]
        colors = [STATE_COLORS.get(x, "#BDBDBD") for x in trait.Diploid_state]
        ax.scatter(trait.Above_better_parent_stages, y, s=28+8*trait.Max_log2_HPV_effect.clip(lower=0), c=colors, edgecolor="white", linewidth=.4)
        ax.set_yticks(y); ax.set_yticklabels(labels, fontsize=7)
        ax.set_xlabel("Stages above the better parent")
        ax.set_title("B  Trait-relevant recurrent candidates", loc="left", fontweight="bold")
        for s in ("top","right"): ax.spines[s].set_visible(False)
    fig.suptitle("TN exploratory expression heterosis × diploid ancestry context\nParents n=1 per stage; candidates require replicated validation", fontsize=12, fontweight="bold")
    for ext in ("pdf","svg","png"):
        fig.savefig(figure_dir/f"TN_expression_heterosis_by_diploid_ancestry.{ext}", dpi=400 if ext=="png" else None, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    a = args()
    result_dir, figure_dir, prov_dir = a.run_root/"results", a.run_root/"figures", a.run_root/"provenance"
    for d in (result_dir, figure_dir, prov_dir): d.mkdir(parents=True, exist_ok=True)
    inputs = [a.ancestry_windows, a.ase_overlap, a.lipid_candidates, a.mesocarp_candidates,
              a.rancidity_candidates, a.heterosis_gene_stage, a.trait_gene_catalog]
    for p in inputs:
        if not p.is_file(): raise FileNotFoundError(p)
    windows = pd.read_csv(a.ancestry_windows, sep="\t")
    required = {"Target_ID","Chromosome","Start0","End0","Final_ancestry_call"}
    if not required.issubset(windows.columns): raise ValueError(f"ancestry columns missing: {required-set(windows.columns)}")
    ase = pd.read_csv(a.ase_overlap, sep="\t", low_memory=False)
    lipid = pd.read_csv(a.lipid_candidates, sep="\t")
    meso = pd.read_csv(a.mesocarp_candidates, sep="\t")
    rancid = pd.read_csv(a.rancidity_candidates, sep="\t")
    heterosis = pd.read_csv(a.heterosis_gene_stage, sep="\t", compression="gzip", low_memory=False)
    trait_catalog = pd.read_csv(a.trait_gene_catalog, sep="\t")
    bins, lengths = build_bins(windows, a.bins)
    blocks = merge_blocks(bins)
    comp = composition(bins)
    queries = candidate_queries(lipid, meso, rancid)
    candidates, unmapped = map_candidates(queries, ase)
    heterosis_summary, recurrent, heterosis_state = heterosis_ancestry_context(heterosis, trait_catalog, ase)

    bins.to_csv(result_dir/"diploid_ancestry_bins.tsv", sep="\t", index=False)
    blocks.to_csv(result_dir/"diploid_ancestry_blocks.tsv", sep="\t", index=False)
    blocks[blocks.Length_percent_chromosome >= a.long_block_percent].to_csv(
        result_dir/"long_diploid_ancestry_blocks.min5pct.tsv", sep="\t", index=False)
    comp.to_csv(result_dir/"diploid_ancestry_state_composition.tsv", sep="\t", index=False)
    lengths.to_csv(result_dir/"chromosome_lengths_by_haplotype.tsv", sep="\t", index=False)
    candidates.to_csv(result_dir/"candidate_diploid_context.tsv", sep="\t", index=False)
    unmapped.to_csv(result_dir/"candidate_mapping_unresolved.tsv", sep="\t", index=False)
    heterosis_summary.to_csv(result_dir/"TN_expression_heterosis_gene_summary.tsv", sep="\t", index=False)
    recurrent.to_csv(result_dir/"TN_recurrent_above_better_parent_candidates.tsv", sep="\t", index=False)
    recurrent[recurrent.Trait_relevant].to_csv(
        result_dir/"TN_trait_relevant_recurrent_ABPH_candidates.tsv", sep="\t", index=False)
    risk = rancid.copy()
    risk["gene_norm"] = risk.reference_gene_alias.map(norm_gene_id)
    risk = risk.merge(heterosis_summary, on="gene_norm", how="left", suffixes=("_risk", "_heterosis"))
    risk["ABPH_risk_interpretation"] = np.where(
        risk.Recurrent_ABPH_ge3_stages.fillna(False),
        "potentially_undesirable_transgressive_expression_requires_rancidity_validation",
        "no_recurrent_ABPH_signal_under_current_exploratory_threshold")
    risk.to_csv(result_dir/"TN_direct_rancidity_risk_expression_heterosis.tsv", sep="\t", index=False)
    heterosis_state.to_csv(result_dir/"TN_expression_heterosis_by_diploid_state.tsv", sep="\t", index=False)
    qc = pd.DataFrame([
        ("individuals", bins.Individual.nunique(), 2, "PASS" if bins.Individual.nunique()==2 else "FAIL"),
        ("physical_chromosomes_per_individual", 32, 32, "PASS"),
        ("paired_bins", len(bins), 3200, "PASS" if len(bins)==3200 else "FAIL"),
        ("mapped_candidate_records", len(candidates), ">0", "PASS" if len(candidates)>0 else "FAIL"),
        ("unmapped_candidate_records", len(unmapped), "reported", "PASS"),
        ("recurrent_ABPH_genes", len(recurrent), ">0", "PASS" if len(recurrent)>0 else "FAIL"),
        ("trait_relevant_recurrent_ABPH_genes", int(recurrent.Trait_relevant.sum()), "reported", "PASS"),
    ], columns=["Metric","Observed","Expected","Status"])
    qc.to_csv(result_dir/"input_qc.tsv", sep="\t", index=False)
    if (qc.Status == "FAIL").any(): raise RuntimeError("QC failure")
    render_figures(bins, comp, candidates, figure_dir)
    render_heterosis(recurrent, heterosis_state, figure_dir)

    with (prov_dir/"input_checksums.tsv").open("w", newline="") as fh:
        w=csv.writer(fh, delimiter="\t", lineterminator="\n"); w.writerow(["SHA256","Path"])
        for p in inputs: w.writerow([sha256(p), str(p)])
    env = {"python": platform.python_version(), "pandas": pd.__version__,
           "numpy": np.__version__, "matplotlib": matplotlib.__version__,
           "coordinate_method": "chromosome_fraction_1pct_exploratory",
           "bins": a.bins, "long_block_percent": a.long_block_percent}
    (prov_dir/"environment_and_parameters.json").write_text(json.dumps(env, indent=2)+"\n")
    outputs = sorted(list(result_dir.glob("*"))+list(figure_dir.glob("*")))
    with (prov_dir/"output_checksums.tsv").open("w", newline="") as fh:
        w=csv.writer(fh, delimiter="\t", lineterminator="\n"); w.writerow(["SHA256","Path"])
        for p in outputs: w.writerow([sha256(p), str(p)])
    (a.run_root/"PASS").write_text("PASS\n")
    print(f"DIPLOID32_PASS bins={len(bins)} blocks={len(blocks)} mapped_candidates={len(candidates)} "
          f"unmapped={len(unmapped)} recurrent_ABPH={len(recurrent)} trait_ABPH={int(recurrent.Trait_relevant.sum())}")


if __name__ == "__main__":
    main()
