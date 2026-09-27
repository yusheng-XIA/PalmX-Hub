#!/usr/bin/env python3
"""Render a favorable-only D/P/M plus actual-haplotype rescue blueprint."""
from __future__ import annotations

import argparse
import math
import os
import re
import sys
from collections import Counter
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Rectangle

sys.path.insert(0, str(Path(__file__).resolve().parent))
from importlib import import_module

baseplot = import_module("05_render_P9")


ANCESTRY_COLORS = baseplot.ANCESTRY_COLORS
GENOTYPE_COLORS = baseplot.GENOTYPE_COLORS
BASE_ACTION_STYLE = {
    "INTRODUCE_OR_TUNE": ("^", "#159D82", "Introduce / tune favorable D/P/M allele"),
    "RETAIN_FL": ("^", "#4F70B5", "Retain favorable FL D/P/M allele"),
    "TIMING_SCREEN": (">", "#E3A018", "Screen favorable timing / dosage allele"),
}
RESCUE_COLOR = "#8A5FB2"
ACTUAL_HAP_COLORS = {
    "TN_h1": "#00897B", "TN_h2": "#6BC5BA",
    "FL_HapA": "#7651A8", "FL_HapB": "#B091CF",
}


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--base-targets", required=True, type=Path)
    p.add_argument("--candidates", required=True, type=Path)
    p.add_argument("--ase", required=True, type=Path)
    p.add_argument("--bins", required=True, type=Path)
    p.add_argument("--output-dir", required=True, type=Path)
    return p.parse_args()


def boolish(x) -> bool:
    return str(x).strip().lower() in {"true", "yes", "1", "t"}


def num(x, default=0.0) -> float:
    try:
        v = float(x)
        return default if math.isnan(v) else v
    except (TypeError, ValueError):
        return default


def reference_gene(pair: str) -> str:
    xs = str(pair).split("|")
    ref = next((x for x in xs if re.search(r"\.chr\d+B\.", x)), xs[0] if xs else "")
    return ref.replace("evm.model.", "evm.TU.")


def select_actual_haplotypes(ase: pd.DataFrame) -> pd.DataFrame:
    d = ase.copy()
    d["gene_id"] = d["Allele_pair"].map(reference_gene)
    d["bias_n"] = pd.to_numeric(d["ASE_bias_toward_this_haplotype_count"], errors="coerce").fillna(0)
    g = d.groupby(["gene_id", "Individual", "Target_ID"], as_index=False)["bias_n"].max()
    rows = []
    for (gid, individual), sub in g.groupby(["gene_id", "Individual"], sort=False):
        sub = sub.sort_values("bias_n", ascending=False)
        top = sub.iloc[0]
        second = float(sub.iloc[1]["bias_n"]) if len(sub) > 1 else 0.0
        selected = str(top["Target_ID"]) if float(top["bias_n"]) >= 3 and float(top["bias_n"]) > second else ""
        rows.append({"gene_id": gid, "Individual": individual, "Actual_haplotype": selected,
                     "Top_bias_stage_n": int(top["bias_n"]), "Second_bias_stage_n": int(second)})
    wide = pd.DataFrame(rows).pivot(index="gene_id", columns="Individual", values="Actual_haplotype").reset_index()
    return wide.rename(columns={"TN": "TN_ASE_hap", "FL": "FL_ASE_hap"})


def build_rescue(candidates: pd.DataFrame, ase: pd.DataFrame):
    tn_resolved = candidates["TN_resolved_genotype"].fillna("").ne("")
    fl_resolved = candidates["FL_resolved_genotype"].fillna("").ne("")
    high = candidates[
        (candidates["breeding_tier"] == "Tier_C") &
        (candidates["breeding_evidence_score"] >= 7) &
        (candidates["molecular_support_class_n"] >= 2) &
        (~(tn_resolved & fl_resolved))
    ].copy()
    high = high.merge(select_actual_haplotypes(ase), on="gene_id", how="left")
    high["TN_ASE_hap"] = high["TN_ASE_hap"].fillna("")
    high["FL_ASE_hap"] = high["FL_ASE_hap"].fillna("")
    needs_tn = high["TN_resolved_genotype"].fillna("").eq("")
    needs_fl = high["FL_resolved_genotype"].fillna("").eq("")
    complete = (~needs_tn | high["TN_ASE_hap"].ne("")) & (~needs_fl | high["FL_ASE_hap"].ne(""))
    rescued = high[complete].copy()
    risk_mask = rescued["module_ids"].fillna("").str.split(";").map(lambda xs: "M09" in xs)
    rescue_risk = rescued[risk_mask].copy()
    favorable = rescued[~risk_mask].copy()

    actual_action, design_label, favorable_source = [], [], []
    for _, r in favorable.iterrows():
        tn_hap, fl_hap = str(r["TN_ASE_hap"]), str(r["FL_ASE_hap"])
        tn_label = tn_hap or f"TN[{r['TN_resolved_genotype']}]"
        fl_label = fl_hap or f"FL[{r['FL_resolved_genotype']}]"
        avg = num(r.get("prot_current_average_TN_FL_log2FC", 0), 0)
        is_timing = str(r.get("primary_module", "")) == "Fruit development and lipid-window timing"
        quality = str(r.get("primary_module", "")) in {
            "Antioxidant and postharvest stability", "Fatty-acid elongation and composition"
        }
        if is_timing:
            action, source = "FAVORABLE_TIMING_HAPLOTYPE_SCREEN", tn_hap or fl_hap
        elif quality and fl_hap and avg <= 0:
            action, source = "RETAIN_FL_ACTUAL_HAPLOTYPE", fl_hap
        elif tn_hap and (avg >= 0 or boolish(r.get("exploratory_parent_F1_ABPH_support", False))):
            action, source = "INTRODUCE_TN_ACTUAL_HAPLOTYPE", tn_hap
        elif fl_hap:
            action, source = "RETAIN_FL_ACTUAL_HAPLOTYPE", fl_hap
        else:
            action, source = "FAVORABLE_ACTUAL_HAPLOTYPE_SCREEN", tn_hap
        actual_action.append(action)
        design_label.append(f"{tn_label} × {fl_label}")
        favorable_source.append(source)
    favorable["rescue_action"] = actual_action
    favorable["actual_haplotype_design"] = design_label
    favorable["favorable_actual_haplotype"] = favorable_source
    favorable = favorable.sort_values(["breeding_evidence_score", "molecular_support_class_n", "Chromosome", "Start0"], ascending=[False, False, True, True]).reset_index(drop=True)
    favorable.insert(0, "Rescue_Code", [f"R{i:03d}" for i in range(1, len(favorable) + 1)])
    return high, favorable, rescue_risk


def label_indices(base: pd.DataFrame, rescue: pd.DataFrame):
    chosen_base = []
    for module in sorted(base["primary_module"].dropna().unique()):
        chosen_base.extend(base[base["primary_module"] == module].sort_values("breeding_evidence_score", ascending=False).head(1).index)
    chosen_base = list(dict.fromkeys(chosen_base))[:12]
    chosen_rescue = list(rescue.sort_values("breeding_evidence_score", ascending=False).head(12).index)
    return set(chosen_base), set(chosen_rescue)


def render(base: pd.DataFrame, rescue: pd.DataFrame, bins: pd.DataFrame, out: Path) -> None:
    mpl.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 8.0, "axes.linewidth": 0.6,
        "pdf.fonttype": 42, "ps.fonttype": 42, "hatch.linewidth": 0.45, "svg.fonttype": "none",
    })
    fl = bins[bins["Individual"] == "FL"].copy()
    for c in ("Hap1_start0", "Hap1_end0", "Hap2_start0", "Hap2_end0"):
        fl[c] = pd.to_numeric(fl[c], errors="coerce")
    max_mb = max(fl["Hap1_end0"].max(), fl["Hap2_end0"].max()) / 1e6
    fig = plt.figure(figsize=(20.5, 17.0), facecolor="white")
    gs = fig.add_gridspec(1, 2, width_ratios=[4.75, 1.5], wspace=0.035)
    ax, side = fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1])
    height, pitch = 0.23, 1.38
    y_centers = {f"chr{i:02d}": (16 - i) * pitch for i in range(1, 17)}
    pair_len = {}
    for i in range(1, 17):
        chrom = f"chr{i:02d}"; rows = fl[fl["Chromosome"] == chrom].sort_values("Bin_index"); y = y_centers[chrom]
        l1 = baseplot.draw_haplotype(ax, rows, 1, y + 0.165, height)
        l2 = baseplot.draw_haplotype(ax, rows, 2, y - 0.165, height)
        pair_len[chrom] = max(l1, l2)
        ax.text(-8.8, y, str(i), ha="right", va="center", fontsize=8, fontweight="bold")
        ax.text(-7.1, y + 0.165, "H1", ha="right", va="center", fontsize=5.7, color="#555555")
        ax.text(-7.1, y - 0.165, "H2", ha="right", va="center", fontsize=5.7, color="#555555")

    base_labels, rescue_labels = label_indices(base, rescue)
    for chrom in y_centers:
        bsub = base[base["Chromosome"] == chrom]
        rsub = rescue[rescue["Chromosome"] == chrom]
        combined = pd.concat([
            bsub.assign(_kind="base", _frac=pd.to_numeric(bsub["chromosome_fraction"], errors="coerce")),
            rsub.assign(_kind="rescue", _frac=pd.to_numeric(rsub["chromosome_fraction"], errors="coerce")),
        ], sort=False)
        lanes = [-10.0] * 6
        lane_map = {}
        for idx, r in combined.sort_values("_frac").iterrows():
            f = float(r["_frac"]); lane = next((j for j, last in enumerate(lanes) if f - last >= 0.026), int(np.argmin(lanes)))
            lanes[lane] = f; lane_map[(r["_kind"], idx)] = lane
        for idx, r in bsub.iterrows():
            f = float(r["chromosome_fraction"]); x = f * pair_len[chrom]
            y = y_centers[chrom] + 0.39 + lane_map[("base", idx)] * 0.105
            marker, color, _ = BASE_ACTION_STYLE[r["breeding_action"]]
            size = 30 if r["breeding_tier"] == "Tier_A" else 17
            ax.scatter(x, y, marker=marker, s=size, facecolor=color, edgecolor="white", linewidth=0.35, zorder=6)
            alleles = str(r["target_DPM_genotype"]).split("/")
            if len(alleles) == 2:
                ax.scatter([x - 0.65, x + 0.65], [y - 0.075, y - 0.075], marker="s", s=7,
                           c=[GENOTYPE_COLORS.get(alleles[0], "#BBB"), GENOTYPE_COLORS.get(alleles[1], "#BBB")],
                           edgecolor="#333", linewidth=0.15, zorder=7)
            if idx in base_labels:
                ax.text(x + 1.0, y, r["Target_Code"], fontsize=5, fontweight="bold", va="center",
                        bbox=dict(boxstyle="round,pad=0.10", facecolor="white", edgecolor="#D8D8D8", lw=0.25, alpha=0.88))
        for idx, r in rsub.iterrows():
            f = float(r["chromosome_fraction"]); x = f * pair_len[chrom]
            y = y_centers[chrom] + 0.39 + lane_map[("rescue", idx)] * 0.105
            ax.scatter(x, y, marker="^", s=34, facecolor=RESCUE_COLOR, edgecolor="#3F2950", linewidth=0.45, zorder=7)
            haps = [h for h in (str(r["TN_ASE_hap"]), str(r["FL_ASE_hap"])) if h]
            if haps:
                xs = np.linspace(x - 0.65, x + 0.65, len(haps)) if len(haps) > 1 else [x]
                ax.scatter(xs, [y - 0.075] * len(haps), marker="s", s=7,
                           c=[ACTUAL_HAP_COLORS[h] for h in haps], edgecolor="#333", linewidth=0.15, zorder=8)
            if idx in rescue_labels:
                ax.text(x + 1.0, y, r["Rescue_Code"], fontsize=5, fontweight="bold", color="#5A3673", va="center",
                        bbox=dict(boxstyle="round,pad=0.10", facecolor="#F7F0FA", edgecolor="#C9B1D5", lw=0.3, alpha=0.92))

    ymin, ymax = min(y_centers.values()) - 0.7, max(y_centers.values()) + 0.9
    ax.set_xlim(-11.5, max_mb * 1.03); ax.set_ylim(ymin, ymax); ax.set_yticks([])
    ax.xaxis.set_ticks_position("top"); ax.xaxis.set_label_position("top")
    ax.set_xlabel("Physical position (Mb)", fontsize=8.5, labelpad=7)
    ax.set_xticks(np.arange(0, math.ceil(max_mb / 25) * 25 + 1, 25)); ax.tick_params(axis="x", labelsize=7, length=3, width=0.6)
    for s in ("left", "right", "bottom"): ax.spines[s].set_visible(False)
    ax.spines["top"].set_color("#777777")
    ax.text(-0.028, 1.03, "a", transform=ax.transAxes, fontsize=14, fontweight="bold", va="bottom")
    ax.text(0.0, 1.035, "FL genomic chassis — favorable-only 32-chromosome breeding blueprint",
            transform=ax.transAxes, fontsize=12.8, fontweight="bold", va="bottom")
    ax.text(0.0, 1.009, "D/P/M-resolved favorable targets plus ASE-resolved TN_h1/TN_h2/FL_HapA/FL_HapB rescue targets; risk markers are excluded",
            transform=ax.transAxes, fontsize=7.4, color="#555555", va="bottom")
    ax.text(-8.8, ymax - 0.2, "Chr.", ha="right", fontsize=7, color="#555555")

    side.axis("off"); side.text(0.0, 1.03, "b", transform=side.transAxes, fontsize=14, fontweight="bold", va="bottom")
    side.text(0.08, 1.032, "Favorable breeding key", transform=side.transAxes, fontsize=11.2, fontweight="bold", va="bottom")
    y = 0.987; side.text(0.0, y, "Primitive ancestry background", fontsize=8.7, fontweight="bold", transform=side.transAxes); y -= 0.028
    for anc in ("Dura", "Pisifera", "Meizhou4"):
        side.add_patch(Rectangle((0, y - 0.009), 0.06, 0.018, transform=side.transAxes, facecolor=ANCESTRY_COLORS[anc], edgecolor="#333", lw=0.35))
        side.text(0.075, y, anc, transform=side.transAxes, va="center", fontsize=7.4); y -= 0.026
    side.add_patch(Rectangle((0, y - 0.009), 0.06, 0.018, transform=side.transAxes, facecolor="white", edgecolor="#999", hatch="////", lw=0.35))
    side.text(0.075, y, "Background ancestry pending", transform=side.transAxes, va="center", fontsize=7.4)
    y -= 0.05; side.text(0.0, y, "Favorable target action", fontsize=8.7, fontweight="bold", transform=side.transAxes); y -= 0.029
    for action, (marker, color, label) in BASE_ACTION_STYLE.items():
        side.scatter([0.025], [y], transform=side.transAxes, marker=marker, s=36, color=color, edgecolor="white", lw=0.4)
        side.text(0.075, y, label, transform=side.transAxes, va="center", fontsize=6.9); y -= 0.03
    side.scatter([0.025], [y], transform=side.transAxes, marker="^", s=39, color=RESCUE_COLOR, edgecolor="#3F2950", lw=0.45)
    side.text(0.075, y, "Favorable actual-haplotype rescue target", transform=side.transAxes, va="center", fontsize=6.9)
    y -= 0.046; side.text(0.0, y, "Actual haplotype chips", fontsize=8.7, fontweight="bold", transform=side.transAxes); y -= 0.027
    for hap in ("TN_h1", "TN_h2", "FL_HapA", "FL_HapB"):
        side.scatter([0.025], [y], transform=side.transAxes, marker="s", s=25, color=ACTUAL_HAP_COLORS[hap], edgecolor="#333", lw=0.25)
        side.text(0.075, y, hap, transform=side.transAxes, va="center", fontsize=7.0); y -= 0.025

    y -= 0.03; side.text(0.0, y, "Design totals", fontsize=8.7, fontweight="bold", transform=side.transAxes); y -= 0.03
    side.text(0.0, y, f"Favorable D/P/M targets: {len(base)}", transform=side.transAxes, fontsize=7.0); y -= 0.024
    side.text(0.0, y, f"Favorable actual-haplotype rescue targets: {len(rescue)}", transform=side.transAxes, fontsize=7.0); y -= 0.024
    side.text(0.0, y, f"Total favorable plotted targets: {len(base)+len(rescue)}", transform=side.transAxes, fontsize=7.0, fontweight="bold"); y -= 0.024
    side.text(0.0, y, "All risk-exclusion markers removed from this design", transform=side.transAxes, fontsize=6.7, color="#A33A33")

    y -= 0.055; side.text(0.0, y, "Priority D/P/M and rescue labels", fontsize=8.7, fontweight="bold", transform=side.transAxes); y -= 0.029
    base_key = base.loc[list(base_labels)].sort_values("Target_Code")
    rescue_key = rescue.loc[list(rescue_labels)].sort_values("Rescue_Code")
    entries = [(r.Target_Code, baseplot.short_family(r.primary_family), r.target_DPM_genotype) for r in base_key.itertuples()]
    entries += [(r.Rescue_Code, baseplot.short_family(r.primary_family), r.actual_haplotype_design) for r in rescue_key.itertuples()]
    for j, (code, fam, target) in enumerate(entries[:24]):
        col, row = divmod(j, 12); x0 = 0 if col == 0 else 0.51; yy = y - row * 0.023
        fam = fam[:10] + "…" if len(fam) > 11 else fam
        target = target.replace("TN[", "T[").replace("FL[", "F[")
        if len(target) > 14: target = target[:13] + "…"
        side.text(x0, yy, f"{code} {fam} {target}", transform=side.transAxes, fontsize=5.3, family="DejaVu Sans Mono")
    side.text(0.0, 0.015, "Only favorable, retained-quality and timing-screen targets are shown.\nActual-haplotype rescue does not claim primitive D/P/M ancestry.",
              transform=side.transAxes, fontsize=6.1, color="#555555", va="bottom")

    for ext in ("pdf", "svg", "png"):
        fig.savefig(out / f"FL_favorable_only_32chrom_DPM_plus_haplotype_rescue.{ext}", dpi=450 if ext == "png" else None,
                    bbox_inches="tight", facecolor="white")
    plt.close(fig)


def main() -> None:
    a = parse_args()
    if a.output_dir.exists(): raise SystemExit(f"[ERROR] output exists: {a.output_dir}")
    stage = a.output_dir.with_name(f"{a.output_dir.name}.building.{os.getpid()}"); stage.mkdir(parents=True)
    base_all = pd.read_csv(a.base_targets, sep="\t", low_memory=False)
    candidates = pd.read_csv(a.candidates, sep="\t", low_memory=False)
    ase = pd.read_csv(a.ase, sep="\t", low_memory=False)
    bins = pd.read_csv(a.bins, sep="\t", low_memory=False)
    base_risk = base_all[base_all["breeding_action"] == "EXCLUDE_TN_RISK"].copy()
    base = base_all[base_all["breeding_action"] != "EXCLUDE_TN_RISK"].copy().reset_index(drop=True)
    high, rescue, rescue_risk = build_rescue(candidates, ase)
    if len(base_all) != 370 or len(base_risk) != 77 or len(base) != 293 or len(high) != 278:
        raise SystemExit(f"[ERROR] unexpected counts base_all={len(base_all)} base={len(base)} risk={len(base_risk)} high={len(high)}")
    if rescue["module_ids"].fillna("").str.contains(r"(^|;)M09(;|$)", regex=True).any():
        raise SystemExit("[ERROR] risk module remains in rescue targets")
    if (base["breeding_action"] == "EXCLUDE_TN_RISK").any():
        raise SystemExit("[ERROR] risk action remains in base targets")
    base_risk.to_csv(stage / "excluded_base_risk_candidates.tsv", sep="\t", index=False)
    rescue_risk.to_csv(stage / "excluded_rescue_risk_candidates.tsv", sep="\t", index=False)
    rescue.to_csv(stage / "favorable_haplotype_rescue_candidates.tsv", sep="\t", index=False)
    base.to_csv(stage / "favorable_DPM_targets.tsv", sep="\t", index=False)
    combined_min = pd.concat([
        base.assign(Plot_Code=base["Target_Code"], Plot_class="DPM_resolved", Plot_target=base["target_DPM_genotype"]),
        rescue.assign(Plot_Code=rescue["Rescue_Code"], Plot_class="actual_haplotype_rescue", Plot_target=rescue["actual_haplotype_design"]),
    ], ignore_index=True, sort=False)
    combined_min.to_csv(stage / "favorable_only_combined_plot_targets.tsv", sep="\t", index=False)
    render(base, rescue, bins, stage)
    summary = [
        "# Favorable-only drawing attempt003", "",
        f"- Original D/P/M targets: {len(base_all)}", f"- Excluded original risk markers: {len(base_risk)}",
        f"- Favorable D/P/M targets retained: {len(base)}", f"- High-evidence ancestry-incomplete candidates reviewed: {len(high)}",
        f"- Excluded rescue risk markers: {len(rescue_risk)}", f"- Favorable actual-haplotype rescue targets added: {len(rescue)}",
        f"- Total favorable targets plotted: {len(combined_min)}", "",
        "Risk markers are retained only in excluded audit tables and are not plotted or proposed for introgression.",
    ]
    (stage / "FAVORABLE_ONLY_SUMMARY.md").write_text("\n".join(summary) + "\n", encoding="utf-8")
    stage.rename(a.output_dir)
    print(f"FAVORABLE_RESCUE_RENDER_PASS base={len(base)} rescue={len(rescue)} total={len(combined_min)} excluded_risk={len(base_risk)+len(rescue_risk)}")


if __name__ == "__main__":
    main()
