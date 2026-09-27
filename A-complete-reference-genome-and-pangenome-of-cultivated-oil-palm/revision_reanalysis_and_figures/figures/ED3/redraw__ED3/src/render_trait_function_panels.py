#!/usr/bin/env python3
"""Render the comprehensive trait-family atlas and FL functional ASE panel."""

from pathlib import Path
import math

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch, Rectangle
import numpy as np
import pandas as pd


OUT = Path(__file__).resolve().parent
ASE_SOURCE = (
    OUT.parent / "01_ASE" / "00_shared" / "runs" /
    "RUN-ASE-HETEROSIS-DOWNSTREAM-V2-001" / "output" /
    "trait_haplotype_ASE.tsv.gz"
)

FAMILY_GROUPS = [
    ("Fatty-acid synthesis", [
        "ACCase", "ACP", "FabD (MCAT)", "FabG (KAR)", "FabI (ENR)",
        "KASIII", "KAS I/II", "FATA/B", "LACS",
    ]),
    ("Unsaturation and elongation", [
        "SAD", "FAD2", "FAD3", "FAD6", "FAD7/FAD8", "KCS",
    ]),
    ("TAG assembly and storage", [
        "GPAT", "LPAT", "PAP", "DGAT", "PDAT", "PDCT",
    ]),
    ("Rancidity and antioxidant defence", [
        "LOX9", "VTE1", "APX7",
    ]),
]
FAMILIES = [family for _, families in FAMILY_GROUPS for family in families]
GROUP_OF = {family: group for group, families in FAMILY_GROUPS for family in families}
GROUP_SHORT = {
    "Fatty-acid synthesis": "FA SYNTHESIS",
    "Unsaturation and elongation": "UNSATURATION",
    "TAG assembly and storage": "TAG / STORAGE",
    "Rancidity and antioxidant defence": "OXIDATION",
}
FAMILY_ROLE = {
    "ACCase": "malonyl-CoA supply", "ACP": "acyl carrier",
    "FabD (MCAT)": "malonyl transfer", "FabG (KAR)": "chain reduction",
    "FabI (ENR)": "chain reduction", "KASIII": "FA initiation",
    "KAS I/II": "FA elongation", "FATA/B": "acyl termination",
    "LACS": "acyl activation", "SAD": "stearate desaturation",
    "FAD2": "oleate to linoleate", "FAD3": "linolenate formation",
    "FAD6": "plastid desaturation", "FAD7/FAD8": "plastid omega-3 desaturation",
    "KCS": "very-long-chain FA", "GPAT": "glycerolipid entry",
    "LPAT": "LPA acylation", "PAP": "DAG formation",
    "DGAT": "TAG synthesis", "PDAT": "TAG synthesis",
    "PDCT": "PC-DAG exchange", "LOX9": "lipid peroxidation",
    "VTE1": "tocopherol synthesis", "APX7": "ROS detoxification",
}
STAGE_GROUPS = ["Days 0–65", "Days 80–140", "Days 155–185", "Hours 12–72"]

COLORS = {
    "A recurrent": "#E76F61",
    "B recurrent": "#4C9ED1",
    "Switching recurrent": "#F3C85B",
    "No/stage-limited ASE": "#D2D5D7",
}
FL_A = "Africa hap2"
FL_B = "American hap1"


def setup() -> None:
    mpl.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 8.2,
        "axes.linewidth": 0.9, "axes.labelsize": 9.5,
        "axes.titlesize": 11.5, "xtick.labelsize": 7.8,
        "ytick.labelsize": 7.2, "pdf.fonttype": 42,
        "ps.fonttype": 42, "svg.fonttype": "none",
        "savefig.facecolor": "white",
    })


def save(fig: plt.Figure, stem: str) -> None:
    for suffix in ("pdf", "svg", "png"):
        kwargs = {"bbox_inches": "tight", "facecolor": "white"}
        if suffix == "png":
            kwargs["dpi"] = 600
        fig.savefig(OUT / f"{stem}.{suffix}", **kwargs)
    plt.close(fig)


def load_unique_gene_stage() -> pd.DataFrame:
    data = pd.read_csv(ASE_SOURCE, sep="\t", low_memory=False)
    data = data[data["family"].isin(FAMILIES) & data["eligible"].fillna(False)].copy()
    # The upstream catalogue can assign one gene to more than one trait module.
    # ASE statistics are identical, so count a gene-stage observation once.
    data = data.sort_values(["analysis", "family", "gene_id", "stage", "trait_module"])
    data = data.drop_duplicates(["analysis", "family", "gene_id", "stage"], keep="first")
    return data


def classify_genes(data: pd.DataFrame) -> pd.DataFrame:
    records = []
    for (analysis, family, gene_id), frame in data.groupby(
        ["analysis", "family", "gene_id"], observed=True
    ):
        robust = frame[frame["robust_ase"]]
        tested = int(frame["stage"].nunique())
        robust_n = int(robust["stage"].nunique())
        a_n = int((robust["ase_call"] == "Allele_A_biased").sum())
        b_n = int((robust["ase_call"] == "Allele_B_biased").sum())
        threshold = max(2, int(math.ceil(tested / 2)))
        recurrent = robust_n >= threshold
        share = max(a_n, b_n) / robust_n if robust_n else np.nan
        if not recurrent:
            status = "No/stage-limited ASE"
        elif share < 0.75:
            status = "Switching recurrent"
        elif a_n > b_n:
            status = "A recurrent"
        else:
            status = "B recurrent"
        preferred = str(frame["preferred_name"].dropna().iloc[0]) if frame["preferred_name"].notna().any() else ""
        description = str(frame["description"].dropna().iloc[0]) if frame["description"].notna().any() else ""
        records.append({
            "analysis": analysis, "pathway_group": GROUP_OF[family],
            "family": family, "gene_id": gene_id, "preferred_name": preferred,
            "description": description, "tested_stages": tested,
            "robust_stages": robust_n, "allele_A_robust_stages": a_n,
            "allele_B_robust_stages": b_n, "recurrent_threshold": threshold,
            "direction_share": share, "ASE_status": status,
            "allele_A_label": frame["allele_A_label"].dropna().iloc[0],
            "allele_B_label": frame["allele_B_label"].dropna().iloc[0],
        })
    result = pd.DataFrame(records)
    result["family_order"] = result["family"].map({x: i for i, x in enumerate(FAMILIES)})
    result = result.sort_values(
        ["analysis", "family_order", "ASE_status", "robust_stages", "gene_id"],
        ascending=[True, True, True, False, True],
    ).reset_index(drop=True)
    labels = []
    for (_, family), frame in result.groupby(["analysis", "family"], sort=False):
        counters = {}
        for row in frame.itertuples():
            base = row.preferred_name.strip()
            if not base or base in {"-", "nan", "None"}:
                base = family.replace(" (MCAT)", "").replace(" (KAR)", "").replace(" (ENR)", "")
            base = base.replace("/", "")[:9]
            counters[base] = counters.get(base, 0) + 1
            labels.append(base if counters[base] == 1 else f"{base}{counters[base]}")
    result["display"] = labels
    return result.drop(columns="family_order")


def render_02d(atlas: pd.DataFrame) -> None:
    # Preserve the historical source name while adding an explicitly complete atlas.
    atlas.to_csv(OUT / "source_02D_trait_gene_atlas_all.tsv", sep="\t", index=False)
    atlas.to_csv(OUT / "source2_D_gene_tiles_reference.tsv", sep="\t", index=False)
    max_genes = int(atlas.groupby(["analysis", "family"]).size().max())
    nrows = len(FAMILIES) * 2
    fig, ax = plt.subplots(figsize=(12.4, 11.4))
    fig.subplots_adjust(left=0.185, right=0.91, bottom=0.055, top=0.91)
    ax.set_xlim(-1.5, max_genes + 1.15)
    ax.set_ylim(-1.0, nrows + 0.2)
    ax.invert_yaxis()
    ax.axis("off")

    tile_w = 0.86
    for fi, family in enumerate(FAMILIES):
        y0 = fi * 2
        if fi % 2:
            ax.add_patch(Rectangle(
                (-1.45, y0 - 0.48), max_genes + 2.45, 1.92,
                facecolor="#F7F8F8", edgecolor="none", zorder=0,
            ))
        ax.text(-1.62, y0 + 0.45, family, ha="right", va="center",
                fontsize=7.2, fontweight="bold" if family in {"FAD2", "DGAT", "LOX9", "VTE1"} else "normal")
        ax.text(-1.62, y0 + 0.85, FAMILY_ROLE[family], ha="right", va="center",
                fontsize=5.1, color="#6B7478")
        for offset, analysis in enumerate(["FL", "TN"]):
            y = y0 + offset
            subset = atlas[(atlas["analysis"] == analysis) & (atlas["family"] == family)].copy()
            subset["status_order"] = subset["ASE_status"].map({
                "A recurrent": 0, "B recurrent": 1,
                "Switching recurrent": 2, "No/stage-limited ASE": 3,
            })
            subset = subset.sort_values(["status_order", "robust_stages", "gene_id"], ascending=[True, False, True])
            ax.text(-0.72, y, analysis, ha="right", va="center", fontsize=5.8,
                    color="#50585C", fontweight="bold")
            for j, row in enumerate(subset.itertuples()):
                ax.add_patch(Rectangle(
                    (j, y - 0.36), tile_w, 0.72,
                    facecolor=COLORS[row.ASE_status], edgecolor="white", lw=0.55,
                ))
                text_color = "white" if row.ASE_status in {"A recurrent", "B recurrent"} else "#34383A"
                ax.text(j + tile_w / 2, y, str(row.display)[:10], ha="center", va="center",
                        fontsize=4.15, color=text_color)
            ax.text(max_genes + 0.15, y, f"n={len(subset)}", ha="left", va="center",
                    fontsize=5.2, color="#5B6468")
        if fi < len(FAMILIES) - 1:
            ax.plot([-1.45, max_genes + 0.95], [y0 + 1.48, y0 + 1.48],
                    color="#DADFE1", lw=0.45)

    # Pathway strips mirror the pathway blocks in the supplied reference panel.
    start = 0
    for group, families in FAMILY_GROUPS:
        end = start + len(families) * 2
        ax.add_patch(Rectangle(
            (max_genes + 0.72, start - 0.48), 0.62, end - start - 0.04,
            facecolor="#D8DADB", edgecolor="#3F4446", lw=0.65, clip_on=False,
        ))
        ax.text(max_genes + 1.03, (start + end - 1) / 2, GROUP_SHORT[group],
                rotation=-90, ha="center", va="center", fontsize=5.7, fontweight="bold")
        start = end

    ax.text(-1.45, -0.8, "D", fontsize=21, fontweight="bold", va="bottom")
    ax.text(-0.65, -0.72, "Oil-production, unsaturation and rancidity ASE atlas",
            fontsize=13.0, fontweight="bold", va="bottom", color="#273238")
    legend = [
        Patch(facecolor=COLORS["A recurrent"], edgecolor="none", label="Recurrent allele A bias"),
        Patch(facecolor=COLORS["B recurrent"], edgecolor="none", label="Recurrent allele B bias"),
        Patch(facecolor=COLORS["Switching recurrent"], edgecolor="none", label="Recurrent, direction switching"),
        Patch(facecolor=COLORS["No/stage-limited ASE"], edgecolor="none", label="No/stage-limited ASE"),
    ]
    fig.legend(handles=legend, ncol=4, frameon=False, loc="upper center",
               bbox_to_anchor=(0.56, 0.958), fontsize=7.0, columnspacing=1.15, handlelength=1.35)
    fig.text(
        0.55, 0.022,
        "All eligible genes in the named curated families are shown.  "
        "FL: A = Africa hap2, B = American hap1; TN: A = Dura/TK-like, B = Pisifera/NS-like.",
        ha="center", fontsize=6.5, color="#5F696D",
    )
    save(fig, "02D_trait_family_matrix")


def summarise_fl_families(data: pd.DataFrame) -> pd.DataFrame:
    fl = data[data["analysis"] == "FL"].copy()
    records = []
    for family in FAMILIES:
        frame = fl[fl["family"] == family]
        if frame.empty:
            continue
        group_values = []
        group_rows = []
        for group in STAGE_GROUPS:
            part = frame[frame["stage_group"] == group]
            robust = part[part["robust_ase"]]
            value = float(robust["log2_allele_ratio"].median()) if not robust.empty else np.nan
            group_values.append(value)
            group_rows.append(int(len(robust)))
        valid = np.asarray([x for x in group_values if np.isfinite(x)], dtype=float)
        robust_all = frame[frame["robust_ase"]]
        same_direction = bool(len(valid) == 4 and (np.all(valid > 0) or np.all(valid < 0)))
        overall = float(robust_all["log2_allele_ratio"].median()) if not robust_all.empty else np.nan
        records.append({
            "pathway_group": GROUP_OF[family], "family": family,
            "functional_role": FAMILY_ROLE[family],
            "eligible_gene_stage_rows": int(len(frame)),
            "robust_gene_stage_rows": int(len(robust_all)),
            "robust_ASE_percentage": 100 * len(robust_all) / len(frame),
            "robust_median_log2_A_over_B": overall,
            "stage_group_min": float(np.nanmin(valid)) if len(valid) else np.nan,
            "stage_group_max": float(np.nanmax(valid)) if len(valid) else np.nan,
            "stage_groups_with_robust_ASE": int(len(valid)),
            "same_direction_all_four_groups": same_direction,
            "dominant_haplotype": FL_A if overall > 0 else FL_B,
            **{f"median_{group}": value for group, value in zip(STAGE_GROUPS, group_values)},
            **{f"robust_rows_{group}": value for group, value in zip(STAGE_GROUPS, group_rows)},
        })
    result = pd.DataFrame(records)
    result["family_order"] = result["family"].map({x: i for i, x in enumerate(FAMILIES)})
    return result.sort_values("family_order").drop(columns="family_order")


def render_06c(summary: pd.DataFrame) -> None:
    summary.to_csv(OUT / "source_06C_FL_functional_complementarity.tsv", sep="\t", index=False)
    plot = summary[np.isfinite(summary["robust_median_log2_A_over_B"])].copy().reset_index(drop=True)
    plot["display"] = plot["family"] + "  ·  " + plot["functional_role"]
    y = np.arange(len(plot))[::-1]
    fig, ax = plt.subplots(figsize=(9.2, 8.2))
    fig.subplots_adjust(left=0.335, right=0.96, bottom=0.19, top=0.84)

    ax.axvspan(-4.2, 0, color="#EEF5FA", zorder=0)
    ax.axvspan(0, 4.2, color="#FCEFEB", zorder=0)
    for i in range(len(plot)):
        if i % 2:
            ax.axhspan(y[i] - 0.46, y[i] + 0.46, color="#F7F8F8", alpha=0.7, zorder=1)
    ax.axvline(0, color="#30383B", lw=0.9, zorder=2)

    display_clip = 4.0
    for yi, row in zip(y, plot.itertuples()):
        colour = "#3E86B5" if row.robust_median_log2_A_over_B < 0 else "#D86656"
        range_min = float(np.clip(row.stage_group_min, -display_clip, display_clip))
        range_max = float(np.clip(row.stage_group_max, -display_clip, display_clip))
        point_x = float(np.clip(row.robust_median_log2_A_over_B, -display_clip, display_clip))
        ax.plot([range_min, range_max], [yi, yi],
                color="#899397", lw=1.05, zorder=3)
        ax.plot([range_min, range_min], [yi - 0.10, yi + 0.10],
                color="#899397", lw=0.85, zorder=3)
        ax.plot([range_max, range_max], [yi - 0.10, yi + 0.10],
                color="#899397", lw=0.85, zorder=3)
        size = 24 + max(0, row.robust_ASE_percentage - 35) * 2.4
        if row.same_direction_all_four_groups:
            ax.scatter(point_x, yi, s=size, marker="o",
                       facecolor=colour, edgecolor="white", lw=0.8, zorder=5)
        else:
            ax.scatter(point_x, yi, s=size * 0.78, marker="D",
                       facecolor="white", edgecolor=colour, lw=1.25, zorder=5)

    ax.set_yticks(y, plot["display"])
    limit = 4.2
    ax.set_xlim(-limit, limit)
    ax.set_ylim(-0.8, len(plot) - 0.2)
    ax.set_xlabel("Robust ASE median  log$_2$(Africa hap2 / American hap1)", labelpad=9)
    ax.grid(axis="x", color="#DDE2E3", lw=0.55, zorder=1)
    ax.tick_params(axis="y", length=0, pad=8)
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.spines["bottom"].set_linewidth(1.05)

    # Pathway headings on the left make the single axis readable as a functional narrative.
    cursor = 0
    for group, families in FAMILY_GROUPS:
        present = [f for f in families if f in set(plot["family"])]
        if not present:
            continue
        indices = [int(plot.index[plot["family"] == f][0]) for f in present]
        ymid = float(np.mean([y[i] for i in indices]))
        ax.text(-0.50, ymid, GROUP_SHORT[group], transform=ax.get_yaxis_transform(),
                rotation=90, ha="center", va="center", fontsize=6.0,
                color="#677276", fontweight="bold", clip_on=False)
        cursor += len(present)

    ax.text(0.25, 1.025, "American hap1 contribution", transform=ax.transAxes,
            ha="center", va="bottom", color="#3279A8", fontsize=8.3, fontweight="bold")
    ax.text(0.75, 1.025, "Africa hap2 contribution", transform=ax.transAxes,
            ha="center", va="bottom", color="#C95B4D", fontsize=8.3, fontweight="bold")
    fig.text(0.035, 0.97, "C", fontsize=23, fontweight="bold", va="top")
    fig.text(0.105, 0.957, "FL functional haplotype complementarity",
             fontsize=13.2, fontweight="bold", va="top", color="#263238")
    fig.text(
        0.105, 0.915,
        "Filled circles: the same allelic direction in all four stage groups; "
        "open diamonds: stage-dependent direction. Horizontal ranges show stage-group medians (display clipped at |4|).",
        fontsize=7.3, color="#647075", va="top",
    )
    shape_handles = [
        Line2D([0], [0], marker="o", color="none", markerfacecolor="#6B879A",
               markeredgecolor="white", markersize=7, label="4/4 groups concordant"),
        Line2D([0], [0], marker="D", color="none", markerfacecolor="white",
               markeredgecolor="#6B879A", markersize=6, label="Direction switching"),
    ]
    size_handles = [
        plt.scatter([], [], s=24 + max(0, p - 35) * 2.4, color="#C7CDD0",
                    edgecolor="white", label=f"{p}%") for p in (50, 70, 90)
    ]
    first = ax.legend(handles=shape_handles, frameon=False, ncol=2,
                      loc="lower left", bbox_to_anchor=(-0.01, -0.19),
                      fontsize=7.0, columnspacing=1.0)
    ax.add_artist(first)
    ax.legend(handles=size_handles, title="Robust ASE", frameon=False, ncol=3,
              loc="lower left", bbox_to_anchor=(0.48, -0.205), fontsize=6.8,
              title_fontsize=6.8, columnspacing=0.7, handletextpad=0.25)
    fig.text(
        0.53, 0.018,
        "Functional families are curated from the ASE catalogue. GO trends provide exploratory context only; no displayed GO term passed BH FDR < 0.05.",
        ha="center", fontsize=6.4, color="#687277",
    )
    save(fig, "06C_FL_parent_independent_complementarity")


def main() -> None:
    setup()
    data = load_unique_gene_stage()
    atlas = classify_genes(data)
    render_02d(atlas)
    fl_summary = summarise_fl_families(data)
    render_06c(fl_summary)
    stable = fl_summary[fl_summary["same_direction_all_four_groups"]]
    print(f"02D genes: {len(atlas)} across {len(FAMILIES)} families")
    print(f"06C families: {len(fl_summary)}; concordant across four groups: {len(stable)}")
    for direction in [FL_B, FL_A]:
        names = stable.loc[stable["dominant_haplotype"] == direction, "family"].tolist()
        print(f"{direction}: {', '.join(names) if names else 'none'}")


if __name__ == "__main__":
    main()
