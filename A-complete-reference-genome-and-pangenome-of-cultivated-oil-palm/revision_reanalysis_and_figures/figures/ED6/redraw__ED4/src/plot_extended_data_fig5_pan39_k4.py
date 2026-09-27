#!/usr/bin/env python3
"""Build the PAN39/K4 review candidate for Extended Data Fig. 5."""

from __future__ import annotations

import csv
import hashlib
import io
import json
import math
import sys
from pathlib import Path

import Bio
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from Bio import Phylo
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.patches import Patch, Rectangle


RUN = Path(__file__).resolve().parents[1]
ATTEMPT = RUN / "attempt4_citation_order"
RESULTS = ATTEMPT / "results"
PANELS = RESULTS / "panels"
QC = ATTEMPT / "qc"
PROVENANCE = ATTEMPT / "provenance"

TREE = Path("${ANALYSIS_DIR}/14_pan_genome/11_new_pan/runs/RUN-PAN39-OF315-20260811-001/results/orthofinder_recovery_attempt2/Results_Aug12/Species_Tree/SpeciesTree_rooted.txt")
MANIFEST = Path("${ANALYSIS_DIR}/14_pan_genome/11_new_pan/config/sample_manifest.tsv")
ASSEMBLY_K4 = Path("${ANALYSIS_DIR}/22_answer_reviews/00_ms/03_V3/04_figure4/assembly34_material_class_scheme_20260709/assembly34_recommended_material_classes.tsv")
K4_ASSIGNMENTS = Path("${ANALYSIS_DIR}/22_answer_reviews/00_ms/03_V3/04_figure4/K3_K4_pi_fst_pca_groups/metadata/K4_dominant_assignments.tsv")
LD_POINTS = Path("${ANALYSIS_DIR}/22_answer_reviews/00_ms/03_V3/04_figure4/K4_LD_decay_20260709/quick_thin2kb/tables_smooth/K4_LD_decay_binned_smooth_points.tsv")
LD50 = Path("${ANALYSIS_DIR}/22_answer_reviews/00_ms/03_V3/04_figure4/K4_LD_decay_20260709/quick_thin2kb/tables_smooth/K4_LD50_summary_smooth.tsv")
PI = Path("${ANALYSIS_DIR}/22_answer_reviews/00_ms/03_V3/04_figure4/K3_K4_pi_fst_pca_groups/summary/K4_PI_weighted.tsv")
FST = Path("${ANALYSIS_DIR}/22_answer_reviews/00_ms/03_V3/04_figure4/K3_K4_pi_fst_pca_groups/summary/K4_FST_weighted_matrix.tsv")

K4_ORDER = ["K4_Pop1", "K4_Pop2", "K4_Pop3", "K4_Pop4"]
K4_LABEL = {k: k.replace("K4_Pop", "K4-Pop") for k in K4_ORDER}
K4_COLORS = {
    "K4_Pop1": "#70B4E7",
    "K4_Pop2": "#F3766A",
    "K4_Pop3": "#82C7B8",
    "K4_Pop4": "#BBA8D8",
}

mpl.rcParams.update(
    {
        "font.family": "DejaVu Sans",
        "font.size": 7,
        "axes.linewidth": 0.65,
        "xtick.major.width": 0.65,
        "ytick.major.width": 0.65,
        "xtick.major.size": 3,
        "ytick.major.size": 3,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "svg.fonttype": "none",
        "savefig.facecolor": "white",
    }
)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def display_name(sample: str) -> str:
    explicit = {
        "Africa_hap2": "FL-Hap2",
        "American_hap1": "FL-Hap1",
        "bk_hap1": "TN-Hap1",
        "bk_hap2": "TN-Hap2",
        "dura_hap1": "TK-Hap1",
        "dura_hap2": "TK-Hap2",
        "pisifera_hap1": "NS-Hap1",
        "pisifera_hap2": "NS-Hap2",
        "nrly_hap1": "Nigerian-Hap1",
        "nrly_hap2": "Nigerian-Hap2",
        "meizhou4_hap1": "E. oleifera-Hap1",
        "meizhou4_hap2": "E. oleifera-Hap2",
        "houke_pa_genome": "Houke",
    }
    if sample in explicit:
        return explicit[sample]
    if sample.endswith("_pa_genome"):
        return f"EG-{sample.removesuffix('_pa_genome')}"
    return sample


def mapping_source(sample: str) -> str | None:
    if sample.startswith("dura_hap"):
        return "EG_dura"
    if sample.startswith("pisifera_hap"):
        return "EG_pisifera"
    if sample.startswith("nrly_hap"):
        return "EG_niriliya"
    if sample == "houke_pa_genome":
        return "EG_houke"
    if sample.endswith("_pa_genome"):
        numeric_id = int(sample.removesuffix("_pa_genome"))
        return f"EG_{numeric_id:03d}"
    if sample.startswith("meizhou4_hap"):
        return None
    return sample


def read_tree_and_mapping():
    tree = Phylo.read(io.StringIO(TREE.read_text().strip()), "newick")
    tips = [tip.name for tip in tree.get_terminals()]
    manifest = pd.read_csv(MANIFEST, sep="\t")
    manifest_samples = set(manifest["sample_id"])
    assert len(tips) == 39 and len(set(tips)) == 39
    assert set(tips) == manifest_samples

    source = pd.read_csv(ASSEMBLY_K4, sep="\t").set_index("Sample")
    rows = []
    for tip in tips:
        source_id = mapping_source(tip)
        row = {
            "tree_tip": tip,
            "display_name": display_name(tip),
            "k4_mapping_source": source_id or "external_E_oleifera",
            "k4_status": "reviewed_material_mapping" if source_id else "external_outside_308_panel",
        }
        if source_id:
            assert source_id in source.index, (tip, source_id)
            s = source.loc[source_id]
            for i in range(1, 5):
                row[f"Q{i}"] = float(s[f"K4_Q{i}"])
            row["dominant_group"] = f"K4_Pop{int(np.argmax([row[f'Q{i}'] for i in range(1, 5)])) + 1}"
            row["max_Q"] = max(row[f"Q{i}"] for i in range(1, 5))
        else:
            for i in range(1, 5):
                row[f"Q{i}"] = np.nan
            row["dominant_group"] = "outside_K4_panel"
            row["max_Q"] = np.nan
        rows.append(row)
    mapping = pd.DataFrame(rows)
    assert mapping["k4_status"].eq("reviewed_material_mapping").sum() == 37
    assert mapping["k4_status"].eq("external_outside_308_panel").sum() == 2
    mapping.to_csv(RESULTS / "PAN39_tree_tip_K4_mapping.tsv", sep="\t", index=False, na_rep="NA")
    return tree, mapping.set_index("tree_tip")


def tree_coordinates(tree):
    terminals = tree.get_terminals()
    y = {tip: float(len(terminals) - idx) for idx, tip in enumerate(terminals)}
    x = tree.depths()
    if max(x.values()) == 0:
        x = tree.depths(unit_branch_lengths=True)

    def assign_y(clade):
        if clade in y:
            return y[clade]
        child_y = [assign_y(child) for child in clade.clades]
        y[clade] = (min(child_y) + max(child_y)) / 2
        return y[clade]

    assign_y(tree.root)
    return x, y


def draw_tree(ax, tree, mapping):
    x, y = tree_coordinates(tree)
    max_x = max(x.values())
    label_x = max_x + max_x * 0.018
    for clade in tree.find_clades(order="preorder"):
        if clade.clades:
            child_y = [y[c] for c in clade.clades]
            ax.plot([x[clade], x[clade]], [min(child_y), max(child_y)], color="#4A4A4A", lw=0.6)
            for child in clade.clades:
                ax.plot([x[clade], x[child]], [y[child], y[child]], color="#4A4A4A", lw=0.6)
    for tip in tree.get_terminals():
        label = mapping.loc[tip.name, "display_name"]
        style = "italic" if label.startswith("E. oleifera") else "normal"
        ax.text(label_x, y[tip], label, va="center", ha="left", fontsize=5.8, fontstyle=style)
    scale = 0.005
    scale_x0 = max_x * 0.02
    scale_y = 0.2
    ax.plot([scale_x0, scale_x0 + scale], [scale_y, scale_y], color="black", lw=0.8, clip_on=False)
    ax.text(scale_x0 + scale / 2, scale_y - 0.55, "0.005 substitutions per site", ha="center", va="top", fontsize=5.5)
    ax.set_xlim(-max_x * 0.015, max_x * 1.52)
    ax.set_ylim(-0.9, len(tree.get_terminals()) + 1.6)
    ax.axis("off")


def draw_qbars(ax, tree, mapping):
    terminals = tree.get_terminals()
    _, y = tree_coordinates(tree)
    for tip in terminals:
        yy = y[tip]
        record = mapping.loc[tip.name]
        if record["k4_status"] == "external_outside_308_panel":
            ax.barh(yy, 1, height=0.68, color="#F2F2F2", edgecolor="#666666", linewidth=0.5)
            ax.text(0.5, yy, "outside K4", ha="center", va="center", fontsize=4.4, color="#444444")
            continue
        left = 0.0
        for idx, group in enumerate(K4_ORDER, start=1):
            value = float(record[f"Q{idx}"])
            ax.barh(yy, value, left=left, height=0.68, color=K4_COLORS[group], edgecolor="none")
            left += value
    ax.set_xlim(0, 1)
    ax.set_ylim(-0.9, len(terminals) + 1.6)
    ax.set_yticks([])
    ax.set_xticks([0, 0.5, 1.0])
    ax.set_xticklabels(["0", "0.5", "1"])
    ax.tick_params(axis="x", labelsize=5.5, pad=1)
    ax.spines[["left", "right", "top"]].set_visible(False)
    ax.set_xlabel("Ancestry proportion", fontsize=6, labelpad=1)
    ax.set_title("K = 4", fontsize=7, pad=3)


def draw_panel_tree(tree_ax, q_ax, tree, mapping, add_label=True, add_legend=True, panel_label="c"):
    draw_tree(tree_ax, tree, mapping)
    draw_qbars(q_ax, tree, mapping)
    if add_label:
        tree_ax.text(-0.04, 1.01, panel_label, transform=tree_ax.transAxes,
                     fontsize=11, fontweight="bold", va="bottom")
    handles = [Patch(facecolor=K4_COLORS[g], label=K4_LABEL[g]) for g in K4_ORDER]
    handles.append(Patch(facecolor="#F2F2F2", edgecolor="#666666",
                         label=r"$E.\ oleifera$ (outside K4 panel)"))
    if add_legend:
        tree_ax.legend(handles=handles, loc="lower left", bbox_to_anchor=(0.0, 1.005), ncol=5,
                       frameon=False, handlelength=1.2, columnspacing=0.8, fontsize=6)


def draw_panel_ld(ax, ld, ld50, add_label=True, panel_label="a"):
    for group in K4_ORDER:
        data = ld.loc[ld["group_id"] == group].sort_values("dist_kb")
        info = ld50.loc[ld50["group_id"] == group].iloc[0]
        ld50_value = str(info["LD50_kb_smooth"])
        ld50_text = f"LD50{ld50_value}" if ld50_value.startswith(">") else f"LD50={ld50_value}"
        label = f"{K4_LABEL[group]} (n={int(info['sample_n'])}; {ld50_text} kb)"
        ax.plot(data["dist_kb"], data["r2_smooth"], color=K4_COLORS[group], lw=1.35, label=label)
    ax.set_xlim(0, 500)
    ax.set_ylim(bottom=0)
    ax.set_xlabel("Distance (kb)")
    ax.set_ylabel("Mean $r^2$")
    ax.grid(axis="y", color="#D8D8D8", lw=0.45, alpha=0.8)
    ax.legend(frameon=False, fontsize=5.5, loc="upper right", handlelength=2.0)
    if add_label:
        ax.text(-0.16, 1.04, panel_label, transform=ax.transAxes,
                fontsize=11, fontweight="bold", va="bottom")


def draw_panel_diversity(pi_ax, fst_ax, pi, fst, add_label=True, panel_label="b"):
    values = pi.set_index("group").loc[K4_ORDER, "PI_weighted"]
    x = np.arange(4)
    pi_ax.bar(x, values.values, color=[K4_COLORS[g] for g in K4_ORDER], width=0.72)
    pi_ax.set_xticks(x, ["Pop1", "Pop2", "Pop3", "Pop4"], rotation=30, ha="right")
    pi_ax.set_ylabel(r"Weighted nucleotide diversity ($\pi$)")
    pi_ax.set_ylim(0, max(values) * 1.22)
    pi_ax.ticklabel_format(axis="y", style="sci", scilimits=(-3, -3))
    pi_ax.grid(axis="y", color="#D8D8D8", lw=0.4)
    if add_label:
        pi_ax.text(-0.35, 1.04, panel_label, transform=pi_ax.transAxes,
                   fontsize=11, fontweight="bold", va="bottom")

    matrix = fst.loc[K4_ORDER, K4_ORDER].astype(float)
    cmap = LinearSegmentedColormap.from_list("fst", ["#F7FBFF", "#9ECAE1", "#3F6C8E"])
    norm = Normalize(vmin=0.04, vmax=0.13)
    for i in range(4):
        for j in range(4):
            facecolor = "#F2F2F2" if i == j else cmap(norm(matrix.iloc[i, j]))
            fst_ax.add_patch(Rectangle((j - 0.5, i - 0.5), 1, 1,
                                       facecolor=facecolor, edgecolor="white", linewidth=0.7))
    fst_ax.set_xlim(-0.5, 3.5)
    fst_ax.set_ylim(3.5, -0.5)
    fst_ax.set_aspect("equal")
    fst_ax.set_xticks(x, ["Pop1", "Pop2", "Pop3", "Pop4"], rotation=30, ha="right")
    fst_ax.set_yticks(x, ["Pop1", "Pop2", "Pop3", "Pop4"])
    for i in range(4):
        for j in range(4):
            text = "-" if i == j else f"{matrix.iloc[i, j]:.3f}"
            color = "white" if matrix.iloc[i, j] > 0.095 and i != j else "#222222"
            fst_ax.text(j, i, text, ha="center", va="center", fontsize=5.5, color=color)
    fst_ax.set_title("Pairwise weighted $F_{ST}$", fontsize=7, pad=3)


def save_figure(fig, stem: Path, tight=True):
    kwargs = {"bbox_inches": "tight"} if tight else {}
    fig.savefig(stem.with_suffix(".pdf"), **kwargs)
    fig.savefig(stem.with_suffix(".svg"), **kwargs)
    fig.savefig(stem.with_suffix(".png"), dpi=600, **kwargs)


def build_panels(tree, mapping, ld, ld50, pi, fst):
    fig, ax = plt.subplots(figsize=(4.5, 3.3))
    draw_panel_ld(ax, ld, ld50, add_label=False)
    fig.tight_layout()
    save_figure(fig, PANELS / "EDFig5a_K4_LD_decay")
    plt.close(fig)

    fig = plt.figure(figsize=(5.5, 3.3))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.0, 1.25], wspace=0.62)
    draw_panel_diversity(fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1]),
                         pi, fst, add_label=False)
    save_figure(fig, PANELS / "EDFig5b_K4_pi_FST")
    plt.close(fig)

    fig = plt.figure(figsize=(8.0, 9.1))
    gs = fig.add_gridspec(1, 2, width_ratios=[4.8, 1.45], wspace=0.03)
    draw_panel_tree(fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1]),
                    tree, mapping, add_label=False)
    save_figure(fig, PANELS / "EDFig5c_PAN39_tree_K4_ancestry")
    plt.close(fig)

def build_composite(tree, mapping, ld, ld50, pi, fst):
    fig = plt.figure(figsize=(8.27, 11.69))
    outer = fig.add_gridspec(2, 1, height_ratios=[2.85, 7.55], hspace=0.10,
                             left=0.075, right=0.97, top=0.97, bottom=0.06)
    top = outer[0].subgridspec(1, 2, width_ratios=[1.08, 1.0], wspace=0.34)
    draw_panel_ld(fig.add_subplot(top[0, 0]), ld, ld50, panel_label="a")
    diversity = top[0, 1].subgridspec(1, 2, width_ratios=[0.85, 1.25], wspace=0.63)
    draw_panel_diversity(fig.add_subplot(diversity[0, 0]), fig.add_subplot(diversity[0, 1]),
                         pi, fst, panel_label="b")

    bottom = outer[1].subgridspec(1, 2, width_ratios=[4.85, 1.35], wspace=0.025)
    draw_panel_tree(fig.add_subplot(bottom[0, 0]), fig.add_subplot(bottom[0, 1]),
                    tree, mapping, panel_label="c", add_legend=False)

    save_figure(fig, RESULTS / "Extended_Data_Fig_5_PAN39_K4_citation_order_review_candidate",
                tight=False)
    plt.close(fig)


def write_caption():
    caption = (
        "Extended Data Fig. 5 | Phylogenomic and population-genetic context of the updated oil palm pangenome. "
        "a, Linkage-disequilibrium decay, measured as smoothed mean r2 against physical distance, for 308 resequenced "
        "accessions assigned to their dominant K4 component (Pop1, n = 89; Pop2, n = 23; Pop3, n = 153; Pop4, n = 43). "
        "LD50 is the distance at which the smoothed curve reaches half its initial value; Pop2 did not reach LD50 within "
        "500 kb. b, Weighted nucleotide diversity (pi; left) and pairwise weighted FST (right) for the four K4 groups. "
        "c, Rooted species tree inferred by OrthoFinder v3.1.5 from 39 haplotype-resolved proteomes representing "
        "33 biological materials. Branch lengths denote substitutions per site. Stacked bars show the four "
        "ADMIXTURE ancestry components (K = 4) assigned through the reviewed material mapping; the independent "
        "E. oleifera Meizhou4 haplotypes are additional assemblies outside the documented 308-accession K4 panel "
        "and are explicitly shown as external to the K4 analysis."
    )
    (RESULTS / "Extended_Data_Fig_5_caption_EN.txt").write_text(caption + "\n")
    caption_cn = (
        "扩展数据图5 | 更新油棕泛基因组的系统发育与群体遗传背景。"
        "a，308份重测序材料按其占比最高的K4组分分组后的连锁不平衡衰减曲线，纵轴为平滑后的平均r2，"
        "横轴为物理距离（Pop1，n = 89；Pop2，n = 23；Pop3，n = 153；Pop4，n = 43）。"
        "LD50定义为平滑曲线降至初始值一半时的距离；Pop2在500 kb范围内未达到LD50。"
        "b，4个K4群体的加权核苷酸多样性pi（左）及两两加权FST（右）。"
        "c，基于代表33份生物学材料的39套单倍型解析蛋白组，使用OrthoFinder v3.1.5推断的有根物种树；"
        "分支长度表示每个位点的替换数。堆叠条表示通过已审核材料映射获得的4个ADMIXTURE祖源组分（K = 4）。"
        "独立的E. oleifera Meizhou4两个单倍型是额外加入的组装，不属于已有记录的308份材料K4群体，"
        "因此明确标记为K4分析之外的外群。"
    )
    (RESULTS / "Extended_Data_Fig_5_caption_CN.txt").write_text(caption_cn + "\n")


def main():
    for directory in (RESULTS, PANELS, QC, PROVENANCE):
        directory.mkdir(parents=True, exist_ok=True)

    inputs = [TREE, MANIFEST, ASSEMBLY_K4, K4_ASSIGNMENTS, LD_POINTS, LD50, PI, FST]
    missing = [str(path) for path in inputs if not path.is_file()]
    assert not missing, missing

    tree, mapping = read_tree_and_mapping()
    ld = pd.read_csv(LD_POINTS, sep="\t")
    ld50 = pd.read_csv(LD50, sep="\t", dtype={"LD50_kb_smooth": str})
    pi = pd.read_csv(PI, sep="\t")
    fst = pd.read_csv(FST, sep="\t", index_col=0)
    assignments = pd.read_csv(K4_ASSIGNMENTS, sep="\t")

    group_sizes = assignments["dominant_group"].value_counts().reindex(K4_ORDER).to_dict()
    assert group_sizes == {"K4_Pop1": 89, "K4_Pop2": 23, "K4_Pop3": 153, "K4_Pop4": 43}
    assert set(ld["group_id"]) == set(K4_ORDER)
    assert set(pi["group"]) == set(K4_ORDER)
    assert list(fst.index) == K4_ORDER and list(fst.columns) == K4_ORDER

    build_panels(tree, mapping, ld, ld50, pi, fst)
    build_composite(tree, mapping, ld, ld50, pi, fst)
    write_caption()

    with (PROVENANCE / "input_sha256.tsv").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["input", "sha256"])
        for path in inputs:
            writer.writerow([path, sha256(path)])

    with (PROVENANCE / "output_sha256.tsv").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["output", "sha256"])
        for path in sorted(p for p in RESULTS.rglob("*") if p.is_file()):
            writer.writerow([path.relative_to(RUN), sha256(path)])
    software = {
        "python": sys.version.split()[0],
        "biopython": Bio.__version__,
        "matplotlib": mpl.__version__,
        "numpy": np.__version__,
        "pandas": pd.__version__,
    }
    (PROVENANCE / "software_versions.json").write_text(json.dumps(software, indent=2) + "\n")
    (PROVENANCE / "command.txt").write_text(f"{sys.executable} {Path(__file__).resolve()}\n")

    qc = {
        "tree_tips": 39,
        "tree_tips_with_reviewed_K4": 37,
        "tree_tips_external_to_K4": 2,
        "K4_group_sizes": group_sizes,
        "resequenced_accession_total": int(sum(group_sizes.values())),
        "outputs": sorted(str(p.relative_to(RUN)) for p in RESULTS.rglob("*") if p.is_file()),
    }
    (QC / "data_consistency.json").write_text(json.dumps(qc, indent=2) + "\n")


if __name__ == "__main__":
    main()
