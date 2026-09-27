#!/usr/bin/env python3
"""Integrate new ED candidates and build clean FAD2 wet-lab qPCR panels.

Incoming composite artwork is kept outside the formal Extended_Data and
Supplementary directories unless it adds non-redundant evidence. Figure 1
microsynteny loci are retained as selectable candidates, Figure 4's PAV panel is
marked as a main-Fig. 4g duplicate, and the useful Figure 5 dSV composite is split
into independent vector panels without embedded panel letters.
"""

from __future__ import annotations

import shutil
import subprocess
import zipfile
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import wilcoxon


ROOT = Path("${ANALYSIS_DIR}")
FINAL = ROOT / "22_answer_reviews/00_ms/04_最终正文/01_revision_ms/02_supplementary_figures_tables_final"
INCOMING = ROOT / "22_answer_reviews/00_ms/04_最终正文/02_Extended_Data_Figures"
FAD2 = ROOT / "22_answer_reviews/01_FAD2/03_FAD2湿实验设计"

ZIP1 = INCOMING / "Supple-fig1.zip"
ZIP4 = INCOMING / "Supple_fig4.zip"
ZIP5 = INCOMING / "Supple-fig5.zip"
FAD_ZIPS = [FAD2 / f"上交补充实验_{i}.zip" for i in (1, 2, 3)]

F2_SUPP = FINAL / "Figure2/Supplementary"
F2_SRC = FINAL / "Figure2/Source_Data"
F2_TAB = FINAL / "Figure2/Tables"

STAGES = ["0d", "15d", "35d", "50d", "65d", "80d", "95d", "110d", "125d", "140d", "155d", "170d", "185d", "12h", "24h", "36h", "48h", "60h", "72h"]
SAMPLES = ["FL", "TN", "NS", "TK"]
COLORS = {"FL": "#159A9C", "TN": "#E76F51", "NS": "#8E6C9E", "TK": "#6FJOINT2"}


def configure() -> None:
    mpl.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 8, "axes.labelsize": 8,
        "xtick.labelsize": 7, "ytick.labelsize": 7, "legend.fontsize": 7,
        "axes.linewidth": 0.7, "pdf.fonttype": 42, "ps.fonttype": 42,
        "svg.fonttype": "none", "savefig.transparent": True,
    })


def clean_axis(ax: plt.Axes) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(length=3, width=0.6, color="#555555")


def save_vector(fig: plt.Figure, stem: Path) -> None:
    stem.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(stem.with_suffix(".pdf"), bbox_inches="tight", pad_inches=0.04)
    fig.savefig(stem.with_suffix(".svg"), bbox_inches="tight", pad_inches=0.04)
    plt.close(fig)


def zip_member(zpath: Path, suffix: str) -> bytes:
    with zipfile.ZipFile(zpath) as zf:
        hits = [n for n in zf.namelist() if n.endswith(suffix)]
        if len(hits) != 1:
            raise RuntimeError(f"Expected one {suffix} in {zpath}, found {hits}")
        return zf.read(hits[0])


def write_and_svg(data: bytes, pdf: Path) -> None:
    pdf.parent.mkdir(parents=True, exist_ok=True)
    pdf.write_bytes(data)
    svg = pdf.with_suffix(".svg")
    subprocess.run(["pdftocairo", "-svg", str(pdf), str(svg)], check=True)


def integrate_figure1_candidates(rows: list[list[str]]) -> None:
    dest = FINAL / "Figure1/Candidate_Selection_From_02_Extended_Data_Figures/Recommended"
    loci = [
        ("PDAT1_microsyn.pdf", "F1_CAND_01_PDAT1_multihaplotype_microsynteny.pdf"),
        ("Polyketide_cyclase-like_microsynteny.pdf", "F1_CAND_02_polyketide_cyclase_multihaplotype_microsynteny.pdf"),
        ("GDSL_lipase_microsyn.pdf", "F1_CAND_03_GDSL_lipase_multihaplotype_microsynteny.pdf"),
        ("HACD_microsyn.pdf", "F1_CAND_04_HACD_multihaplotype_microsynteny.pdf"),
    ]
    for src, name in loci:
        out = dest / name
        write_and_svg(zip_member(ZIP1, src), out)
        rows.append(["Figure1", name.removesuffix(".pdf"), "candidate_recommended", "Non-redundant locus-level PAV/microsynteny example; select only the minimum needed."])
    tables = FINAL / "Figure1/Candidate_Selection_From_02_Extended_Data_Figures/Tables/Supplementary_Tables_17-22.xlsx"
    tables.parent.mkdir(parents=True, exist_ok=True)
    tables.write_bytes(zip_member(ZIP1, "Supplementary_Tables_17-22.xlsx"))
    rows.append(["Figure1", "FabI-like_microsyn", "excluded_duplicate", "Already represented by formal F1_ED_04."])
    rows.append(["Figure1", "supplefigx_PAV_allhaplotypes", "excluded_composite", "Four-letter composite replaced by independent editable candidate panels."])


def integrate_figure4_duplicate(rows: list[list[str]]) -> None:
    dest = FINAL / "Figure4/Candidate_Selection_From_02_Extended_Data_Figures/Not_Recommended_Main_Figure_Duplicate"
    out = dest / "F4_NOT_RECOMMENDED_pan_gene_PAV_heatmap_duplicates_main_Fig4g.pdf"
    write_and_svg(zip_member(ZIP4, "Supplefigx_pav.pdf"), out)
    rows.append(["Figure4", out.stem, "not_recommended_main_duplicate", "This is the same PAV content already shown in main Fig. 4g."])


def split_dsv_vector(rows: list[list[str]]) -> None:
    dest = FINAL / "Figure5/Candidate_Selection_From_02_Extended_Data_Figures/Recommended"
    dest.mkdir(parents=True, exist_ok=True)
    table_dir = FINAL / "Figure5/Tables"
    by_chr = pd.read_csv(table_dir / "fig5_dSV_dSNP__dsv_by_chromosome.tsv", sep="\t")
    candidates = pd.read_csv(table_dir / "fig5_dSV_dSNP__dsv_v5_candidates.tsv", sep="\t", low_memory=False)
    chroms = [f"chr{i:02d}B" for i in range(1, 17)]

    fig, ax = plt.subplots(figsize=(7.15, 3.1))
    sub = by_chr.set_index("Chrom").reindex(chroms)
    bars = ax.bar(np.arange(16), sub.DSV_Count, color="#4F9DA6", edgecolor="none", width=0.72)
    for bar, value in zip(bars, sub.DSV_Count):
        ax.text(bar.get_x() + bar.get_width()/2, value + max(sub.DSV_Count)*0.018, f"{int(value)}",
                ha="center", va="bottom", fontsize=6.2, color="#4E545B")
    ax.set_xticks(np.arange(16), chroms, rotation=45, ha="right")
    ax.set_ylabel("Derived structural variants")
    ax.set_xlabel("Chromosome")
    clean_axis(ax)
    fig.subplots_adjust(left=0.1, right=0.99, bottom=0.24, top=0.96)
    save_vector(fig, dest / "F5_CAND_01_dSV_count_by_chromosome")

    fig, ax = plt.subplots(figsize=(7.15, 4.4))
    ymap = {c: 15-i for i, c in enumerate(chroms)}
    xmax = by_chr.Chrom_Length_bp.max() / 1e6
    for c in chroms:
        length = float(by_chr.loc[by_chr.Chrom == c, "Chrom_Length_bp"].iloc[0]) / 1e6
        ax.plot([0, length], [ymap[c], ymap[c]], color="#D9DDE1", lw=5.8, solid_capstyle="round", zorder=1)
    type_colors = {"DEL": "#D95F02", "INS": "#3182BD", "INV": "#1B9E77", "DUP": "#756BB1", "TRA": "#E6AB02"}
    for svtype, color in type_colors.items():
        s = candidates[candidates.SVTYPE_Group.astype(str).str.upper() == svtype]
        if s.empty: continue
        ax.scatter(s.Pos / 1e6, s.Chrom.map(ymap), marker="|", s=28, linewidths=0.75,
                   color=color, label=svtype, zorder=2)
    ax.set_yticks([ymap[c] for c in chroms], chroms)
    ax.set_xlim(-2, xmax + 3); ax.set_ylim(-0.7, 15.7)
    ax.set_xlabel("Genomic position (Mb)"); ax.set_ylabel("Chromosome")
    ax.legend(frameon=False, ncol=5, loc="upper center", bbox_to_anchor=(0.5, 1.08), handletextpad=0.3)
    clean_axis(ax)
    fig.subplots_adjust(left=0.12, right=0.99, bottom=0.12, top=0.91)
    save_vector(fig, dest / "F5_CAND_02_dSV_genomic_distribution")
    rows.extend([
        ["Figure5", "F5_CAND_01_dSV_count_by_chromosome", "candidate_recommended", "Chromosome-level dSV counts redrawn from the final source table."],
        ["Figure5", "F5_CAND_02_dSV_genomic_distribution", "candidate_recommended", "Genome-wide derived-SV distribution redrawn from the final candidate catalogue."],
    ])

    tdest = FINAL / "Figure5/Candidate_Selection_From_02_Extended_Data_Figures/Tables"
    tdest.mkdir(parents=True, exist_ok=True)
    for member in [
        "Supplementary_Tables_S1-S6_dSV_dSNP_colocalization.xlsx",
        "Supplementary_Tables_S7-S10_IPHs_design.xlsx",
        "Supplementary_Tables_S11-S13_dSV_genome_distribution.xlsx",
    ]:
        (tdest / member).write_bytes(zip_member(ZIP5, member))
    rows.extend([
        ["Figure5", "Supplefig_Zoom_c1-16", "excluded_duplicate", "Chromosome IPH content already supplied as 16 independent formal vector panels."],
        ["Figure5", "Supplefigx_SV", "excluded_main_overlap", "Most subpanels repeat main Fig. 5b–g; not copied into formal figures."],
    ])


def bh_adjust(pvalues: list[float]) -> np.ndarray:
    p = np.asarray(pvalues, dtype=float)
    n = len(p); order = np.argsort(p); ranked = p[order]
    q = np.minimum.accumulate((ranked * n / np.arange(1, n + 1))[::-1])[::-1]
    out = np.empty(n); out[order] = np.clip(q, 0, 1)
    return out


def fad2_data() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    tidy = pd.read_csv(FAD2 / "tables/qpcr_all3zip_tidy_V2_20260622.tsv", sep="\t")
    tidy["sample"] = tidy["sample_original"].replace({"FL": "FL"})
    tidy["stage_index"] = tidy["timepoint"].map({s: i for i, s in enumerate(STAGES)})
    stage = (tidy.groupby(["phase", "timepoint", "stage_index", "sample"], as_index=False)
             .agg(mean_log2_expression=("log2_relative_expression", "mean"),
                  median_log2_expression=("log2_relative_expression", "median"),
                  assay_count=("assay", "nunique")))
    tests = []
    phase_sets = [
        ("all_19_stages", stage),
        ("fruit_development", stage[stage.phase == "fruit_development"]),
        ("postharvest", stage[stage.phase == "postharvest"]),
    ]
    for phase, sub in phase_sets:
        wide = sub.pivot(index="timepoint", columns="sample", values="mean_log2_expression")
        phase_rows = []
        for other in ["TN", "NS", "TK"]:
            paired = pd.concat([wide["FL"], wide[other]], axis=1, keys=["FL", other]).dropna()
            diff = paired["FL"] - paired[other]
            stat, pval = wilcoxon(diff, alternative="two-sided", zero_method="wilcox")
            phase_rows.append({
                "phase": phase, "contrast": f"FL - {other}", "matched_stages": len(diff),
                "median_paired_difference_log2": diff.median(), "mean_paired_difference_log2": diff.mean(),
                "wilcoxon_W": stat, "P_value_two_sided": pval,
            })
        q = bh_adjust([r["P_value_two_sided"] for r in phase_rows])
        for r, qi in zip(phase_rows, q):
            r["BH_FDR_within_phase"] = qi
            tests.append(r)
    return tidy, stage, pd.DataFrame(tests)


def fad2_timecourse_panel(stage: pd.DataFrame) -> None:
    fig, ax = plt.subplots(figsize=(7.15, 3.25))
    for sample in SAMPLES:
        sub = stage[stage["sample"] == sample].sort_values("stage_index")
        ax.plot(sub.stage_index, sub.mean_log2_expression, lw=1.25, color=COLORS[sample],
                marker="o", ms=3.6, markeredgecolor="white", markeredgewidth=0.35, label=sample)
    ax.axhline(0, color="#A8ADB4", lw=0.7)
    ax.axvline(12.5, color="#A8ADB4", lw=0.8, ls=(0, (2, 2)))
    ax.text(6, 1.02, "Development", transform=ax.get_xaxis_transform(), ha="center", va="bottom", color="#6FJOINT2")
    ax.text(15.5, 1.02, "Post-harvest", transform=ax.get_xaxis_transform(), ha="center", va="bottom", color="#6FJOINT2")
    ax.set_xticks(np.arange(len(STAGES)), STAGES, rotation=45, ha="right")
    ax.set_ylabel("Mean log2 relative FAD2 RNA")
    ax.set_xlabel("Stage")
    ax.legend(frameon=False, ncol=4, loc="upper left")
    clean_axis(ax)
    fig.subplots_adjust(left=0.1, right=0.99, bottom=0.24, top=0.9)
    save_vector(fig, F2_SUPP / "F2_SUPP_14_FAD2_wetlab_qPCR_timecourse")


def sig_label(q: float) -> str:
    if q < 0.001: return "***"
    if q < 0.01: return "**"
    if q < 0.05: return "*"
    return "ns"


def fad2_significance_panel(stage: pd.DataFrame, tests: pd.DataFrame) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(7.15, 3.35), sharey=True)
    phases = [("fruit_development", "Development"), ("postharvest", "Post-harvest")]
    rng = np.random.default_rng(20260806)
    for ax, (phase, label) in zip(axes, phases):
        sub = stage[stage.phase == phase]
        vals = [sub[sub["sample"] == s].mean_log2_expression.to_numpy() for s in SAMPLES]
        bp = ax.boxplot(vals, positions=np.arange(4), widths=0.55, patch_artist=True, showfliers=False,
                        medianprops={"color": "#222222", "lw": 1.0},
                        whiskerprops={"color": "#6FJOINT2", "lw": 0.7},
                        capprops={"color": "#6FJOINT2", "lw": 0.7})
        for box, sample in zip(bp["boxes"], SAMPLES):
            box.set(facecolor=COLORS[sample], alpha=0.28, edgecolor=COLORS[sample], linewidth=0.8)
        for xi, (sample, y) in enumerate(zip(SAMPLES, vals)):
            jitter = rng.uniform(-0.12, 0.12, size=len(y))
            ax.scatter(xi + jitter, y, s=12, color=COLORS[sample], alpha=0.75, edgecolor="white", linewidth=0.25)
        ax.set_xticks(np.arange(4), SAMPLES)
        ax.set_xlabel(label)
        ax.axhline(0, color="#B0B4BA", lw=0.6)
        phase_tests = tests[tests.phase == phase].set_index("contrast")
        y0 = max(max(v) for v in vals) + 0.55
        step = 0.55
        for j, other in enumerate(["TN", "NS", "TK"]):
            x2 = SAMPLES.index(other); y = y0 + j * step
            ax.plot([0, 0, x2, x2], [y - .10, y, y, y - .10], color="#555B63", lw=0.65)
            q = float(phase_tests.loc[f"FL - {other}", "BH_FDR_within_phase"])
            ax.text((0 + x2) / 2, y + .05, f"{sig_label(q)}  q={q:.3g}", ha="center", va="bottom", fontsize=6.5)
        ax.set_ylim(top=y0 + 3 * step + .25)
        clean_axis(ax)
    axes[0].set_ylabel("Stage-level mean log2 relative FAD2 RNA")
    fig.subplots_adjust(left=0.1, right=0.99, bottom=0.17, top=0.97, wspace=0.22)
    save_vector(fig, F2_SUPP / "F2_SUPP_15_FAD2_qPCR_matched_stage_significance")


def integrate_fad2(rows: list[list[str]]) -> None:
    tidy, stage, tests = fad2_data()
    fad2_timecourse_panel(stage)
    fad2_significance_panel(stage, tests)
    tidy.to_csv(F2_SRC / "F2_SUPP_14_FAD2_qPCR_assay_level_tidy.tsv", sep="\t", index=False)
    stage.to_csv(F2_SRC / "F2_SUPP_14_FAD2_qPCR_stage_summary.tsv", sep="\t", index=False)
    tests.to_csv(F2_SRC / "F2_SUPP_15_FAD2_qPCR_matched_stage_tests.tsv", sep="\t", index=False)
    tests.to_csv(F2_TAB / "Table_Fig2_32_FAD2_qPCR_matched_stage_statistics.tsv", sep="\t", index=False)
    shutil.copy2(FAD2 / "tables/fad2_qpcr_primer_reliability_summary.tsv", F2_TAB / "Table_Fig2_33_FAD2_qPCR_primer_reliability.tsv")
    shutil.copy2(FAD2 / "tables/qpcr_all3zip_submission_inventory_V2_20260622.tsv", F2_TAB / "Table_Fig2_34_FAD2_qPCR_submission_inventory.tsv")
    zdir = F2_SRC / "FAD2_wetlab_submission_archives"
    zdir.mkdir(parents=True, exist_ok=True)
    for z in FAD_ZIPS:
        shutil.copy2(z, zdir / z.name)
    rows.extend([
        ["Figure2", "F2_SUPP_14_FAD2_wetlab_qPCR_timecourse", "formal_supplementary", "Wet-lab FAD2 RNA trajectory from the final 19-stage workbook."],
        ["Figure2", "F2_SUPP_15_FAD2_qPCR_matched_stage_significance", "formal_supplementary", "Stage-paired two-sided Wilcoxon tests with within-phase BH correction."],
    ])


def write_index(rows: list[list[str]]) -> None:
    manifest = pd.DataFrame(rows, columns=["figure", "artifact", "status", "reason"])
    manifest.to_csv(FINAL / "CANDIDATE_SELECTION_MANIFEST.tsv", sep="\t", index=False)
    text = """# Figure 1–5 supplementary selection index

The five top-level `Figure1`–`Figure5` directories are the selection units. Formal
recommendations remain in `Extended_Data` and `Supplementary`; newly imported
alternatives are under `Candidate_Selection_From_02_Extended_Data_Figures`.

- Figure 1: four non-redundant microsynteny candidates; FabI is excluded because it already exists as F1_ED_04.
- Figure 2: two formal FAD2 wet-lab qPCR panels plus assay-level source data, corrected tests and all three submission archives.
- Figure 3: no new files in the incoming ED archives; the completed ASE/protein support remains unchanged.
- Figure 4: the imported PAV heatmap is retained only in a not-recommended directory because it duplicates main Fig. 4g.
- Figure 5: the useful dSV overview is split into two independent vector candidates; main-figure-overlapping SV and 16-chromosome composites are excluded.

All newly generated formal panels and recommended candidates are PDF+SVG and have
no preset panel letters. Statistical interpretation limits are recorded in
`CLAIM_BOUNDARIES.md`.
"""
    (FINAL / "CANDIDATE_SELECTION_INDEX.md").write_text(text, encoding="utf-8")


def main() -> None:
    configure()
    rows: list[list[str]] = []
    integrate_figure1_candidates(rows)
    integrate_figure4_duplicate(rows)
    split_dsv_vector(rows)
    integrate_fad2(rows)
    write_index(rows)
    print("Integrated Figure 1/4/5 candidates and built two formal FAD2 qPCR panels with corrected tests.")


if __name__ == "__main__":
    main()
