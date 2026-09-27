#!/usr/bin/env python3
"""Compact, publication-oriented redraw of Figure 2d/e/f.

Layout decisions (2026-09-19, corrected LATEST114 handoff):

d - narrow vertical strip: one enzyme class per row, the four focal palms
    (Coconut, Oleifera, Dura, Pisifera) as four coloured points joined by a
    short connector.  Identical rows are banded teal, variable rows banded
    pale red with a Delta label.  A compact ancestral +/- /net column sits on
    the right (FabG/FabI ancestral +2 gains are highlighted).

e - unchanged overall structure (four axis trajectories on top, FL / TN / and
    FL-TN delta heatmaps below) but the trajectories keep the shaded
    FL - TN band, the largest divergence is marked with a diamond and a text
    label, and the delta heatmap highlights the single most divergent cell.

f - story-oriented: the upper panel is the FL - TN delta RNA heatmap with the
    peak symbol placed by a black diamond-free marker drawn from LATEST114,
    right-hand columns give FL/TN RNA peaks and Astral protein detection,
    the lower panel gives the z>=1 activity windows with peak dots, and a
    story card on the right summarises the two protein-supported highlights
    (FAD2 170 d FL/TN = 1.74x; OLE16 185 d FL/TN = 53.3x) computed directly
    from the stage-detection table.
"""

from __future__ import annotations

import csv
import hashlib
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle


BASE = Path("${ANALYSIS_DIR}/22_answer_reviews/00_ms")
OUT = BASE / "05_MS/0918_revision"
OUT.mkdir(parents=True, exist_ok=True)

D_TABLE = BASE / "03_V3/02_figure/05_Figure2_evolution_multiomics_panels_flat_20260727/Fig2c_core_pathway_enzyme_copy_number_long.tsv"
D_ANCT = BASE / "03_V3/02_figure/03_Figure2_panels_flat_20260727/Fig2b_core_lipid_copy_number_rigorous_v2_ancestral_changes.tsv"
SRC = BASE / "05_MS/MS_revision_3/Final_20260823_VectorRevision/8月29日/Revised_Submission_Candidates_20260901/source_data"
E_AXIS = SRC / "Figure2e_axis_summary_19stage.tsv"
E_MET = SRC / "Figure2e_representative_metabolites_19stage.tsv"
E_COV = SRC / "Figure2e_axis_coverage.tsv"
F_RNA = BASE / "05_MS/MS_revision_3/Final_20260823_VectorRevision/8月29日/10_Table1_ST8_AI_Cleanup_20260902/ST9_latest_reanalysis_audit/Figure2f_RNA_latest114_joint_source.tsv"
F_MARK = BASE / "05_MS/MS_revision_3/Final_20260823_VectorRevision/8月29日/10_Table1_ST8_AI_Cleanup_20260902/ST9_latest_reanalysis_audit/Figure2f_LATEST114_marker_summary.tsv"
F_PROT = SRC / "Figure2f_Astral114_stage_detection.tsv"

FL_C = "#238B8B"
TN_C = "#D95F4A"
INK = "#1F2933"
MUTED = "#66737D"
GRID = "#D9DEE2"
MODULE_C = {
    "Plastid fatty-acid synthesis": "#C98A2C",
    "FA export and modification": "#E0B25C",
    "VLCFA side branch": "#3E8E7E",
    "ER TAG assembly": "#3E6FA3",
    "PC–DAG exchange and desaturation": "#7FA7C9",
}
CM_DELTA = LinearSegmentedColormap.from_list(
    "fig2_delta", ["#2D6EA5", "#F7F7F5", "#C94D3C"]
)
CM_BLUE = LinearSegmentedColormap.from_list("fig2_blue", ["#FFFFFF", "#2D6EA5"])

STAGES = [
    "0d", "15d", "35d", "50d", "65d", "80d", "95d", "110d", "125d",
    "140d", "155d", "170d", "185d", "12h", "24h", "36h", "48h", "60h", "72h",
]
MARKERS = [
    "ACCase", "KASIII", "ENR", "FATA/B", "LACS", "GPAT", "FAD2-like",
    "DGAT", "OLE16", "LOX (loss risk)",
]
GROUPS = {"ACCase": "Push", "KASIII": "Push", "ENR": "Push", "FATA/B": "Push",
          "LACS": "Pull", "GPAT": "Pull", "FAD2-like": "Pull",
          "DGAT": "Package / protect", "OLE16": "Package / protect", "LOX (loss risk)": "Package / protect"}

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 7.0,
    "axes.linewidth": 0.65,
    "axes.edgecolor": "#33383D",
    "axes.labelcolor": INK,
    "axes.titlecolor": INK,
    "xtick.color": INK,
    "ytick.color": INK,
    "savefig.dpi": 600,
})


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def save_pair(fig: plt.Figure, stem: str) -> None:
    fig.savefig(OUT / f"{stem}.pdf", bbox_inches="tight", facecolor="white")
    fig.savefig(OUT / f"{stem}.png", dpi=600, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def stage_info() -> tuple[list[str], int]:
    rows = read_tsv(E_AXIS)
    order = sorted({int(r["stage_index"]): (r["stage"], r["phase"]) for r in rows}.items())
    stages = [stage for _, (stage, _) in order]
    n_dev = sum(phase == "development" for _, (_, phase) in order)
    if stages != STAGES:
        raise AssertionError(f"Unexpected stage order: {stages}")
    return stages, n_dev


DISPLAY_D = {
    "ACCase–BCCP1/2": "ACCase–BCCP1/2", "ACCase–BC / CAC2": "ACCase–BC",
    "ACCase–CTα / CAC3": "ACCase–CTα", "MCAT / FabD": "MCAT/FabD",
    "KASII / FAB1": "KASII/FAB1", "KAR / FabG": "FabG (KAR)",
    "ENR / FabI (MOD1)": "FabI (ENR)", "LACS9-like": "LACS9-like",
    "KCS / FAE (VLCFA side branch)": "KCS/FAE", "PAH1/2 (PAP)": "PAH1/2",
    "FAD2/6 (omega-6 FAD family)": "FAD2/6", "FAD3/7/8 (omega-3 FAD family)": "FAD3/7/8",
    "PDCT / ROD1": "PDCT",
}


def build_d() -> None:
    rows = read_tsv(D_TABLE)
    order_mod = [
        "Plastid fatty-acid synthesis", "FA export and modification",
        "VLCFA side branch", "ER TAG assembly", "PC–DAG exchange and desaturation",
    ]
    genomes = ["Calamus", "Daemonorops", "Nypa_fruticans", "Phoenix_dactylifera",
               "Areca_catechu", "Cocos_nucifera", "American_hap1", "Dura", "Pisifera"]
    labels = {"Calamus": "Calamus", "Daemonorops": "D. draco", "Nypa_fruticans": "Nypa",
              "Phoenix_dactylifera": "Phoenix", "Areca_catechu": "Areca",
              "Cocos_nucifera": "Coconut", "American_hap1": "Oleifera",
              "Dura": "Dura", "Pisifera": "Pisifera"}
    focal = ["Cocos_nucifera", "American_hap1", "Dura", "Pisifera"]
    counts: dict[tuple[str, str], int] = {}
    module: dict[str, str] = {}
    for row in rows:
        counts[(row["Enzyme"], row["Genome"])] = int(row["Curated_gene_locus_count"])
        module[row["Enzyme"]] = row["Pathway_section"]
    enzymes: list[str] = []
    for mod in order_mod:
        enzymes.extend(dict.fromkeys(r["Enzyme"] for r in rows if r["Pathway_section"] == mod))
    identical = sum(len({counts[(e, g)] for g in focal}) == 1 for e in enzymes)

    anc = {}
    for row in read_tsv(D_ANCT):
        anc[row["Enzyme"]] = (
            int(row["Total_inferred_gains"]),
            int(row["Total_inferred_losses"]),
            int(row["Net_change"]),
        )
    anc_map = {
        "ACCase–BCCP1/2": "ACCase", "ACCase–BC / CAC2": "ACCase", "ACCase–CTα / CAC3": "ACCase",
        "MCAT / FabD": "FabD (MCAT)",
        "KAR / FabG": "FabG (KAR)", "ENR / FabI (MOD1)": "FabI (ENR)",
        "FATA": "FATA/B", "FATB": "FATA/B", "LACS9-like": "LACS",
        "KCS / FAE (VLCFA side branch)": "KCS", "GPAT9": "GPAT", "LPAT2": "LPAT",
        "PAH1/2 (PAP)": "PAP", "DGAT1": "DGAT", "DGAT2": "DGAT", "PDCT / ROD1": "PDCT",
        "FAD2/6 (omega-6 FAD family)": "FAD2", "FAD3/7/8 (omega-3 FAD family)": "FAD3/7/8",
    }
    anc_rows = {}
    shown: set[str] = set()
    for enzyme in enzymes:
        source = anc_map.get(enzyme)
        if source and source in anc and source not in shown:
            anc_rows[enzyme] = anc[source]
            shown.add(source)

    focal_colors = {"Cocos_nucifera": "#4E9B99", "American_hap1": "#3978A7",
                    "Dura": "#C5862D", "Pisifera": "#D95F4A"}
    n_rows = len(enzymes)
    max_count = max(counts[(e, g)] for e in enzymes for g in focal)

    fig = plt.figure(figsize=(3.55, 7.05))
    ax = fig.add_axes([0.285, 0.115, 0.415, 0.735])
    ax2 = fig.add_axes([0.775, 0.115, 0.205, 0.735])

    # Main distribution strip.
    for i, enzyme in enumerate(enzymes):
        values = [counts[(enzyme, genome)] for genome in focal]
        is_same = len(set(values)) == 1
        lo, hi = min(values), max(values)
        ax.add_patch(Rectangle((-3.4, i - 0.46), max_count + 7.2, 0.92,
                               color="#E7F1EF" if is_same else "#FCEFEB", zorder=0, lw=0))
        ax.add_patch(Rectangle((-1.95, i - 0.46), 0.16, 0.92,
                               color=MODULE_C[module[enzyme]], zorder=2, lw=0))
        ax.plot([lo, hi], [i, i], color="#9AA5AA" if is_same else "#C66B5E",
                lw=1.4 if not is_same else 1.0, solid_capstyle="round", zorder=2)
        for genome, value in zip(focal, values):
            ax.plot(value, i, marker="o", ms=3.6, color=focal_colors[genome],
                    mec="white", mew=0.5, zorder=4)
        if not is_same:
            ax.text(hi + 0.65, i, f"Δ{hi - lo}", ha="left", va="center",
                    fontsize=4.9, fontweight="bold", color="#B43A32")

    ax.set_xlim(-3.4, max_count + 2.6)
    ax.set_ylim(n_rows - 0.5, -1.75)
    ax.set_xticks(np.arange(0, max_count + 1, 10))
    ax.set_xlabel("locus copies", fontsize=5.6, labelpad=2.5)
    ax.set_yticks(range(n_rows))
    ax.set_yticklabels([DISPLAY_D.get(e, e.replace(" / ", "/")) for e in enzymes], fontsize=4.8)
    ax.tick_params(length=1.8, pad=1.2, labelsize=4.8)
    ax.grid(axis="x", color=GRID, lw=0.4, zorder=0)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)

    # Compact ancestral column: gains / losses / net.
    ax2.set_xlim(0, 3.0)
    ax2.set_ylim(n_rows - 0.5, -2.6)
    ax2.axis("off")
    ax2.set_xticks([])
    ax2.set_yticks([])
    ax2.text(0.95, -1.28, "anc. +/– /net", fontsize=5.4, color=MUTED, ha="center",
             fontweight="bold")
    legend_g = Line2D([], [], marker="o", color="none", markerfacecolor="#4E9B99",
                      markeredgecolor="white", markersize=4.2, label="Coconut")
    legend_o = Line2D([], [], marker="o", color="none", markerfacecolor="#3978A7",
                      markeredgecolor="white", markersize=4.2, label="Oleifera")
    legend_d = Line2D([], [], marker="o", color="none", markerfacecolor="#C5862D",
                      markeredgecolor="white", markersize=4.2, label="Dura")
    legend_p = Line2D([], [], marker="o", color="none", markerfacecolor="#D95F4A",
                      markeredgecolor="white", markersize=4.2, label="Pisifera")
    for i, enzyme in enumerate(enzymes):
        if enzyme not in anc_rows:
            continue
        gain, loss, net = anc_rows[enzyme]
        g_hi = (enzyme in ("KAR / FabG", "ENR / FabI (MOD1)")) and gain > 0
        ax2.text(0.65, i, f"+{gain}" if gain else "·",
                 ha="center", va="center", fontsize=5.0,
                 color="#B43A32" if gain else "#A8AFB4",
                 fontweight="bold" if g_hi else "normal")
        ax2.text(1.50, i, f"{loss}" if loss else "·",
                 ha="center", va="center", fontsize=5.0,
                 color="#2D6EA5" if loss else "#A8AFB4")
        ax2.text(2.35, i, f"{net:+d}" if net else "·",
                 ha="center", va="center", fontsize=5.0,
                 color="#B43A32" if net > 0 else ("#2D6EA5" if net < 0 else MUTED))
    ax2.text(0.65, -1.85, "gains", ha="center", va="center", fontsize=4.5, color=MUTED)
    ax2.text(1.50, -1.85, "loss", ha="center", va="center", fontsize=4.5, color=MUTED)
    ax2.text(2.35, -1.85, "net", ha="center", va="center", fontsize=4.5, color=MUTED)

    fig.legend(handles=[legend_g, legend_o, legend_d, legend_p], loc="lower left",
               bbox_to_anchor=(0.150, 0.012), frameon=False, fontsize=4.9,
               ncol=4, handletextpad=0.2, columnspacing=0.5)
    fig.text(0.030, 0.012, "FabG & FabI: +2 ancestral gains",
             fontsize=4.8, color="#B43A32", fontweight="bold")
    fig.text(0.030, 0.965, "d  Core lipid-pathway dosage is conserved",
             fontsize=8.8, fontweight="bold", color=INK, va="top")
    fig.text(0.030, 0.928,
             f"{identical}/24 enzyme classes have identical locus counts across Coconut, "
             "Oleifera, Dura and Pisifera; separated points (red band) show the minority "
             "dosage changes; right column = ancestral copy-number changes.",
             fontsize=5.7, color="#38434A", va="top", wrap=True)
    save_pair(fig, "Figure2d_dosage_compact_LATEST")


def build_e() -> None:
    stages, n_dev = stage_info()
    axis_rows = read_tsv(E_AXIS)
    coverage = {r["axis_id"]: r for r in read_tsv(E_COV)}
    metabolites = read_tsv(E_MET)
    axes_order = ["P02", "P01", "P03", "P04"]
    titles = {"P02": "Storage lipids", "P01": "Oleic balance",
              "P03": "Hydrolytic deterioration", "P04": "Oxidative deterioration"}
    series: dict[tuple[str, str], dict[str, float]] = {}
    all_values = []
    for row in axis_rows:
        series.setdefault((row["axis_id"], row["genotype"]), {})[row["stage"]] = float(row["mean"])
        all_values.append(float(row["mean"]))
    y_min = min(all_values) - 0.04
    y_max = max(all_values) + 0.04
    x = np.arange(len(stages))

    fig = plt.figure(figsize=(10.8, 7.4))
    gs = fig.add_gridspec(2, 4, height_ratios=(1.02, 1.65), hspace=0.55, wspace=0.28,
                          left=0.075, right=0.955, top=0.885, bottom=0.09)
    for k, axis_id in enumerate(axes_order):
        ax = fig.add_subplot(gs[0, k])
        fl = np.array([series[(axis_id, "FL")][stage] for stage in stages])
        tn = np.array([series[(axis_id, "TN")][stage] for stage in stages])
        ax.axvspan(n_dev - 0.5, len(stages) - 0.5, color="#F1F2F2", zorder=0)
        ax.fill_between(x, fl, tn, where=fl >= tn, color=FL_C, alpha=0.24, lw=0)
        ax.fill_between(x, fl, tn, where=fl < tn, color=TN_C, alpha=0.24, lw=0)
        ax.plot(x, fl, color=FL_C, lw=1.35, marker="o", ms=2.4, zorder=3)
        ax.plot(x, tn, color=TN_C, lw=1.35, marker="o", ms=2.4, zorder=3)
        difference = fl - tn
        largest = int(np.argmax(np.abs(difference)))
        winner = "FL" if difference[largest] > 0 else "TN"
        ax.plot(largest, fl[largest] if winner == "FL" else tn[largest], marker="D", ms=3.6,
                color=FL_C if winner == "FL" else TN_C, mec="white", mew=0.55, zorder=5)
        ax.text(0.98, 0.93, f"{winner} higher · {stages[largest]}", transform=ax.transAxes,
                ha="right", va="top", fontsize=5.8, color=FL_C if winner == "FL" else TN_C,
                fontweight="bold")
        ax.axhline(0, color="#B7BEC2", lw=0.55)
        ax.axvline(n_dev - 0.5, color="#7B858B", lw=0.7, ls=(0, (3, 2)))
        ax.set_ylim(y_min, y_max)
        ax.set_xlim(-0.5, len(stages) - 0.5)
        ax.set_title(titles[axis_id], fontsize=8.2, pad=12)
        c = coverage[axis_id]
        ax.text(0.0, 1.02, f"{c['level2_compound_n']} Level-2 + {c['ms1_mass_proxy_compound_n']} MS1 proxies",
                transform=ax.transAxes, fontsize=5.7, color=MUTED)
        ax.set_xticks(x[::2])
        ax.set_xticklabels([stages[i] for i in range(0, len(stages), 2)], rotation=45,
                           ha="right", fontsize=5.6)
        ax.tick_params(axis="y", labelsize=5.8)
        if k == 0:
            ax.set_ylabel("axis score", fontsize=6.8)
        else:
            ax.set_yticklabels([])
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
        if k == 3:
            ax.legend(handles=[Line2D([], [], color=FL_C, marker="o", ms=3, lw=1.2, label="FL"),
                               Line2D([], [], color=TN_C, marker="o", ms=3, lw=1.2, label="TN")],
                      frameon=False, fontsize=6.0, loc="lower right", ncol=2,
                      handlelength=1.4, columnspacing=0.7)

    met_order: list[tuple[str, str]] = []
    for axis_id in axes_order:
        names = list(dict.fromkeys(r["display_name"] for r in metabolites if r["axis_id"] == axis_id))
        met_order.extend((axis_id, name) for name in names)
    mvals = {(r["display_name"], r["genotype"], r["stage"]):
             (float(r["mean_consensus_z"]) if r["mean_consensus_z"] not in ("", "NA") else np.nan)
             for r in metabolites}

    def matrix(genotype: str) -> np.ndarray:
        matrix_ = np.full((len(met_order), len(stages)), np.nan)
        for i, (_, name) in enumerate(met_order):
            for j, stage in enumerate(stages):
                matrix_[i, j] = mvals.get((name, genotype, stage), np.nan)
        return matrix_

    mf, mt = matrix("FL"), matrix("TN")
    md = mf - mt
    sub = gs[1, :].subgridspec(1, 3, wspace=0.10)
    hm_axes = [fig.add_subplot(sub[0, i]) for i in range(3)]
    for ax, mat, title, color in ((hm_axes[0], mf, "FL", FL_C),
                                  (hm_axes[1], mt, "TN", TN_C),
                                  (hm_axes[2], md, f"Δ = FL − TN  highest |Δ|: "
                                                  f"{STAGES[np.unravel_index(np.nanargmax(np.abs(md)), md.shape)[1]]}",
                                   INK)):
        if ax is hm_axes[2]:
            image = ax.imshow(mat, cmap=CM_DELTA, vmin=-1.6, vmax=1.6, aspect="auto")
            i_hi, j_hi = np.unravel_index(np.nanargmax(np.abs(md)), md.shape)
            ax.add_patch(Rectangle((j_hi - 0.5, i_hi - 0.5), 1.0, 1.0, fill=False,
                                   edgecolor="#22282A", lw=1.1, zorder=5))
        else:
            ax.imshow(mat, cmap=CM_BLUE, vmin=0, vmax=2.6, aspect="auto")
        ax.set_title(title, fontsize=7.6, color=color, pad=4)
        ax.axvline(n_dev - 0.5, color="#626C72", lw=0.75, ls=(0, (3, 2)))
        ax.set_xticks(range(len(stages)))
        ax.set_xticklabels(stages, rotation=90, fontsize=5.0)
        ax.set_yticks(range(len(met_order)))
        ax.set_yticklabels([name for _, name in met_order] if ax is hm_axes[0] else [], fontsize=5.9)
        ax.set_xticks(np.arange(-0.5, len(stages), 1), minor=True)
        ax.grid(which="minor", color="white", lw=0.4)
        ax.tick_params(which="minor", length=0)
        for boundary in (1.5, 3.5, 5.5):
            ax.axhline(boundary, color="#7E898F", lw=0.45, zorder=4)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
    cbar = fig.colorbar(image, ax=hm_axes[2], fraction=0.045, pad=0.02)
    cbar.set_label("consensus z difference", fontsize=6.0)
    cbar.ax.tick_params(labelsize=5.5)
    fig.text(0.076, 0.975, "e  Molecular-axis trajectories expose the FL–TN divergence",
             fontsize=9.2, fontweight="bold", color=INK, va="top")
    fig.text(0.076, 0.950,
             "19 measured stages, no interpolation; shaded band = FL − TN direction; "
             "◆ = stage of largest |Δ|.",
             fontsize=6.6, color="#38434A", va="top")
    fig.text(0.076, 0.928,
             "Constructed axes: Level-2 + MS1 proxies (see coverage under each title); "
             "dashed line = transition to post-harvest.",
             fontsize=6.2, color=MUTED, va="top")
    fig.text(0.076, 0.525, "Representative compounds", fontsize=8.0, fontweight="bold",
             color=INK, va="bottom")
    fig.text(0.076, 0.505,
             "13 developmental stages (0–185 d) + 6 post-harvest stages (12–72 h); "
             "TG 42:0: FL 95 d n=1, TN 36 h n=0/NA.",
             fontsize=5.7, color=MUTED, va="bottom")
    save_pair(fig, "Figure2e_axis_delta_LATEST")


def build_f() -> None:
    stages, n_dev = stage_info()
    rna_rows = read_tsv(F_RNA)
    mark_rows = read_tsv(F_MARK)
    prot_rows = read_tsv(F_PROT)
    summary = {(r["marker"], r["genotype"]): r for r in mark_rows}
    z = {(r["marker"], r["genotype"], r["stage"]): float(r["within_track_zscore"]) for r in rna_rows}
    for marker in MARKERS:
        for genotype in ("FL", "TN"):
            values = np.array([z[(marker, genotype, stage)] for stage in stages], dtype=float)
            peak = stages[int(np.argmax(values))]
            if summary[(marker, genotype)]["RNA_peak_stage"] != peak:
                raise AssertionError(f"LATEST114 peak mismatch for {marker} {genotype}: {peak}")

    protein = {}
    for row in prot_rows:
        protein[(row["family"], row["genotype"], row["stage"])] = int(
            row["detected_at_least_two_replicates"]
        )

    delta = np.array([[z[(m, "FL", s)] - z[(m, "TN", s)] for s in stages] for m in MARKERS])
    flz = np.array([[z[(m, "FL", s)] for s in stages] for m in MARKERS])
    tnz = np.array([[z[(m, "TN", s)] for s in stages] for m in MARKERS])

    fig = plt.figure(figsize=(11.0, 7.3))
    gs = fig.add_gridspec(2, 1, height_ratios=(1.06, 1.28), hspace=0.42,
                          left=0.125, right=0.715, top=0.865, bottom=0.09)
    ax = fig.add_subplot(gs[0])
    image = ax.imshow(delta, cmap=CM_DELTA, vmin=-2.4, vmax=2.4, aspect="auto")
    ax.axvline(n_dev - 0.5, color="#5C666B", lw=0.75)
    ax.axvspan(n_dev - 0.5, len(stages) - 0.5, color="#555555", alpha=0.07, zorder=2)
    ax.set_xticks(range(len(stages)))
    ax.set_xticklabels(stages, rotation=90, fontsize=5.9)
    ax.set_yticks(range(len(MARKERS)))
    ax.set_yticklabels(MARKERS, fontsize=6.8)
    ax.set_xticks(np.arange(-0.5, len(stages), 1), minor=True)
    ax.grid(which="minor", color="white", lw=0.45)
    ax.tick_params(which="minor", length=0)
    bounds = [("Push", 0, 3), ("Pull", 4, 6), ("Package / protect", 7, 9)]
    for name, start, end in bounds[:-1]:
        ax.axhline(end + 0.5, color="#4C565C", lw=0.8)
    for name, start, end in bounds:
        ax.text(-0.135, (start + end) / 2, name.replace(" / ", "/"),
                transform=ax.get_yaxis_transform(), rotation=90, ha="center",
                va="center", fontsize=6.4, color=MUTED)
    for row in (6, 8):
        ax.add_patch(Rectangle((-0.5, row - 0.5), len(stages), 1.0, fill=False,
                               edgecolor="#B77C2B", lw=0.9, zorder=5))
    ax.text(5.5, -1.10, "development", ha="center", va="center", fontsize=5.5, color=MUTED)
    ax.text(15.5, -1.10, "post-harvest", ha="center", va="center", fontsize=5.5, color=MUTED)
    for i, marker in enumerate(MARKERS):
        fl = summary[(marker, "FL")]
        tn = summary[(marker, "TN")]
        protein_fl = float(fl["detection_percent"])
        protein_tn = float(tn["detection_percent"])
        ax.text(18.2, i, f"{fl['RNA_peak_stage']} / {tn['RNA_peak_stage']}",
                ha="left", va="center", fontsize=5.7, color=INK)
        ax.text(21.2, i, f"{protein_fl:.0f}% / {protein_tn:.0f}%",
                ha="left", va="center", fontsize=5.7, color="#4C565C")
    ax.set_xlim(-0.5, 24.9)
    ax.set_title("Δ RNA = FL − TN", loc="left", fontsize=8.6, pad=5)
    ax.text(18.2, -1.10, "RNA peak FL/TN", fontsize=5.8, color=MUTED, ha="left")
    ax.text(21.2, -1.10, "Prot. %", fontsize=5.8, color=MUTED, ha="left")

    ax2 = fig.add_subplot(gs[1])
    for i, marker in enumerate(MARKERS):
        for values, color, yoffset in ((flz[i], FL_C, -0.18), (tnz[i], TN_C, 0.18)):
            on = values >= 1.0
            start = 0
            while start < len(stages):
                if not on[start]:
                    start += 1
                    continue
                end = start
                while end + 1 < len(stages) and on[end + 1]:
                    end += 1
                ax2.plot([start - 0.36, end + 0.36], [i + yoffset, i + yoffset],
                         color=color, lw=3.8, solid_capstyle="round", zorder=3)
                start = end + 1
            peak = int(np.argmax(values))
            ax2.plot(peak, i + yoffset, marker="o", ms=3.2, mec="#1B252B", mew=0.55,
                     color=color, zorder=4)
    ax2.axvline(n_dev - 0.5, color="#5C666B", lw=0.75)
    ax2.axvspan(n_dev - 0.5, len(stages) - 0.5, color="#555555", alpha=0.07, zorder=0)
    for _, start, end in bounds[:-1]:
        ax2.axhline(end + 0.5, color="#4C565C", lw=0.8)
    ax2.set_xlim(-0.5, len(stages) - 0.5)
    ax2.set_ylim(len(MARKERS) - 0.5, -0.5)
    ax2.set_xticks(range(len(stages)))
    ax2.set_xticklabels(stages, rotation=90, fontsize=5.9)
    ax2.set_yticks(range(len(MARKERS)))
    ax2.set_yticklabels(MARKERS, fontsize=6.8)
    ax2.set_title("RNA activity windows (z ≥ 1); upper lane FL, lower lane TN; ● = RNA peak",
                  loc="left", fontsize=8.1, pad=5)
    for side in ("top", "right"):
        ax2.spines[side].set_visible(False)
    ax2.legend(handles=[Line2D([], [], color=FL_C, lw=3.5, label="FL"),
                        Line2D([], [], color=TN_C, lw=3.5, label="TN")],
               frameon=False, fontsize=6.2, loc="lower right", ncol=2, handlelength=1.7)

    cbar = fig.colorbar(image, ax=ax, fraction=0.033, pad=0.05)
    cbar.set_label("within-track z difference", fontsize=6.2)
    cbar.ax.tick_params(labelsize=5.7)

    fad2_170 = _protein_ratio(prot_rows, "FAD2", "170d")
    ole16_185 = _protein_ratio(prot_rows, "OLE16", "185d")
    fad2_fl = float(summary[("FAD2-like", "FL")]["detection_percent"])
    fad2_tn = float(summary[("FAD2-like", "TN")]["detection_percent"])
    ole16_fl_peak, ole16_tn_peak = _protein_peak_stage(prot_rows, "OLE16")

    fig.text(0.035, 0.965,
             "f  Marker-family timing separates FL and TN without a single uniform switch point",
             fontsize=9.2, fontweight="bold", color=INK, va="top")
    fig.text(0.035, 0.932,
             "Joint DESeq2 normalization across 114 RNA samples; red = higher in FL,  "
             "blue = higher in TN;  protein detection: ≥2/3 Astral-114 replicates.",
             fontsize=6.7, color="#38434A", va="top")

    story = (
        "STORY AT A GLANCE\n"
        "\n"
        "TN activates the upstream Push\n"
        "programme late (onset ~110–125 d):\n"
        "ACCase, KASIII, ENR, FATA/B stay\n"
        "aligned through late development\n"
        "and post-harvest.\n"
        "\n"
        f"FAD2 (chr08): RNA 155d/170d; protein\n"
        f"{fad2_fl:.0f}% vs {fad2_tn:.0f}% stages;\n"
        f"170d FL/TN = {fad2_170:.2f}×.\n"
        "\n"
        f"OLE16: RNA peak 170d both; protein\n"
        f"185d FL/TN = {ole16_185:.1f}×;\n"
        f"protein peak {ole16_fl_peak} FL vs {ole16_tn_peak} TN.\n"
        "\n"
        "LOX stays FL-timed (peak 80 d) as the\n"
        "panel’s loss-risk marker.\n"
    )
    fig.text(0.75, 0.73, "STORY AT A GLANCE", fontsize=8.1, fontweight="bold", color=INK)
    fig.text(0.75, 0.695, story.split("\n", 1)[1], fontsize=6.7, color="#38434A",
             va="top", linespacing=1.45)
    fig.text(0.75, 0.385, "Reading the panel", fontsize=8.1, fontweight="bold", color=INK)
    fig.text(0.75, 0.350,
             "Top: stage-wise FL−TN RNA contrast.\n"
             "Bottom: each family’s positive RNA\n"
             "window with the maximum marked by\n"
             "a dot. Peak labels are from the\n"
             "corrected LATEST114 summary and\n"
             "match the RNA table 20/20.",
             fontsize=6.7, color="#38434A", va="top", linespacing=1.45)
    save_pair(fig, "Figure2f_program_delta_LATEST114")


def _protein_ratio(rows: list[dict[str, str]], family: str, stage: str) -> float:
    values = {}
    for r in rows:
        if r["family"] == family and r["stage"] == stage:
            if r["mean_detected_replicate_abundance"] not in ("", "NA"):
                values[r["genotype"]] = float(r["mean_detected_replicate_abundance"])
    if "FL" not in values or "TN" not in values or values["TN"] <= 0:
        raise AssertionError(f"Cannot compute protein ratio {family} {stage}")
    return values["FL"] / values["TN"]


def _protein_peak_stage(rows: list[dict[str, str]], family: str) -> tuple[str, str]:
    out = {}
    for genotype in ("FL", "TN"):
        best, best_stage = -1.0, None
        for r in rows:
            if r["family"] == family and r["genotype"] == genotype:
                if r["mean_detected_replicate_abundance"] in ("", "NA"):
                    continue
                value = float(r["mean_detected_replicate_abundance"])
                if value > best:
                    best, best_stage = value, r["stage"]
        if best_stage is None:
            raise AssertionError(f"No protein abundance for {family} {genotype}")
        out[genotype] = best_stage
    return out["FL"], out["TN"]


def write_manifest() -> None:
    inputs = [D_TABLE, D_ANCT, E_AXIS, E_MET, E_COV, F_RNA, F_MARK, F_PROT]
    with (OUT / "Figure2_def_latest114_manifest.tsv").open("w") as handle:
        handle.write("path\tsha256\tsize_bytes\n")
        for path in inputs:
            digest = hashlib.sha256(path.read_bytes()).hexdigest()
            handle.write(f"{path}\t{digest}\t{path.stat().st_size}\n")


if __name__ == "__main__":
    build_d()
    build_e()
    build_f()
    write_manifest()
    print(f"Wrote redraws to {OUT}")