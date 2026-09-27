#!/usr/bin/env python3
"""Restore the original 0918 d/e/f panel styles (as in Figure2-A4 assembly),
but keep the corrected LATEST114 data scope from the 2026-09-19 handoff.

Styles (identical to build_0918_panels.py so the printed assembly does not change):
  d  = dot-matrix of curated locus counts (24 enzyme classes x 9 genomes),
       grey row band = identical across the four focal palms, pink cell = deviates
       from coconut, right blocks = Delta-max + ancestral gains/losses/net.
  e  = four constructed-axis trajectories (top) + FL / TN / delta heatmaps (bottom).
  f  = FL-TN delta RNA heatmap (top) + RNA activity windows as gantt lanes,
       RNA peaks taken from the LATEST114 marker summary and validated 20/20
       against the within-track z argmax.
"""

from __future__ import annotations

import csv
from pathlib import Path

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle

BASE = Path("${ANALYSIS_DIR}/22_answer_reviews/00_ms")
OUT = BASE / "05_MS" / "0918_revision"
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

FL_C = "#2E8B8B"
TN_C = "#D9604A"
MODULE_C = {
    "Plastid fatty-acid synthesis": "#C98A2C",
    "FA export and modification": "#E0B25C",
    "VLCFA side branch": "#3E8E7E",
    "ER TAG assembly": "#3E6FA3",
    "PC–DAG exchange and desaturation": "#7FA7C9",
}
CM_DELTA = LinearSegmentedColormap.from_list("delta", ["#2C6FA6", "#F3F5F7", "#C8442E"])
CM_BLUE = LinearSegmentedColormap.from_list("blues2", ["#FFFFFF", "#2C6FA6"])

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 7.2,
    "axes.linewidth": 0.6,
    "xtick.major.width": 0.6,
    "ytick.major.width": 0.6,
    "text.color": "black",
    "axes.labelcolor": "black",
    "xtick.color": "black",
    "ytick.color": "black",
    "savefig.dpi": 600,
})


def load_stages():
    rows = list(csv.DictReader(open(E_AXIS), delimiter="\t"))
    idx = {}
    for r in rows:
        idx[int(r["stage_index"])] = (r["stage"], r["phase"])
    return [idx[i] for i in sorted(idx)]


def save_pair(fig, stem):
    fig.savefig(OUT / f"{stem}.pdf", facecolor="white")
    fig.savefig(OUT / f"{stem}.png", dpi=600, facecolor="white")
    plt.close(fig)


# ----------------------------------------------------------------------------
# Panel d - dot matrix (original style)
# ----------------------------------------------------------------------------
def build_d():
    rows = list(csv.DictReader(open(D_TABLE), delimiter="\t"))
    order_mod = ["Plastid fatty-acid synthesis", "FA export and modification",
                 "VLCFA side branch", "ER TAG assembly", "PC–DAG exchange and desaturation"]
    genomes = ["Calamus", "Daemonorops", "Nypa_fruticans", "Phoenix_dactylifera",
               "Areca_catechu", "Cocos_nucifera", "American_hap1", "Dura", "Pisifera"]
    label = {"Calamus": "Calamus", "Daemonorops": "D. draco", "Nypa_fruticans": "Nypa",
             "Phoenix_dactylifera": "Phoenix", "Areca_catechu": "Areca",
             "Cocos_nucifera": "Coconut", "American_hap1": "Oleifera",
             "Dura": "Dura", "Pisifera": "Pisifera"}
    focal = ["Cocos_nucifera", "American_hap1", "Dura", "Pisifera"]

    cnt, module = {}, {}
    for r in rows:
        cnt[(r["Enzyme"], r["Genome"])] = int(r["Curated_gene_locus_count"])
        module[r["Enzyme"]] = r["Pathway_section"]
    enzymes = []
    for m in order_mod:
        enzymes += [e for e in dict.fromkeys(r["Enzyme"] for r in rows) if module[e] == m]

    ident = sum(1 for e in enzymes if len(set(cnt[(e, g)] for g in focal)) == 1)

    anc = {}
    for r in csv.DictReader(open(D_ANCT), delimiter="\t"):
        anc[r["Enzyme"]] = (int(r["Total_inferred_gains"]), int(r["Total_inferred_losses"]),
                            int(r["Net_change"]))
    anc_map = {
        "ACCase–BCCP1/2": "ACCase", "ACCase–BC / CAC2": "ACCase", "ACCase–CTα / CAC3": "ACCase",
        "MCAT / FabD": "FabD (MCAT)", "KASIII": "KASIII", "KASI": "KAS I/II",
        "KASII / FAB1": "KAS I/II", "KAR / FabG": "FabG (KAR)", "ENR / FabI (MOD1)": "FabI (ENR)",
        "SAD family": "SAD", "FATA": "FATA/B", "FATB": "FATA/B", "LACS9-like": "LACS",
        "KCS / FAE (VLCFA side branch)": "KCS", "GPAT9": "GPAT", "LPAT2": "LPAT",
        "PAH1/2 (PAP)": "PAP", "DGAT1": "DGAT", "DGAT2": "DGAT", "PDCT / ROD1": "PDCT",
        "FAD2/6 (omega-6 FAD family)": "FAD2", "FAD3/7/8 (omega-3 FAD family)": "FAD3/7/8",
    }
    anc_shown = set()
    anc_row = {}
    for e in enzymes:
        src = anc_map.get(e)
        if src is None or src not in anc or src in anc_shown:
            anc_row[e] = None
            continue
        anc_shown.add(src)
        anc_row[e] = anc[src]

    fig_w, fig_h = 7.4, 6.3
    fig = plt.figure(figsize=(fig_w, fig_h))
    ax = fig.add_axes([0.30, 0.095, 0.60, 0.735])
    n_r, n_c = len(enzymes), len(genomes)

    for i, e in enumerate(enzymes):
        ax.add_patch(Rectangle((-1.35, i - 0.5), 0.22, 1.0,
                               color=MODULE_C[module[e]], clip_on=False, lw=0))
        row_same = len(set(cnt[(e, g)] for g in focal)) == 1
        if row_same:
            ax.add_patch(Rectangle((-1.1, i - 0.5), n_c + 1.1, 1.0,
                                   color="#F2F3F5", zorder=0, lw=0))
        for j, g in enumerate(genomes):
            v = cnt[(e, g)]
            dev = j >= 5 and v != cnt[(e, "Cocos_nucifera")]
            ax.add_patch(Rectangle((j - 0.5, i - 0.5), 1, 1,
                                   facecolor="#F6D9D2" if dev else "none",
                                   edgecolor="none", zorder=1))
            shade = 0.12 + 0.5 * min(v, 30) / 30
            ax.add_patch(plt.Circle((j, i), 0.30, color=(0.17, 0.44, 0.65), alpha=shade,
                                    zorder=2, lw=0))
            ax.text(j, i, str(v), ha="center", va="center", fontsize=6.6, zorder=3,
                    color="white" if v >= 12 else "black")
        dmax = max(cnt[(e, g)] for g in focal) - min(cnt[(e, g)] for g in focal)
        ax.text(9.7, i, "0" if dmax == 0 else f"+{dmax}",
                ha="center", va="center", fontsize=6.8,
                color="#888888" if dmax == 0 else "#B03A2E")
        av = anc_row[e]
        if av is not None:
            g_, l_, n_ = av
            ax.text(10.6, i, str(g_) if g_ else "·", ha="center", va="center", fontsize=6.8,
                    color="#B03A2E" if g_ else "#C4C4C4")
            ax.text(11.2, i, str(l_) if l_ else "·", ha="center", va="center", fontsize=6.8,
                    color="#2C6FA6" if l_ else "#C4C4C4")
            ax.text(11.8, i, f"{n_:+d}" if n_ else "0", ha="center", va="center", fontsize=6.8,
                    color=("#B03A2E" if n_ > 0 else "#2C6FA6") if n_ else "#888888")

    ax.set_xlim(-1.4, n_c + 3.4)
    ax.set_ylim(n_r - 0.5, -1.95)
    ax.set_xticks(range(n_c))
    ax.set_xticklabels([label[g] for g in genomes], rotation=38, ha="right", fontsize=7.4)
    for t, g in zip(ax.get_xticklabels(), genomes):
        t.set_color("black" if g in focal else "#8A8A8A")
        if g in focal:
            t.set_fontweight("bold")
    ax.set_yticks(range(n_r))
    ax.set_yticklabels([e.replace(" (MOD1)", "").replace(" / ", "/") for e in enzymes], fontsize=7.0)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for j in (5, 6, 7, 8):
        ax.plot([j - 0.5, j - 0.5], [-1.95, n_r - 0.5], color="#C9CDD3", lw=0.7, zorder=1)
    ax.plot([4.5, 4.5], [-1.95, n_r - 0.5], color="black", lw=0.9, zorder=1)
    ax.text(9.7, -1.10, "Δmax", fontsize=7.4, ha="center", va="center")
    ax.text(10.6, -1.10, "+", fontsize=7.4, ha="center", va="center", color="#B03A2E")
    ax.text(11.2, -1.10, "−", fontsize=7.4, ha="center", va="center", color="#2C6FA6")
    ax.text(11.8, -1.10, "net", fontsize=7.4, ha="center", va="center")
    ax.text(11.2, -1.62, "ancestral, oil-palm ancestor", fontsize=6.6, ha="center", va="center",
            color="#555555")
    ax.text(6.5, -1.10, "coconut + 3 Elaeis genomes", ha="center", fontsize=7.6,
            fontweight="bold")
    ax.text(2.0, -1.10, "other palms (context)", ha="center", fontsize=7.4, color="#8A8A8A")

    fig.text(0.02, 0.975, "d  Curated locus counts of 24 core lipid-pathway enzyme classes",
             fontsize=9.2, fontweight="bold", va="top")
    fig.text(0.02, 0.948,
             f"Grey rows: {ident}/24 classes with identical locus counts in coconut and the three Elaeis genomes; pink cells differ from coconut.",
             fontsize=7.0, va="top")
    fig.text(0.02, 0.925,
             "Right blocks: Δmax across the four focal genomes; ancestral reconstruction at the oil-palm ancestor (gains +, losses −, net).",
             fontsize=7.0, va="top", color="#333333")
    fig.text(0.02, 0.899,
             "FabG and FabI each gained two ancestral copies; KCS/FAE is the VLCFA side branch, not part of the core TAG pathway.",
             fontsize=7.0, va="top", color="#333333")
    legend_short = {
        "Plastid fatty-acid synthesis": "Plastid FA synthesis",
        "FA export and modification": "FA export",
        "VLCFA side branch": "VLCFA branch",
        "ER TAG assembly": "ER TAG assembly",
        "PC–DAG exchange and desaturation": "PC–DAG / desaturation",
    }
    lx = 0.05
    for m in order_mod:
        fig.add_artist(Rectangle((lx, 0.862), 0.012, 0.013, color=MODULE_C[m],
                                 transform=fig.transFigure))
        fig.text(lx + 0.017, 0.868, legend_short[m], fontsize=6.4, va="center", color="#333333")
        lx += 0.017 + 0.0062 * len(legend_short[m]) + 0.026
    save_pair(fig, "Figure2d_dosage_dotmatrix_LATEST")
    print("d done")


# ----------------------------------------------------------------------------
# Panel e - original style (top trajectories + FL/TN/delta heatmaps)
# ----------------------------------------------------------------------------
def build_e():
    stages = load_stages()
    stage_names = [s for s, p in stages]
    n_dev = sum(1 for s, p in stages if p == "development")
    ax_rows = list(csv.DictReader(open(E_AXIS), delimiter="\t"))
    cov = {r["axis_id"]: r for r in csv.DictReader(open(E_COV), delimiter="\t")}
    met = list(csv.DictReader(open(E_MET), delimiter="\t"))

    axes_order = ["P02", "P01", "P03", "P04"]
    axes_title = {"P02": "Storage lipids", "P01": "Oleic balance",
                  "P03": "Hydrolytic deterioration", "P04": "Oxidative deterioration"}

    series = {}
    for r in ax_rows:
        series.setdefault((r["axis_id"], r["genotype"]), {})[r["stage"]] = float(r["mean"])

    fig = plt.figure(figsize=(10.6, 7.6))
    gs = fig.add_gridspec(2, 4, height_ratios=[1.02, 1.55], hspace=0.52, wspace=0.30,
                          left=0.105, right=0.975, top=0.885, bottom=0.075)

    for k, aid in enumerate(axes_order):
        a = fig.add_subplot(gs[0, k])
        fl = np.array([series[(aid, "FL")][s] for s in stage_names])
        tn = np.array([series[(aid, "TN")][s] for s in stage_names])
        x = np.arange(len(stage_names))
        a.axvspan(n_dev - 0.5, len(stage_names) - 0.5, color="#F0F0F0", zorder=0)
        a.fill_between(x, fl, tn, where=fl >= tn, color=FL_C, alpha=0.28, lw=0, zorder=1)
        a.fill_between(x, fl, tn, where=fl < tn, color=TN_C, alpha=0.28, lw=0, zorder=1)
        a.plot(x, fl, color=FL_C, lw=1.4, marker="o", ms=2.6, label="FL", zorder=3)
        a.plot(x, tn, color=TN_C, lw=1.4, marker="o", ms=2.6, label="TN", zorder=3)
        a.axhline(0, color="#BBBBBB", lw=0.6, zorder=0)
        a.set_title(axes_title[aid], fontsize=8.6, pad=13)
        c = cov[aid]
        a.text(0.0, 1.015,
               f"{c['level2_compound_n']} Level-2 + {c['ms1_mass_proxy_compound_n']} MS1 proxies",
               transform=a.transAxes, fontsize=6.3, color="#555555")
        a.set_xticks(x[::2])
        a.set_xticklabels([stage_names[i] for i in range(0, len(stage_names), 2)],
                          rotation=45, ha="right", fontsize=6.0)
        a.set_xlim(-0.6, len(stage_names) - 0.4)
        a.tick_params(labelsize=6.4)
        if k == 0:
            a.set_ylabel("Constructed axis score", fontsize=7.4)
        else:
            a.set_yticklabels([])
        if k == 3:
            a.legend(handles=[Line2D([], [], color=FL_C, marker="o", ms=3, lw=1.4, label="FL"),
                              Line2D([], [], color=TN_C, marker="o", ms=3, lw=1.4, label="TN")],
                     loc="lower right", fontsize=6.4, frameon=False, ncol=2,
                     handlelength=1.4, columnspacing=1.0)

    met_order = []
    for aid in axes_order:
        ms = [m for m in met if m["axis_id"] == aid]
        names = list(dict.fromkeys(m["display_name"] for m in ms))
        for nm in names:
            met_order.append((aid, nm))
    n_m = len(met_order)
    gsub = gs[1, :].subgridspec(1, 3, wspace=0.10)
    ax_fl = fig.add_subplot(gsub[0, 0])
    ax_tn = fig.add_subplot(gsub[0, 1])
    ax_dl = fig.add_subplot(gsub[0, 2])

    mvals = {}
    for r in met:
        mvals[(r["display_name"], r["genotype"], r["stage"])] = (
            float(r["mean_consensus_z"]) if r["mean_consensus_z"] not in ("NA", "") else np.nan)

    def mat(geno):
        M = np.full((n_m, len(stage_names)), np.nan)
        for i, (aid, nm) in enumerate(met_order):
            for j, s in enumerate(stage_names):
                M[i, j] = mvals.get((nm, geno, s), np.nan)
        return M

    MF, MT = mat("FL"), mat("TN")
    MD = MF - MT

    for a_, M, ttl in ((ax_fl, MF, "FL"), (ax_tn, MT, "TN")):
        a_.imshow(M, cmap=CM_BLUE, vmin=0, vmax=2.6, aspect="auto")
        a_.set_title(ttl, fontsize=8.4, color=FL_C if ttl == "FL" else TN_C, pad=4)
        a_.set_xticks(range(len(stage_names)))
        a_.set_xticklabels(stage_names, rotation=90, fontsize=5.8)
        a_.set_yticks(range(n_m))
        a_.set_yticklabels([nm for aid, nm in met_order], fontsize=6.6)
        a_.axvline(n_dev - 0.5, color="#666666", lw=0.8)
        a_.set_xticks(np.arange(-0.5, len(stage_names), 1), minor=True)
        a_.grid(which="minor", color="white", lw=0.5)
        a_.tick_params(which="minor", length=0)
    ax_tn.set_yticklabels([])
    im = ax_dl.imshow(MD, cmap=CM_DELTA, vmin=-1.6, vmax=1.6, aspect="auto")
    ax_dl.set_title("Δ = FL − TN", fontsize=8.4, pad=4)
    ax_dl.set_xticks(range(len(stage_names)))
    ax_dl.set_xticklabels(stage_names, rotation=90, fontsize=5.8)
    ax_dl.set_yticks([])
    ax_dl.axvline(n_dev - 0.5, color="#666666", lw=0.8)
    ax_dl.set_xticks(np.arange(-0.5, len(stage_names), 1), minor=True)
    ax_dl.grid(which="minor", color="white", lw=0.5)
    ax_dl.tick_params(which="minor", length=0)

    cb = fig.colorbar(im, ax=[ax_fl, ax_tn, ax_dl], fraction=0.02, pad=0.012)
    cb.set_label("Δ consensus z (FL − TN)", fontsize=6.6)
    cb.ax.tick_params(labelsize=6.0)

    ax_fl.text(-0.22, 1.17, "Representative compounds (consensus z-score)",
               transform=ax_fl.transAxes, fontsize=8.0, fontweight="bold")
    ax_fl.text(-0.22, 1.115,
               "stages: 13 development (0–185 d) + 6 post-harvest (12–72 h); TG 42:0 has n=1 (FL 95 d) and n=0 (TN 36 h, NA)",
               transform=ax_fl.transAxes, fontsize=5.7, color="#666666")

    fig.text(0.02, 0.975, "e  Constructed molecular-axis scores and representative compounds differentiate FL and TN",
             fontsize=9.6, fontweight="bold", va="top")
    fig.text(0.02, 0.945,
             "Axes are constructed composite scores dominated by putatively annotated MS1 features; values were not interpolated. "
             "Shaded area in the top row is the FL−TN gap.",
             fontsize=7.0, va="top")
    save_pair(fig, "Figure2e_axis_delta_LATEST")
    print("e done")


# ----------------------------------------------------------------------------
# Panel f - original style with LATEST114 peaks (validated 20/20)
# ----------------------------------------------------------------------------
def build_f():
    stages = load_stages()
    stage_names = [s for s, p in stages]
    n_dev = sum(1 for s, p in stages if p == "development")

    rna = list(csv.DictReader(open(F_RNA), delimiter="\t"))
    summ = {(r["marker"], r["genotype"]): r for r in csv.DictReader(open(F_MARK), delimiter="\t")}
    prot = list(csv.DictReader(open(F_PROT), delimiter="\t"))

    markers, groups = [], {}
    for r in rna:
        markers.append(r["marker"])
        groups[r["marker"]] = r["process"]
    preferred = ["ACCase", "KASIII", "ENR", "FATA/B", "LACS", "GPAT", "FAD2-like",
                 "DGAT", "OLE16", "LOX (loss risk)"]
    markers = [m for m in preferred if m in set(markers)]
    groups_order = ["Push", "Pull", "Package / protect"]

    Z = {}
    for r in rna:
        Z[(r["marker"], r["genotype"], r["stage"])] = float(r["within_track_zscore"])
    # LATEST114 peak validation: argmax of the same RNA table == summary table.
    for m in markers:
        for ge in ("FL", "TN"):
            values = np.array([Z[(m, ge, s)] for s in stage_names], dtype=float)
            peak = stage_names[int(np.argmax(values))]
            if summ[(m, ge)]["RNA_peak_stage"] != peak:
                raise AssertionError(f"LATEST114 peak mismatch {m} {ge}: {peak}")
    P = {}
    alias_back = {"FabI (ENR)": "ENR", "FAD2": "FAD2-like"}
    for r in prot:
        mk = alias_back.get(r["family"], r["family"])
        P[(mk, r["genotype"], r["stage"])] = int(r["detected_at_least_two_replicates"])

    D = np.array([[Z[(m, "FL", s)] - Z[(m, "TN", s)] for s in stage_names] for m in markers])
    FLz = np.array([[Z[(m, "FL", s)] for s in stage_names] for m in markers])
    TNz = np.array([[Z[(m, "TN", s)] for s in stage_names] for m in markers])
    PF = np.array([[P[(m, "FL", s)] for s in stage_names] for m in markers])
    PT = np.array([[P[(m, "TN", s)] for s in stage_names] for m in markers])

    fig = plt.figure(figsize=(10.6, 7.4))
    gs = fig.add_gridspec(2, 1, height_ratios=[1.35, 1.30], hspace=0.30,
                          left=0.150, right=0.985, top=0.855, bottom=0.075)

    ax = fig.add_subplot(gs[0])
    im = ax.imshow(D, cmap=CM_DELTA, vmin=-2.4, vmax=2.4, aspect="auto")
    ax.set_xticks(range(len(stage_names)))
    ax.set_xticklabels(stage_names, rotation=90, fontsize=6.2)
    ax.set_yticks(range(len(markers)))
    ax.set_yticklabels(markers, fontsize=7.2)
    ax.axvline(n_dev - 0.5, color="#444444", lw=1.0)
    ax.set_xticks(np.arange(-0.5, len(stage_names), 1), minor=True)
    ax.grid(which="minor", color="white", lw=0.5)
    ax.tick_params(which="minor", length=0)
    bounds = []
    for gname in groups_order:
        idxs = [i for i, m in enumerate(markers) if groups[m] == gname]
        bounds.append((gname, min(idxs), max(idxs)))
    for gname, i0, i1 in bounds[:-1]:
        ax.axhline(i1 + 0.5, color="#444444", lw=0.8)
    for gname, i0, i1 in bounds:
        ax.text(-0.155, (i0 + i1) / 2, gname.replace(" / ", "/"),
                transform=ax.get_yaxis_transform(), rotation=90, va="center",
                ha="center", fontsize=6.8, color="#555555")
    ax.axvspan(n_dev - 0.5, len(stage_names) - 0.5, color="#8C8C8C", alpha=0.10, zorder=2)
    for i, m in enumerate(markers):
        fl = summ[(m, "FL")]
        tn = summ[(m, "TN")]
        ax.text(len(stage_names) - 0.3, i,
                f"RNA peak: {fl['RNA_peak_stage']} / {tn['RNA_peak_stage']}",
                fontsize=6.0, va="center", ha="left", color="#333333")
        ax.text(len(stage_names) + 5.2, i,
                f"protein {float(fl['detection_percent']):.0f}% / {float(tn['detection_percent']):.0f}%",
                fontsize=6.0, va="center", ha="left", color="#333333")
    ax.set_xlim(-0.5, len(stage_names) + 10.6)
    ax.set_title("Δ RNA timing (FL − TN), family-level z-scores", fontsize=8.8, pad=6)

    ax2 = fig.add_subplot(gs[1])
    for i, m in enumerate(markers):
        for k, (M, pc, col) in enumerate(((FLz, PF, FL_C), (TNz, PT, TN_C))):
            y = i - 0.19 if k == 0 else i + 0.19
            on = M[i] >= 1.0
            j = 0
            while j < len(stage_names):
                if on[j]:
                    j2 = j
                    while j2 + 1 < len(stage_names) and on[j2 + 1]:
                        j2 += 1
                    ax2.plot([j - 0.40, j2 + 0.40], [y, y], color=col, lw=3.6,
                             solid_capstyle="round", zorder=3)
                    j = j2 + 1
                else:
                    j += 1
            pj = int(np.argmax(M[i]))
            ax2.plot(pj, y, marker="o", ms=3.2, color=col, mec="black", mew=0.5, zorder=4)
    ax2.axvspan(n_dev - 0.5, len(stage_names) - 0.5, color="#8C8C8C", alpha=0.10, zorder=0)
    ax2.set_xlim(-0.5, len(stage_names) + 10.6)
    ax2.set_ylim(len(markers) - 0.5, -0.5)
    ax2.set_xticks(range(len(stage_names)))
    ax2.set_xticklabels(stage_names, rotation=90, fontsize=6.2)
    ax2.set_yticks(range(len(markers)))
    ax2.set_yticklabels(markers, fontsize=7.2)
    ax2.axvline(n_dev - 0.5, color="#444444", lw=1.0)
    for gname, i0, i1 in bounds[:-1]:
        ax2.axhline(i1 + 0.5, color="#444444", lw=0.8)
    for side in ("top", "right"):
        ax2.spines[side].set_visible(False)
    ax2.set_title("RNA activity windows (z ≥ 1; upper lane FL, lower lane TN; ● RNA peak). "
                  "Protein detection per family is summarised at right.",
                  fontsize=8.4, pad=6)

    fad2_170 = _protein_ratio(prot, "FAD2", "170d")
    ole16_185 = _protein_ratio(prot, "OLE16", "185d")
    ole16_fl, ole16_tn = _protein_peak_stage(prot, "OLE16")
    yi = markers.index("OLE16")
    ax2.text(19.3, yi - 0.17, f"OLE16: 185 d protein FL/TN = {ole16_185:.1f}×",
             fontsize=6.4, va="center", color="#333333")
    ax2.text(19.3, yi + 0.17, f"protein peak: {ole16_fl} (FL) / {ole16_tn} (TN)",
             fontsize=6.4, va="center", color="#333333")
    yi = markers.index("FAD2-like")
    fad2_fl_n = summ[("FAD2-like", "FL")]["stage_detected_at_least_two_replicates_n"]
    fad2_tn_n = summ[("FAD2-like", "TN")]["stage_detected_at_least_two_replicates_n"]
    ax2.text(19.3, yi - 0.17,
             f"FAD2 (chr08): protein detected {fad2_fl_n}/19 (FL) vs {fad2_tn_n}/19 (TN)",
             fontsize=6.4, va="center", color="#333333")
    ax2.text(19.3, yi + 0.17, f"170 d protein FL/TN = {fad2_170:.2f}×",
             fontsize=6.4, va="center", color="#333333")

    cb = fig.colorbar(im, ax=ax, fraction=0.028, pad=0.30)
    cb.set_label("Δ RNA z (FL − TN)", fontsize=6.8)
    cb.ax.tick_params(labelsize=6.2)

    fig.text(0.02, 0.975, "f  Stage-resolved RNA programmes and protein detection separate FL and TN",
             fontsize=9.6, fontweight="bold", va="top")
    fig.text(0.02, 0.945,
             "Red = higher in FL, blue = higher in TN. Within-track z-scores (joint DESeq2, 114 samples); "
             "protein detection = ≥2 of 3 biological replicates (Astral-114 directLFQ; FAD2 restricted to chr08 orthogroup). "
             "RNA peaks are from the corrected LATEST114 summary (20/20 check passed).",
             fontsize=7.0, va="top")
    save_pair(fig, "Figure2f_program_delta_LATEST114")
    print("f done")


def _protein_ratio(rows, family, stage):
    values = {}
    for r in rows:
        if r["family"] == family and r["stage"] == stage:
            if r["mean_detected_replicate_abundance"] not in ("", "NA"):
                values[r["genotype"]] = float(r["mean_detected_replicate_abundance"])
    if "FL" not in values or "TN" not in values or values["TN"] <= 0:
        raise AssertionError(f"Cannot compute protein ratio {family} {stage}")
    return values["FL"] / values["TN"]


def _protein_peak_stage(rows, family):
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


if __name__ == "__main__":
    build_d()
    build_e()
    build_f()
    print("all done ->", OUT)