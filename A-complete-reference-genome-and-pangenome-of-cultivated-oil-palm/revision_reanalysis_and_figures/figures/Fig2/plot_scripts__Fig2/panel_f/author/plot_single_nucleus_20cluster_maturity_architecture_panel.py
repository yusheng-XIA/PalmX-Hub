#!/usr/bin/env python3
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


BASE = Path("${ANALYSIS_DIR}/22_answer_reviews/00_ms/03_V3/03_figure3/03_single")
CELLS = BASE / "single_nucleus_20cluster_main_panel_cells.tsv"
FOCUS_COMP = BASE / "single_nucleus_20cluster_focus_cluster_composition_185d.tsv"
FAD2 = BASE / "single_nucleus_20cluster_fad2_detection_by_sample.tsv"
OUT_PREFIX = BASE / "single_nucleus_20cluster_maturity_architecture_panel"


mpl.rcParams.update({
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "font.family": "DejaVu Sans",
    "font.size": 8,
    "axes.linewidth": 0.6,
    "axes.labelsize": 8,
    "axes.titlesize": 9,
    "xtick.labelsize": 7,
    "ytick.labelsize": 7,
    "legend.fontsize": 7,
})


def style_umap(ax, xlim=None, ylim=None):
    ax.set_aspect("equal", adjustable="box")
    if xlim is not None:
        ax.set_xlim(xlim)
    if ylim is not None:
        ax.set_ylim(ylim)
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)


def add_panel_label(ax, label):
    ax.text(
        -0.04,
        1.04,
        label,
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        fontsize=11,
        fontweight="bold",
    )


def plot_umap(ax, df, title, subtitle=None, xlim=None, ylim=None):
    background = df[~df["cluster"].isin(FOCUS_CLUSTERS)]
    ax.scatter(
        background["umap_1"],
        background["umap_2"],
        s=0.35,
        c="#d7d7d7",
        alpha=0.35,
        linewidths=0,
        rasterized=True,
    )
    for cluster in FOCUS_CLUSTERS:
        sub = df[df["cluster"] == cluster]
        if sub.empty:
            continue
        ax.scatter(
            sub["umap_1"],
            sub["umap_2"],
            s=0.9,
            c=FOCUS_COLORS[cluster],
            alpha=0.78,
            linewidths=0,
            rasterized=True,
            label=FOCUS_LABELS[cluster],
        )
        x, y = sub[["umap_1", "umap_2"]].median()
        ax.text(
            x,
            y,
            cluster,
            ha="center",
            va="center",
            color="white",
            fontsize=7.5,
            fontweight="bold",
            bbox=dict(boxstyle="round,pad=0.18", facecolor=FOCUS_COLORS[cluster], edgecolor="none", alpha=0.92),
        )
    ax.set_title(title, pad=4, fontweight="bold")
    if subtitle:
        ax.text(0.02, 0.02, subtitle, transform=ax.transAxes, ha="left", va="bottom", color="#4b5563", fontsize=7)
    style_umap(ax, xlim=xlim, ylim=ylim)


FOCUS_CLUSTERS = ["C9", "C6", "C10"]
FOCUS_COLORS = {
    "C9": "#1f6f6a",
    "C6": "#c4472d",
    "C10": "#8f1d2c",
}
FOCUS_LABELS = {
    "C9": "C9 oil-storage",
    "C6": "C6 senescence",
    "C10": "C10 senescence",
}
VAR_COLORS = {"FL": "#1f6f6a", "TN": "#b81e2d"}
VAR_LABELS = {"FL": "SL", "TN": "TN"}


def main():
    cells = pd.read_csv(CELLS, sep="\t")
    comp = pd.read_csv(FOCUS_COMP, sep="\t")
    fad2 = pd.read_csv(FAD2, sep="\t")

    cells["timepoint"] = cells["timepoint"].astype(str).str.replace("d", "", regex=False)
    fad2["timepoint"] = fad2["timepoint"].astype(str).str.replace("d", "", regex=False)
    cells["variety_label"] = cells["variety"].map(VAR_LABELS).fillna(cells["variety"])
    fad2["variety_label"] = fad2["variety"].map(VAR_LABELS).fillna(fad2["variety"])
    pad_x = (cells["umap_1"].max() - cells["umap_1"].min()) * 0.035
    pad_y = (cells["umap_2"].max() - cells["umap_2"].min()) * 0.035
    umap_xlim = (cells["umap_1"].min() - pad_x, cells["umap_1"].max() + pad_x)
    umap_ylim = (cells["umap_2"].min() - pad_y, cells["umap_2"].max() + pad_y)

    fig = plt.figure(figsize=(10.2, 6.2), constrained_layout=False)
    gs = fig.add_gridspec(
        2,
        3,
        height_ratios=[1.18, 0.82],
        width_ratios=[1.0, 1.0, 1.0],
        hspace=0.24,
        wspace=0.12,
    )

    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[0, 1])
    ax_c = fig.add_subplot(gs[0, 2])
    ax_d = fig.add_subplot(gs[1, 0:2])
    ax_e = fig.add_subplot(gs[1, 2])

    plot_umap(
        ax_a,
        cells,
        "Single-nucleus atlas",
        "78,193 nuclei; 20 clusters",
        xlim=umap_xlim,
        ylim=umap_ylim,
    )
    plot_umap(
        ax_b,
        cells[(cells["variety"] == "FL") & (cells["timepoint"] == "185")],
        "SL at 185 d",
        f"n = {((cells['variety'] == 'FL') & (cells['timepoint'] == '185')).sum():,} nuclei",
        xlim=umap_xlim,
        ylim=umap_ylim,
    )
    plot_umap(
        ax_c,
        cells[(cells["variety"] == "TN") & (cells["timepoint"] == "185")],
        "TN at 185 d",
        f"n = {((cells['variety'] == 'TN') & (cells['timepoint'] == '185')).sum():,} nuclei",
        xlim=umap_xlim,
        ylim=umap_ylim,
    )
    add_panel_label(ax_a, "a")
    add_panel_label(ax_b, "b")
    add_panel_label(ax_c, "c")

    # Cluster origin at 185 d, matching the manuscript's C6/C9/C10 statements.
    comp = comp.copy()
    comp["percent"] = comp["cluster_composition_fraction"] * 100
    comp["cluster"] = pd.Categorical(comp["cluster"], categories=FOCUS_CLUSTERS, ordered=True)
    comp = comp.sort_values(["cluster", "variety"])
    y_pos = np.arange(len(FOCUS_CLUSTERS))
    left = np.zeros(len(FOCUS_CLUSTERS))
    for variety in ["FL", "TN"]:
        vals = []
        for cluster in FOCUS_CLUSTERS:
            row = comp[(comp["cluster"] == cluster) & (comp["variety"] == variety)]
            vals.append(float(row["percent"].iloc[0]) if not row.empty else 0)
        ax_d.barh(y_pos, vals, left=left, color=VAR_COLORS[variety], height=0.52, label=VAR_LABELS[variety])
        for i, val in enumerate(vals):
            if val >= 6:
                ax_d.text(left[i] + val / 2, i, f"{val:.1f}%", ha="center", va="center", color="white", fontsize=7, fontweight="bold")
            elif val > 0:
                ax_d.text(left[i] + val + 1.2, i, f"{val:.1f}%", ha="left", va="center", color="#374151", fontsize=7)
        left += np.array(vals)
    ax_d.set_xlim(0, 100)
    ax_d.set_yticks(y_pos)
    ax_d.set_yticklabels(["C9 oil-storage", "C6 senescence", "C10 senescence"])
    ax_d.invert_yaxis()
    ax_d.set_xlabel("Cluster composition at 185 d (%)")
    ax_d.set_title("Mature cell-state composition", fontweight="bold", pad=4)
    ax_d.legend(frameon=False, loc="lower center", bbox_to_anchor=(0.5, -0.38), ncol=2)
    ax_d.grid(axis="x", color="#e5e7eb", linewidth=0.6)
    ax_d.set_axisbelow(True)
    ax_d.spines[["top", "right", "left"]].set_visible(False)
    add_panel_label(ax_d, "d")

    fad2 = fad2.sort_values(["variety", "timepoint"])
    for variety in ["FL", "TN"]:
        sub = fad2[fad2["variety"] == variety].copy()
        sub["timepoint_num"] = sub["timepoint"].astype(int)
        ax_e.plot(
            sub["timepoint_num"],
            sub["fad2_detection_fraction"] * 100,
            marker="o",
            linewidth=1.8,
            markersize=4,
            color=VAR_COLORS[variety],
            label=VAR_LABELS[variety],
        )
        for _, row in sub.iterrows():
            if int(row["timepoint"]) in (125, 185):
                if variety == "TN" and int(row["timepoint"]) == 185:
                    offset = 1.0
                elif variety == "TN" and int(row["timepoint"]) == 125:
                    offset = -1.8
                else:
                    offset = 1.2
                ax_e.text(
                    int(row["timepoint"]),
                    row["fad2_detection_fraction"] * 100 + offset,
                    f"{row['fad2_detection_fraction'] * 100:.1f}%",
                    color=VAR_COLORS[variety],
                    ha="center",
                    va="bottom" if offset > 0 else "top",
                    fontsize=7,
                )
    ax_e.set_xticks([95, 125, 185])
    ax_e.set_xlabel("Developmental stage (d)")
    ax_e.set_ylabel("FAD2+ nuclei (%)")
    ax_e.set_title("FAD2 mRNA detection", fontweight="bold", pad=4)
    ax_e.set_ylim(0, max(24, fad2["fad2_detection_fraction"].max() * 120))
    ax_e.grid(axis="y", color="#e5e7eb", linewidth=0.6)
    ax_e.set_axisbelow(True)
    ax_e.spines[["top", "right"]].set_visible(False)
    ax_e.legend(frameon=False, loc="upper left")
    add_panel_label(ax_e, "e")

    fig.suptitle(
        "Single-nucleus resolution of mature oil-storage and senescence-associated states",
        x=0.5,
        y=0.985,
        fontsize=11,
        fontweight="bold",
    )
    fig.subplots_adjust(left=0.11, right=0.985, top=0.91, bottom=0.13)

    for ext in ("png", "pdf", "svg"):
        out = OUT_PREFIX.with_suffix(f".{ext}")
        if ext == "png":
            fig.savefig(out, dpi=600)
        else:
            fig.savefig(out)
    plt.close(fig)


if __name__ == "__main__":
    main()
