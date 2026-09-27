#!/usr/bin/env python3
"""Build source tables and publication panels Fig. 4d-f without pandas."""

from __future__ import annotations

import csv
import json
import math
import os
from collections import Counter
from pathlib import Path

os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
os.environ.setdefault("MPLBACKEND", "Agg")
os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib-fig4-pan39")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import numpy as np


RUN = Path(__file__).resolve().parents[1]
PAN = Path("${ANALYSIS_DIR}/14_pan_genome/11_new_pan")
OLD_NX = Path(
    "${ANALYSIS_DIR}/20_results/"
    "Figure1/09_new_figure/02_N50/contig_Nx_stats_full.tsv"
)
PAV = PAN / "results/postprocess_attempt1/Orthogroups.PAV.material33.tsv"

NX_COLS = [f"N{x}" for x in range(10, 101, 10)]
PAIR_COLORS = {
    "seedless": "#E65D6D",
    "bk": "#BBA8D8",
    "dura": "#4667AE",
    "pisifera": "#F2B76A",
    "nrly": "#29AFD4",
    "meizhou4": "#A7D9DD",
}
PAIR_DISPLAY = {
    "seedless": "FL",
    "bk": "BK",
    "dura": "Dura",
    "pisifera": "Pisifera",
    "nrly": "NRLY",
    "meizhou4": "Meizhou4",
}
SAMPLE_META = {
    "American_hap1": ("seedless", "hap1"),
    "Africa_hap2": ("seedless", "hap2"),
    "bk_hap1": ("bk", "hap1"),
    "bk_hap2": ("bk", "hap2"),
    "dura_hap1": ("dura", "hap1"),
    "dura_hap2": ("dura", "hap2"),
    "pisifera_hap1": ("pisifera", "hap1"),
    "pisifera_hap2": ("pisifera", "hap2"),
    "nrly_hap1": ("nrly", "hap1"),
    "nrly_hap2": ("nrly", "hap2"),
    "meizhou4_hap1": ("meizhou4", "hap1"),
    "meizhou4_hap2": ("meizhou4", "hap2"),
}
NEW_SAMPLES = [
    "dura_hap1", "dura_hap2", "pisifera_hap1", "pisifera_hap2",
    "nrly_hap1", "nrly_hap2", "meizhou4_hap1", "meizhou4_hap2",
]
OLD_NAME_MAP = {
    "seedless_hap1": "American_hap1",
    "seedless_hap2": "Africa_hap2",
    "tenera_hap1": "bk_hap1",
    "tenera_hap2": "bk_hap2",
    "houke": "houke_pa_genome",
}
DROP_OLD = {"dura", "pisifera", "niriliya"}

F_COLORS = {
    "Core": "#FBA296",
    "Soft-core": "#84C3B8",
    "Shell": "#7BC9F3",
    "Cloud": "#FDE5B0",
}
CLASS_ORDER = ("Core", "Soft-core", "Shell", "Cloud")


def configure_plotting() -> None:
    plt.rcParams.update(
        {
            "font.family": "Liberation Sans",
            "font.size": 10,
            "axes.labelsize": 12,
            "axes.linewidth": 0.7,
            "xtick.labelsize": 9,
            "ytick.labelsize": 9,
            "xtick.major.width": 0.6,
            "ytick.major.width": 0.6,
            "xtick.major.size": 3.0,
            "ytick.major.size": 3.0,
            "legend.fontsize": 8,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
        }
    )


def clean_axis(ax) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_color("#808080")
    ax.spines["bottom"].set_color("#808080")
    ax.tick_params(colors="#202020")


def save_figure(fig, stem: str, dpi: int = 600) -> None:
    out = RUN / "figures"
    fig.savefig(out / f"{stem}.pdf", facecolor="white")
    fig.savefig(out / f"{stem}.svg", facecolor="white")
    fig.savefig(out / f"{stem}_600dpi.png", dpi=dpi, facecolor="white")
    plt.close(fig)


def nx_from_fai(path: Path) -> dict[str, int]:
    lengths = []
    with path.open() as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 2:
                lengths.append(int(fields[1]))
    if not lengths:
        raise ValueError(f"No sequences in {path}")
    lengths.sort(reverse=True)
    total = sum(lengths)
    result: dict[str, int] = {}
    cumulative = 0
    target_idx = 0
    targets = [total * x / 100.0 for x in range(10, 101, 10)]
    for length in lengths:
        cumulative += length
        while target_idx < len(targets) and cumulative >= targets[target_idx]:
            result[NX_COLS[target_idx]] = length
            target_idx += 1
    for name in NX_COLS:
        result.setdefault(name, lengths[-1])
    result["Total"] = total
    result["N_contigs"] = len(lengths)
    result["Largest"] = lengths[0]
    return result


def build_nx_table() -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    with OLD_NX.open() as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            old = row["Sample"]
            old_key = old[1:] if old.startswith("Y") and old[1:].isdigit() else old
            if old_key in DROP_OLD:
                continue
            sample = OLD_NAME_MAP.get(old_key, old_key)
            rows.append(
                {"Sample": sample, **{k: int(row[k]) for k in NX_COLS},
                 "Total": int(row["Total"]), "N_contigs": int(row["N_contigs"]),
                 "Largest": int(row["Largest"]), "Source": "accepted_old_Nx"}
            )
    for sample in NEW_SAMPLES:
        fai = PAN / "inputs/derived_new_haplotypes" / sample / f"{sample}.ragtag.fasta.fai"
        rows.append({"Sample": sample, **nx_from_fai(fai), "Source": str(fai)})
    if len(rows) != 39 or len({str(r["Sample"]) for r in rows}) != 39:
        raise AssertionError(f"Nx rows are not exactly 39: {len(rows)}")
    path = RUN / "tables/Fig4d_Nx_pan39.tsv"
    with path.open("w", newline="") as handle:
        fields = ["Sample", *NX_COLS, "Total", "N_contigs", "Largest", "Source"]
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    return rows


def plot_d(rows: list[dict[str, object]], labelled: bool) -> None:
    x = np.arange(10, 101, 10)
    y_all = np.asarray([[float(r[k]) / 1e6 for k in NX_COLS] for r in rows])
    hi = {str(r["Sample"]): r for r in rows if str(r["Sample"]) in SAMPLE_META}
    bg = [r for r in rows if str(r["Sample"]) not in SAMPLE_META]
    if len(hi) != 12 or len(bg) != 27:
        raise AssertionError(f"Expected 12 highlighted and 27 background rows: {len(hi)}, {len(bg)}")
    fig, ax = plt.subplots(figsize=(9, 6.5))
    med = np.percentile(y_all, 50, axis=0)
    for row in bg:
        ax.plot(x, [float(row[k]) / 1e6 for k in NX_COLS], color="#B3B3B3",
                lw=0.55, alpha=0.48, zorder=1)
    ax.plot(x, med, color="#666666", lw=1.1, ls=(0, (4, 3)), zorder=2)
    for sample, row in hi.items():
        material, hap = SAMPLE_META[sample]
        style = "-" if hap == "hap1" else "--"
        marker = "o" if hap == "hap1" else "s"
        ax.plot(x, [float(row[k]) / 1e6 for k in NX_COLS], color=PAIR_COLORS[material],
                lw=1.65, ls=style, marker=marker, ms=3.1, alpha=0.94, zorder=4)
    colour_handles = [Line2D([0], [0], color=PAIR_COLORS[m], lw=2.0,
                             label=PAIR_DISPLAY[m]) for m in PAIR_COLORS]
    hap_handles = [Line2D([0], [0], color="#555555", lw=1.5, ls="-", marker="o",
                          ms=3, label="hap1"),
                   Line2D([0], [0], color="#555555", lw=1.5, ls="--", marker="s",
                          ms=3, label="hap2")]
    leg1 = ax.legend(handles=colour_handles, loc="upper right", bbox_to_anchor=(0.99, 0.99),
                     frameon=False, handlelength=2.2, borderaxespad=0.1)
    ax.add_artist(leg1)
    ax.legend(handles=hap_handles, loc="upper right", bbox_to_anchor=(0.99, 0.66),
              frameon=False, handlelength=2.2, borderaxespad=0.1)
    ax.set_xlim(8, 102)
    ax.set_xticks(x)
    ax.set_xlabel(r"$N_x$ (%)")
    ax.set_ylabel("Contig length (Mb)")
    ax.set_ylim(bottom=0)
    clean_axis(ax)
    fig.subplots_adjust(left=0.12, right=0.98, bottom=0.14, top=0.96)
    if labelled:
        fig.text(0.012, 0.985, "d", ha="left", va="top", fontsize=17, fontweight="bold")
    save_figure(fig, f"Fig4d_Nx_pan39_{'labelled' if labelled else 'no_label'}")


def read_pav() -> tuple[list[str], np.ndarray, list[str]]:
    ogs: list[str] = []
    matrix: list[list[int]] = []
    with PAV.open() as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader)
        materials = header[1:]
        for row in reader:
            if not row:
                continue
            ogs.append(row[0])
            matrix.append([int(v) for v in row[1:]])
    data = np.asarray(matrix, dtype=np.bool_)
    if data.shape != (48920, 33):
        raise AssertionError(f"Unexpected material PAV shape: {data.shape}")
    return materials, data, ogs


def family_class(freq: int) -> str:
    if freq == 33:
        return "Core"
    if freq == 32:
        return "Soft-core"
    if freq >= 2:
        return "Shell"
    return "Cloud"


def build_pan_core(data: np.ndarray, permutations: int = 1000) -> np.ndarray:
    rng = np.random.default_rng(42)
    n_materials = data.shape[1]
    pans = np.zeros((permutations, n_materials), dtype=np.int32)
    cores = np.zeros((permutations, n_materials), dtype=np.int32)
    for rep in range(permutations):
        order = rng.permutation(n_materials)
        union = np.zeros(data.shape[0], dtype=np.bool_)
        intersection = np.ones(data.shape[0], dtype=np.bool_)
        for i, idx in enumerate(order):
            union |= data[:, idx]
            intersection &= data[:, idx]
            pans[rep, i] = int(union.sum())
            cores[rep, i] = int(intersection.sum())
    summary = np.column_stack(
        [np.arange(1, n_materials + 1), pans.mean(axis=0), pans.std(axis=0, ddof=1),
         cores.mean(axis=0), cores.std(axis=0, ddof=1)]
    )
    path = RUN / "tables/Fig4e_pan_core_material33_1000perm_seed42.tsv"
    np.savetxt(path, summary, delimiter="\t", fmt=["%d", "%.6f", "%.6f", "%.6f", "%.6f"],
               header="N_materials\tPan_mean\tPan_sd\tCore_mean\tCore_sd", comments="")
    return summary


def plot_e(summary: np.ndarray, labelled: bool) -> None:
    n, pan, pan_sd, core, core_sd = summary.T
    fig, ax = plt.subplots(figsize=(10, 6))
    pan_c, core_c = "#E65D6D", "#29AFD4"
    ax.errorbar(n, pan / 1000, yerr=pan_sd / 1000, color=pan_c, lw=1.5,
                elinewidth=0.9, capsize=0, marker="o", ms=4.4, mec=pan_c, mfc=pan_c,
                label="Pan-genome")
    ax.errorbar(n, core / 1000, yerr=core_sd / 1000, color=core_c, lw=1.5,
                elinewidth=0.9, capsize=0, marker="s", ms=4.2, mec=core_c, mfc=core_c,
                label="Core-genome")
    ax.text(n[-1] + 1.0, pan[-1] / 1000, f"{pan[-1] / 1000:.1f}k", color=pan_c,
            fontsize=9, fontweight="bold", va="center")
    ax.text(n[-1] + 1.0, core[-1] / 1000, f"{core[-1] / 1000:.1f}k", color=core_c,
            fontsize=9, fontweight="bold", va="center")
    ax.set_xlabel("Number of genomes")
    ax.set_ylabel(r"Gene families ($\times 10^3$)")
    ax.set_xlim(0.2, 37.5)
    ax.set_xticks(np.arange(5, 36, 5))
    ax.set_ylim(bottom=18)
    ax.legend(frameon=False, loc="upper right", bbox_to_anchor=(0.86, 0.87))
    clean_axis(ax)
    fig.subplots_adjust(left=0.13, right=0.97, bottom=0.15, top=0.96)
    if labelled:
        fig.text(0.012, 0.985, "e", ha="left", va="top", fontsize=17, fontweight="bold")
    save_figure(fig, f"Fig4e_pan_core_material33_{'labelled' if labelled else 'no_label'}")


def build_frequency(data: np.ndarray, ogs: list[str]) -> tuple[np.ndarray, dict[str, int]]:
    freq = data.sum(axis=1).astype(int)
    classes = [family_class(int(x)) for x in freq]
    class_counts = Counter(classes)
    if sum(class_counts.values()) != 48920:
        raise AssertionError("Family class counts do not sum to 48,920")
    with (RUN / "tables/Fig4f_orthogroup_material_frequency.tsv").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["Orthogroup", "Material_frequency", "Family_class"])
        writer.writerows(zip(ogs, freq.tolist(), classes))
    with (RUN / "tables/Fig4f_class_summary.tsv").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["Family_class", "Orthogroups", "Percentage"])
        for category in CLASS_ORDER:
            count = class_counts[category]
            writer.writerow([category, count, f"{count / 48920 * 100:.4f}"])
    return freq, dict(class_counts)


def plot_f(freq: np.ndarray, class_counts: dict[str, int], labelled: bool) -> None:
    counts = Counter(freq.tolist())
    x = np.arange(1, 34)
    y = np.asarray([counts.get(int(i), 0) for i in x])
    colours = [F_COLORS[family_class(int(i))] for i in x]
    fig, ax = plt.subplots(figsize=(14, 5.5))
    ax.bar(x, y, width=0.84, color=colours, edgecolor="white", linewidth=0.25, zorder=1)
    ax.set_xlim(0.35, 33.65)
    ax.set_ylim(bottom=0)
    ax.set_xticks(x)
    ax.tick_params(axis="x", labelsize=7)
    ax.set_xlabel("Number of genomes")
    ax.set_ylabel("Number of gene families")
    clean_axis(ax)

    inset = fig.add_axes([0.49, 0.18, 0.32, 0.70], facecolor="white")
    values = [class_counts[c] for c in CLASS_ORDER]
    labels = [f"{c}\n{class_counts[c]:,} ({class_counts[c] / 48920 * 100:.1f}%)"
              for c in CLASS_ORDER]
    inset.pie(values, labels=labels, colors=[F_COLORS[c] for c in CLASS_ORDER],
              startangle=92, counterclock=False, wedgeprops={"width": 0.42, "edgecolor": "white", "linewidth": 1.1},
              labeldistance=1.08, textprops={"fontsize": 8, "color": "#202020"})
    inset.text(0, 0.05, f"{sum(values):,}", ha="center", va="center",
               fontsize=11, fontweight="bold")
    inset.text(0, -0.11, "gene families", ha="center", va="center", fontsize=8, color="#777777")
    inset.set_axis_off()
    handles = [Patch(facecolor=F_COLORS[c], edgecolor="none", label=c) for c in CLASS_ORDER]
    ax.legend(handles=handles, frameon=False, loc="upper right", bbox_to_anchor=(0.985, 0.985),
              handlelength=1.0, handletextpad=0.4)
    fig.subplots_adjust(left=0.09, right=0.985, bottom=0.17, top=0.96)
    if labelled:
        fig.text(0.012, 0.985, "f", ha="left", va="top", fontsize=17, fontweight="bold")
    save_figure(fig, f"Fig4f_frequency_donut_material33_{'labelled' if labelled else 'no_label'}")


def main() -> None:
    configure_plotting()
    rows = build_nx_table()
    materials, data, ogs = read_pav()
    summary = build_pan_core(data)
    freq, class_counts = build_frequency(data, ogs)
    for labelled in (False, True):
        plot_d(rows, labelled)
        plot_e(summary, labelled)
        plot_f(freq, class_counts, labelled)
    audit = {
        "nx_rows": len(rows),
        "highlighted_haplotypes": len(SAMPLE_META),
        "background_haplotypes": len(rows) - len(SAMPLE_META),
        "pav_shape": list(data.shape),
        "materials": materials,
        "orthogroups": len(ogs),
        "class_counts": {c: class_counts[c] for c in CLASS_ORDER},
        "permutations": 1000,
        "random_seed": 42,
        "status": "PASS",
    }
    (RUN / "provenance/Fig4def_audit.json").write_text(json.dumps(audit, indent=2) + "\n")
    print(json.dumps(audit, indent=2))


if __name__ == "__main__":
    main()
