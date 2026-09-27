#!/usr/bin/env python3
"""Rebuild Fig. 4e/f from the 64,577-family, 39-genome matrix.

The older panels collapsed six phased pairs to 33 materials.  This script
keeps all 39 OrthoFinder genome columns, recomputes frequency classes and
pan/core permutation summaries, and exports publication-ready figures.
"""

from __future__ import annotations

import csv
import hashlib
import json
import os
from pathlib import Path

os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
os.environ.setdefault("MPLBACKEND", "Agg")
os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib-fig4ef-genome39")

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Patch
from scipy.interpolate import CubicSpline, PchipInterpolator


PARENT = Path("${ANALYSIS_DIR}/22_answer_reviews/00_ms/03_V3/04_figure4/Fig4_d_i_pan39_material33_singletons_20260811")
RUN = PARENT / "pan39_genome_redraw_20260812"
FIGURES = RUN / "figures"
TABLES = RUN / "tables"
PROV = RUN / "provenance"
MEMBERS = PARENT / "tables/Orthogroups.members.GO_singletons_64577.tsv"
N_PERMUTATIONS = 1000
SEED = 42
N_GENOMES = 39
N_FAMILIES = 64577

COLORS = {
    "Pan": "#EC5B70",
    "CoreCurve": "#22A9C7",
    "Core": "#F69A8D",
    "Soft-core": "#7EC3B6",
    "Shell": "#71BCE6",
    "Cloud": "#F8DFA5",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def configure_style() -> None:
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Liberation Sans", "DejaVu Sans"],
        "font.size": 8.0,
        "axes.labelsize": 9.0,
        "axes.linewidth": 0.75,
        "xtick.labelsize": 7.4,
        "ytick.labelsize": 7.4,
        "legend.fontsize": 7.2,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "svg.fonttype": "none",
    })


def read_members_as_pav(path: Path):
    """Read only occupancy, avoiding retention of large member strings."""
    families = []
    rows = []
    with path.open(newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader)
        genomes = header[1:]
        if len(genomes) != N_GENOMES or len(set(genomes)) != N_GENOMES:
            raise AssertionError(f"Expected 39 unique genomes; got {len(genomes)}")
        for fields in reader:
            if len(fields) < len(header):
                fields.extend([""] * (len(header) - len(fields)))
            if len(fields) != len(header):
                raise AssertionError(f"Malformed row {fields[0]}: {len(fields)} fields")
            families.append(fields[0])
            rows.append([bool(value.strip()) for value in fields[1:]])
    pav = np.asarray(rows, dtype=bool)
    if pav.shape != (N_FAMILIES, N_GENOMES):
        raise AssertionError(f"Expected {(N_FAMILIES, N_GENOMES)}, got {pav.shape}")
    if len(set(families)) != N_FAMILIES:
        raise AssertionError("Orthogroup/family IDs are not unique")
    if np.any(pav.sum(axis=1) == 0):
        raise AssertionError("At least one family has no genome occupancy")
    return families, genomes, pav


def permutation_summary(pav: np.ndarray) -> pd.DataFrame:
    rng = np.random.default_rng(SEED)
    pan_values = np.empty((N_PERMUTATIONS, N_GENOMES), dtype=np.int32)
    core_values = np.empty_like(pan_values)
    for index in range(N_PERMUTATIONS):
        order = rng.permutation(N_GENOMES)
        ordered = pav[:, order]
        pan_values[index] = np.logical_or.accumulate(ordered, axis=1).sum(axis=0)
        core_values[index] = np.logical_and.accumulate(ordered, axis=1).sum(axis=0)
    ddof = 1
    return pd.DataFrame({
        "N_genomes": np.arange(1, N_GENOMES + 1),
        "Pan_mean": pan_values.mean(axis=0),
        "Pan_sd": pan_values.std(axis=0, ddof=ddof),
        "Core_mean": core_values.mean(axis=0),
        "Core_sd": core_values.std(axis=0, ddof=ddof),
    })


class FlatBlend:
    """C2 body plus quintic tail ending at the observed endpoint, flat."""

    def __init__(self, x, y, transition_start, endpoint):
        self.x0 = float(transition_start)
        self.x1 = float(N_GENOMES)
        self.length = self.x1 - self.x0
        base_knots = np.asarray([1, 2, 3, 4, 6, 8, 11, 14, 17, 20, 23, 26, 29, 32, 35], float)
        knots = np.r_[base_knots[base_knots < self.x0], self.x0]
        values = np.interp(knots, x, y)
        self.body = CubicSpline(knots, values, bc_type=((1, y[1] - y[0]), (2, 0.0)))
        f0 = float(self.body(self.x0))
        d0 = float(self.body(self.x0, 1))
        dd0 = float(self.body(self.x0, 2))
        c0, c1, c2 = f0, d0 * self.length, dd0 * self.length ** 2 / 2
        rhs = np.asarray([endpoint - (c0 + c1 + c2), -(c1 + 2 * c2), -2 * c2])
        system = np.asarray([[1, 1, 1], [3, 4, 5], [6, 12, 20]], float)
        self.coef = np.r_[c0, c1, c2, np.linalg.solve(system, rhs)]

    def __call__(self, values, derivative=0):
        values = np.asarray(values, float)
        output = np.empty_like(values)
        body = values <= self.x0
        output[body] = self.body(values[body], derivative)
        s = (values[~body] - self.x0) / self.length
        coef = self.coef.copy()
        for _ in range(derivative):
            coef = np.asarray([i * coef[i] for i in range(1, len(coef))]) / self.length
        output[~body] = sum(coef[i] * s ** i for i in range(len(coef)))
        return output


def choose_curves(summary: pd.DataFrame):
    x = summary.N_genomes.to_numpy(float)
    pan = summary.Pan_mean.to_numpy(float)
    core = summary.Core_mean.to_numpy(float)
    pan_sd = summary.Pan_sd.to_numpy(float)
    core_sd = summary.Core_sd.to_numpy(float)
    grid = np.linspace(1, N_GENOMES, 3801)
    endpoint_pan = float(pan[-1])
    endpoint_core = float(core[-1])
    candidates = []
    selected = None
    for transition in range(28, 37):
        pan_curve = FlatBlend(x, pan, transition, endpoint_pan)
        core_curve = FlatBlend(x, core, transition, endpoint_core)
        pan_pred, core_pred = pan_curve(x), core_curve(x)
        pan_z = np.max(np.abs(pan_pred[:-1] - pan[:-1]) / np.maximum(pan_sd[:-1], 1))
        core_z = np.max(np.abs(core_pred[:-1] - core[:-1]) / np.maximum(core_sd[:-1], 1))
        monotone = bool(np.min(pan_curve(grid, 1)) >= -1e-6 and np.max(core_curve(grid, 1)) <= 1e-6)
        row = {"transition_start": transition, "pan_max_SD": float(pan_z),
               "core_max_SD": float(core_z), "monotone": monotone}
        candidates.append(row)
        if selected is None and monotone and pan_z <= 1 and core_z <= 1:
            selected = (transition, pan_curve, core_curve)
    if selected is None:
        # A fallback that remains monotone and is still directly audited.
        valid = [r for r in candidates if r["monotone"]]
        if not valid:
            raise AssertionError(f"No monotone flat-end curve: {candidates}")
        best = min(valid, key=lambda r: max(r["pan_max_SD"], r["core_max_SD"]))
        selected = (best["transition_start"], FlatBlend(x, pan, best["transition_start"], endpoint_pan),
                    FlatBlend(x, core, best["transition_start"], endpoint_core))
    return grid, selected, candidates


def save(fig, stem):
    for suffix, kwargs in (("pdf", {}), ("svg", {}), ("png", {"dpi": 600})):
        name = f"{stem}_600dpi.png" if suffix == "png" else f"{stem}.{suffix}"
        fig.savefig(FIGURES / name, facecolor="white", **kwargs)
    plt.close(fig)


def format_k(value: float) -> str:
    return f"{value / 1000:.1f}k"


def draw_e(summary, grid, pan_curve, core_curve):
    x = summary.N_genomes.to_numpy(float)
    pan = summary.Pan_mean.to_numpy(float)
    core = summary.Core_mean.to_numpy(float)
    pan_sd = summary.Pan_sd.to_numpy(float)
    core_sd = summary.Core_sd.to_numpy(float)
    pan_fit, core_fit = pan_curve(grid), core_curve(grid)
    pan_band = np.clip(PchipInterpolator(x, pan_sd)(grid), 0, None)
    core_band = np.clip(PchipInterpolator(x, core_sd)(grid), 0, None)
    for labelled in (False, True):
        # Closely matches the original Fig4-7.8 panel-e aspect ratio (~1.6:1).
        fig, ax = plt.subplots(figsize=(5.20, 3.20), facecolor="white")
        ax.fill_between(grid, (pan_fit-pan_band)/1000, (pan_fit+pan_band)/1000,
                        color=COLORS["Pan"], alpha=.10, linewidth=0)
        ax.fill_between(grid, (core_fit-core_band)/1000, (core_fit+core_band)/1000,
                        color=COLORS["CoreCurve"], alpha=.10, linewidth=0)
        for values, errors, color in ((pan, pan_sd, COLORS["Pan"]),
                                      (core, core_sd, COLORS["CoreCurve"])):
            ax.errorbar(x, values/1000, yerr=errors/1000, fmt="o", ms=2.15,
                        mfc="white", mec=color, mew=.55, ecolor=color,
                        elinewidth=.55, capsize=1.25, capthick=.55, alpha=.45, zorder=2)
        for values, color in ((pan_fit, COLORS["Pan"]), (core_fit, COLORS["CoreCurve"])):
            ax.plot(grid, values/1000, color=color, lw=5.2, alpha=.13,
                    solid_capstyle="round", zorder=2.6)
            ax.plot(grid, values/1000, color=color, lw=2.55,
                    solid_capstyle="round", zorder=3)
        ax.text(38.2, pan[-1]/1000 + 1.0, f"Pan-genome  {format_k(pan[-1])}",
                color=COLORS["Pan"], fontsize=8.6, fontweight="bold", ha="right")
        ax.text(38.2, core[-1]/1000 - 1.0, f"Core-genome  {format_k(core[-1])}",
                color=COLORS["CoreCurve"], fontsize=8.6, fontweight="bold", ha="right", va="top")
        ax.set_xlabel("Number of genomes")
        ax.set_ylabel(r"Gene families ($\times 10^3$)")
        ax.set_xlim(.2, 40.5)
        ax.set_xticks([5, 10, 15, 20, 25, 30, 35, 39])
        upper = np.ceil((pan[-1]/1000 + 5) / 5) * 5
        lower = max(0, np.floor((core[-1]/1000 - 5) / 5) * 5)
        ax.set_ylim(lower, upper)
        ax.grid(axis="y", color="#D9D9D9", lw=.48, alpha=.28)
        ax.spines[["top", "right"]].set_visible(False)
        ax.spines[["left", "bottom"]].set_color("#707070")
        ax.tick_params(width=.7, length=3, colors="#202020")
        fig.subplots_adjust(left=.15, right=.975, bottom=.17, top=.96)
        if labelled:
            fig.text(.012, .985, "e", ha="left", va="top", fontsize=17, fontweight="bold")
        save(fig, f"Fig4e_pan_core_genome39_{'labelled' if labelled else 'no_label'}")


def draw_f(freq_counts):
    frequency = freq_counts.set_index("Genome_frequency")["Number_of_families"].reindex(range(1, 40), fill_value=0)
    class_order = ["Core", "Soft-core", "Shell", "Cloud"]
    class_counts = freq_counts.groupby("Class")["Number_of_families"].sum().reindex(class_order)
    bar_colors = [COLORS["Cloud"] if i == 1 else COLORS["Soft-core"] if i == 38
                  else COLORS["Core"] if i == 39 else COLORS["Shell"] for i in range(1, 40)]
    for labelled in (False, True):
        # Matches the original panel-f wide layout (~1.9:1).
        fig, ax = plt.subplots(figsize=(6.05, 3.20), facecolor="white")
        ax.bar(np.arange(1, 40), frequency.values, color=bar_colors, width=.82, edgecolor="none")
        ax.set_xlabel("Number of genomes")
        ax.set_ylabel("Number of gene families")
        ax.set_xlim(.25, 39.75)
        ax.set_xticks(np.arange(1, 40))
        ax.tick_params(axis="x", labelsize=5.1, length=2.2, pad=1.5)
        ax.tick_params(axis="y", length=3)
        ax.spines[["top", "right"]].set_visible(False)
        ax.spines[["left", "bottom"]].set_color("#707070")
        ax.grid(False)

        inset = ax.inset_axes([.43, .20, .39, .70])
        wedges, _ = inset.pie(class_counts.values, startangle=90, counterclock=False,
                              colors=[COLORS[name] for name in class_order],
                              wedgeprops={"width": .42, "edgecolor": "white", "linewidth": .8})
        inset.text(0, .05, f"{N_FAMILIES:,}", ha="center", va="center", fontsize=9.3, fontweight="bold")
        inset.text(0, -.16, "gene families", ha="center", va="center", fontsize=6.8, color="#777777")
        inset.set_aspect("equal")
        for wedge, name, count in zip(wedges, class_order, class_counts.values):
            angle = np.deg2rad((wedge.theta1 + wedge.theta2) / 2)
            radius = 1.07
            x0, y0 = radius*np.cos(angle), radius*np.sin(angle)
            ha = "left" if x0 >= 0 else "right"
            pct = 100 * count / N_FAMILIES
            inset.text(x0, y0, f"{name}\n{int(count):,} ({pct:.1f}%)", ha=ha, va="center", fontsize=6.5)
        ax.legend(handles=[Patch(facecolor=COLORS[name], label=name) for name in class_order],
                  frameon=False, loc="upper right", borderaxespad=.1, handlelength=.9,
                  handletextpad=.4, labelspacing=.25)
        fig.subplots_adjust(left=.10, right=.99, bottom=.18, top=.96)
        if labelled:
            fig.text(.012, .985, "f", ha="left", va="top", fontsize=17, fontweight="bold")
        save(fig, f"Fig4f_frequency_donut_genome39_{'labelled' if labelled else 'no_label'}")


def main():
    for directory in (FIGURES, TABLES, PROV):
        directory.mkdir(parents=True, exist_ok=True)
    configure_style()
    families, genomes, pav = read_members_as_pav(MEMBERS)
    freq = pav.sum(axis=1)
    classes = np.select([freq == 39, freq == 38, (freq >= 2) & (freq <= 37), freq == 1],
                        ["Core", "Soft-core", "Shell", "Cloud"], default="ERROR")
    if np.any(classes == "ERROR"):
        raise AssertionError("Frequency classification failed")

    pd.DataFrame(pav.astype(np.uint8), index=families, columns=genomes).rename_axis("Orthogroup").to_csv(
        TABLES / "Orthogroups.PAV.genome39.GO_singletons_64577.tsv", sep="\t")
    freq_counts = (pd.DataFrame({"Genome_frequency": freq, "Class": classes})
                   .value_counts(sort=False).rename("Number_of_families").reset_index()
                   .sort_values("Genome_frequency"))
    freq_counts.to_csv(TABLES / "Orthogroups.genome39.frequency_class_counts.GO64577.tsv", sep="\t", index=False)
    class_counts = (pd.DataFrame({"Class": classes}).value_counts(sort=False)
                    .rename("Number_of_families").reset_index())
    class_counts["Percent"] = 100 * class_counts.Number_of_families / N_FAMILIES
    class_counts.to_csv(TABLES / "Orthogroups.genome39.class_counts.GO64577.tsv", sep="\t", index=False)

    summary = permutation_summary(pav)
    if not np.isclose(summary.iloc[-1].Pan_mean, N_FAMILIES) or summary.iloc[-1].Pan_sd != 0:
        raise AssertionError("Pan endpoint is invalid")
    core39 = int(np.sum(freq == 39))
    if not np.isclose(summary.iloc[-1].Core_mean, core39) or summary.iloc[-1].Core_sd != 0:
        raise AssertionError("Core endpoint is invalid")
    summary.to_csv(TABLES / "Fig4e_pan_core_genome39_GO64577_1000perm_seed42.tsv", sep="\t", index=False,
                   float_format="%.6f")

    grid, (transition, pan_curve, core_curve), candidates = choose_curves(summary)
    draw_e(summary, grid, pan_curve, core_curve)
    draw_f(freq_counts)

    audit = {
        "status": "PASS",
        "input": str(MEMBERS),
        "input_sha256": sha256(MEMBERS),
        "dimensions": {"families": int(pav.shape[0]), "genomes": int(pav.shape[1])},
        "genomes": genomes,
        "permutations": N_PERMUTATIONS,
        "seed": SEED,
        "class_definition": {"Core": "39/39", "Soft-core": "38/39", "Shell": "2-37/39", "Cloud": "1/39"},
        "class_counts": {row.Class: int(row.Number_of_families) for row in class_counts.itertuples()},
        "class_count_sum": int(class_counts.Number_of_families.sum()),
        "pan_endpoint": float(summary.iloc[-1].Pan_mean),
        "core_endpoint": float(summary.iloc[-1].Core_mean),
        "flat_blend_transition_start": int(transition),
        "terminal_derivatives": {
            "Pan_first": float(pan_curve(np.asarray([39.]), 1)[0]),
            "Pan_second": float(pan_curve(np.asarray([39.]), 2)[0]),
            "Core_first": float(core_curve(np.asarray([39.]), 1)[0]),
            "Core_second": float(core_curve(np.asarray([39.]), 2)[0]),
        },
        "curve_candidates": candidates,
        "interpretation_boundary": "Smooth curves are visual guides constrained to the observed n=39 endpoints; all 39 permutation means and SD bars are retained.",
        "figure_sizes_inches": {"e": [5.20, 3.20], "f": [6.05, 3.20]},
    }
    if audit["class_count_sum"] != N_FAMILIES:
        raise AssertionError("Class counts do not sum to 64,577")
    (PROV / "AUDIT.json").write_text(json.dumps(audit, indent=2) + "\n")
    (PROV / "README.txt").write_text(
        "Fig4e/f rebuilt at the 39-genome level from the unchanged 64,577-family universe.\n"
        "Old 33-material files were not overwritten.\n"
        "Core=39; Soft-core=38; Shell=2-37; Cloud=1 genome.\n"
        "Curves: 1,000 random genome orders (seed 42), means +/- sample SD; C2 flat-end visual smoother.\n"
    )
    print(json.dumps(audit, indent=2))


if __name__ == "__main__":
    main()
