#!/usr/bin/env python3
"""Back up current figures and redraw panels d-h for the GO64577 universe."""

from __future__ import annotations

import csv
import hashlib
import importlib.util
import json
import os
import platform
import shutil
import sys
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
os.environ.setdefault("MPLBACKEND", "Agg")
os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib-fig4-go64577-full")

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import numpy as np
import pandas as pd
from scipy.interpolate import PchipInterpolator
from scipy.optimize import curve_fit

RUN = Path("${ANALYSIS_DIR}/22_answer_reviews/00_ms/03_V3/04_figure4/Fig4_d_i_pan39_material33_singletons_20260811")
SOURCE = Path("${ANALYSIS_DIR}/22_answer_reviews/00_ms/03_V3/04_figure4/Fig4_d_i_pan39_material33_redraw_20260808")
FIGURES = RUN / "figures"
TABLES = RUN / "tables"
PROV = RUN / "provenance/GO64577_full_redraw"
BACKUP = RUN / "backup/figures_before_GO64577_full_redraw_20260811"
FILTERED_SUBRUN = RUN / "GO_filtered_64577"
SOURCE_PAV = FILTERED_SUBRUN / "tables/Orthogroups.PAV.material33.GO_singletons_64577.tsv"
PAV = TABLES / "Orthogroups.PAV.material33.GO_singletons_64577.tsv"
MEMBERS_SOURCE = TABLES / "Orthogroups.members.with_singletons.tsv"
MEMBERS = TABLES / "Orthogroups.members.GO_singletons_64577.tsv"
CLASS_FILE = TABLES / "Orthogroups.material33.freq_class.GO_singletons_64577.tsv"
CURVE_SOURCE = FILTERED_SUBRUN / "tables/Fig4e_pan_core_material33_GO64577_1000perm_seed42.tsv"
CURVE = TABLES / "Fig4e_pan_core_material33_GO64577_1000perm_seed42.tsv"
TOTAL = 64577
ASSIGNED = 48920
EXPECTED_CLASSES = {"Core": 20736, "Soft-core": 2057, "Shell": 23564, "Cloud": 18220}
CLASS_ORDER = ("Core", "Soft-core", "Shell", "Cloud")
F_COLORS = {"Core": "#FBA296", "Soft-core": "#84C3B8", "Shell": "#7BC9F3", "Cloud": "#FDE5B0"}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def import_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise ImportError(path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def family_class(freq: int) -> str:
    if freq == 33:
        return "Core"
    if freq == 32:
        return "Soft-core"
    if freq >= 2:
        return "Shell"
    return "Cloud"


def backup_figures():
    if not BACKUP.exists():
        BACKUP.parent.mkdir(parents=True, exist_ok=True)
        shutil.copytree(FIGURES, BACKUP)
    if not any(BACKUP.rglob("*.pdf")):
        raise AssertionError("Figure backup is empty")


def prepare_filtered_tables():
    shutil.copy2(SOURCE_PAV, PAV)
    shutil.copy2(CURVE_SOURCE, CURVE)
    retained_ids = []
    matrix = []
    with PAV.open() as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader)
        if len(header) != 34:
            raise AssertionError("Filtered PAV must contain 33 materials")
        for row in reader:
            retained_ids.append(row[0])
            matrix.append([int(value) for value in row[1:]])
    data = np.asarray(matrix, dtype=np.bool_)
    if data.shape != (TOTAL, 33) or len(set(retained_ids)) != TOTAL:
        raise AssertionError(f"Unexpected filtered PAV: {data.shape}")
    retained = set(retained_ids)

    member_rows = 0
    with MEMBERS_SOURCE.open() as source, MEMBERS.open("w", newline="") as target:
        reader = csv.reader(source, delimiter="\t")
        writer = csv.writer(target, delimiter="\t", lineterminator="\n")
        writer.writerow(next(reader))
        for row in reader:
            if row and row[0] in retained:
                writer.writerow(row)
                member_rows += 1
    if member_rows != TOTAL:
        raise AssertionError(f"Filtered members contain {member_rows} rows")

    frequencies = data.sum(axis=1).astype(int)
    classes = [family_class(int(value)) for value in frequencies]
    counts = Counter(classes)
    observed = {name: counts[name] for name in CLASS_ORDER}
    if observed != EXPECTED_CLASSES:
        raise AssertionError(f"Unexpected GO64577 class counts: {observed}")
    with CLASS_FILE.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["Orthogroup", "material_frequency", "class", "family_origin"])
        for index, (orthogroup, frequency, category) in enumerate(zip(retained_ids, frequencies, classes)):
            writer.writerow([orthogroup, int(frequency), category,
                             "assigned_orthogroup" if index < ASSIGNED else "GO_supported_unassigned_singleton"])
    return header[1:], retained_ids, data, frequencies, observed


def finite_endpoint_pan(x, a, b, p, q):
    u = np.clip((33.0 - np.asarray(x)) / 32.0, 0, 1)
    return TOTAL - a * np.power(u, p) - b * np.power(u, q)


def finite_endpoint_core(x, a, b, p, q):
    u = np.clip((33.0 - np.asarray(x)) / 32.0, 0, 1)
    return 20736.0 + a * np.power(u, p) + b * np.power(u, q)


def metrics(y, predicted, parameters):
    residual = y - predicted
    sse = float(np.sum(residual ** 2))
    n = len(y)
    return {"R2": float(1 - sse / np.sum((y - y.mean()) ** 2)),
            "RMSE": float(np.sqrt(sse / n)),
            "AIC": float(n * np.log(sse / n) + 2 * parameters)}


def configure_style():
    plt.rcParams.update({
        "font.family": "sans-serif", "font.sans-serif": ["Liberation Sans", "Arial", "DejaVu Sans"],
        "font.size": 8.0, "axes.labelsize": 9.0, "axes.linewidth": 0.7,
        "xtick.labelsize": 7.5, "ytick.labelsize": 7.5,
        "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
    })


def save(fig, stem):
    fig.savefig(FIGURES / f"{stem}.pdf", facecolor="white")
    fig.savefig(FIGURES / f"{stem}.svg", facecolor="white")
    fig.savefig(FIGURES / f"{stem}_600dpi.png", dpi=600, facecolor="white")
    plt.close(fig)


def draw_e():
    frame = pd.read_csv(CURVE, sep="\t")
    x = frame["N_materials"].to_numpy(float)
    pan = frame["Pan_mean"].to_numpy(float)
    pan_sd = frame["Pan_sd"].to_numpy(float)
    core = frame["Core_mean"].to_numpy(float)
    core_sd = frame["Core_sd"].to_numpy(float)
    if not np.array_equal(frame.iloc[-1].to_numpy(float), [33, 64577, 0, 20736, 0]):
        raise AssertionError("Curve endpoint is not exact")
    pan_par, _ = curve_fit(
        finite_endpoint_pan, x, pan, p0=[27000, 6500, 1.4, 18.0],
        bounds=([0, 0, 1.4, 1.4], [1e6, 1e6, 30, 30]), maxfev=1000000,
    )
    core_par, _ = curve_fit(
        finite_endpoint_core, x, core, p0=[4000, 6500, 1.6, 24],
        bounds=([0, 0, 1.001, 1.001], [1e6, 1e6, 30, 30]), maxfev=1000000,
    )
    grid = np.linspace(1, 33, 1800)
    pan_fit = finite_endpoint_pan(grid, *pan_par) / 1000
    core_fit = finite_endpoint_core(grid, *core_par) / 1000
    pan_band = np.clip(PchipInterpolator(x, pan_sd / 1000)(grid), 0, None)
    core_band = np.clip(PchipInterpolator(x, core_sd / 1000)(grid), 0, None)
    pan_c, core_c = "#E65D6D", "#29AFD4"
    configure_style()
    for labelled in (False, True):
        fig, ax = plt.subplots(figsize=(11.0, 5.2), facecolor="white")
        ax.fill_between(grid, pan_fit - pan_band, pan_fit + pan_band, color=pan_c, alpha=0.10, linewidth=0)
        ax.fill_between(grid, core_fit - core_band, core_fit + core_band, color=core_c, alpha=0.10, linewidth=0)
        ax.errorbar(x, pan / 1000, yerr=pan_sd / 1000, fmt="o", ms=2.25, mfc="white", mec=pan_c,
                    mew=0.6, ecolor=pan_c, elinewidth=0.62, capsize=1.45, capthick=0.62,
                    alpha=0.48, zorder=2)
        ax.errorbar(x, core / 1000, yerr=core_sd / 1000, fmt="o", ms=2.25, mfc="white", mec=core_c,
                    mew=0.6, ecolor=core_c, elinewidth=0.62, capsize=1.45, capthick=0.62,
                    alpha=0.48, zorder=2)
        ax.plot(grid, pan_fit, color=pan_c, lw=3.8, solid_capstyle="round", solid_joinstyle="round", zorder=3)
        ax.plot(grid, core_fit, color=core_c, lw=3.8, solid_capstyle="round", solid_joinstyle="round", zorder=3)
        ax.text(34.0, 64.577, "Pan-genome  64.6k", color=pan_c, fontsize=9, fontweight="bold", va="center")
        ax.text(34.0, 20.736, "Core-genome  20.7k", color=core_c, fontsize=9, fontweight="bold", va="center")
        ax.set_xlabel("Number of materials")
        ax.set_ylabel(r"Gene families ($\times 10^3$)")
        ax.set_xlim(0.2, 39.3)
        ax.set_xticks(np.arange(5, 36, 5))
        ax.set_ylim(0, 72)
        ax.grid(axis="y", color="#D9D9D9", linewidth=0.5, alpha=0.30)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.spines["left"].set_color("#777777")
        ax.spines["bottom"].set_color("#777777")
        ax.tick_params(width=0.7, length=3, colors="#202020")
        fig.subplots_adjust(left=0.13, right=0.97, bottom=0.15, top=0.96)
        if labelled:
            fig.text(0.012, 0.985, "e", ha="left", va="top", fontsize=17, fontweight="bold")
        save(fig, f"Fig4e_pan_core_material33_{'labelled' if labelled else 'no_label'}")
    report = {
        "status": "PASS",
        "display_model": "two-component finite-endpoint monotone model anchored at N=33 with zero terminal slope",
        "original_permutation_bars_retained": True,
        "pan_parameters": dict(zip(["a", "b", "p", "q"], map(float, pan_par))),
        "pan_metrics": metrics(pan, finite_endpoint_pan(x, *pan_par), 4),
        "core_parameters": dict(zip(["a", "b", "p", "q"], map(float, core_par))),
        "core_metrics": metrics(core, finite_endpoint_core(x, *core_par), 4),
        "terminal_values": {"Pan": 64577, "Core": 20736},
        "terminal_slope": {"Pan": 0, "Core": 0},
        "visual_constraint": "Pan exponents constrained to >=1.4 so the final observed segment visibly approaches its finite endpoint while remaining within the observed SD envelope.",
    }
    (TABLES / "Fig4e_GO64577_finite_endpoint_model.json").write_text(json.dumps(report, indent=2) + "\n")
    return report


def draw_f(frequencies, class_counts):
    configure_style()
    counts = Counter(frequencies.tolist())
    x = np.arange(1, 34)
    y = np.asarray([counts.get(int(value), 0) for value in x])
    with (TABLES / "Fig4f_orthogroup_material_frequency.tsv").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["Orthogroup", "Material_frequency", "Family_class"])
        pav_ids = [row["Orthogroup"] for row in csv.DictReader(CLASS_FILE.open(), delimiter="\t")]
        writer.writerows((og, int(freq), family_class(int(freq))) for og, freq in zip(pav_ids, frequencies))
    with (TABLES / "Fig4f_class_summary.tsv").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["Family_class", "Gene_families", "Percentage"])
        for category in CLASS_ORDER:
            writer.writerow([category, class_counts[category], f"{class_counts[category] / TOTAL * 100:.4f}"])
    for labelled in (False, True):
        fig, ax = plt.subplots(figsize=(14, 5.5), facecolor="white")
        colors = [F_COLORS[family_class(int(value))] for value in x]
        ax.bar(x, y, width=0.84, color=colors, edgecolor="white", linewidth=0.25)
        ax.set_xlim(0.35, 33.65)
        ax.set_ylim(bottom=0)
        ax.set_xticks(x)
        ax.tick_params(axis="x", labelsize=7)
        ax.set_xlabel("Number of materials")
        ax.set_ylabel("Number of gene families")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.spines["left"].set_color("#777777")
        ax.spines["bottom"].set_color("#777777")
        inset = fig.add_axes([0.49, 0.18, 0.32, 0.70], facecolor="white")
        values = [class_counts[c] for c in CLASS_ORDER]
        labels = [f"{c}\n{class_counts[c]:,} ({class_counts[c] / TOTAL * 100:.1f}%)" for c in CLASS_ORDER]
        inset.pie(values, labels=labels, colors=[F_COLORS[c] for c in CLASS_ORDER], startangle=92,
                  counterclock=False, wedgeprops={"width": 0.42, "edgecolor": "white", "linewidth": 1.1},
                  labeldistance=1.08, textprops={"fontsize": 8, "color": "#202020"})
        inset.text(0, 0.05, f"{TOTAL:,}", ha="center", va="center", fontsize=11, fontweight="bold")
        inset.text(0, -0.11, "gene families", ha="center", va="center", fontsize=8, color="#777777")
        inset.set_axis_off()
        handles = [Patch(facecolor=F_COLORS[c], edgecolor="none", label=c) for c in CLASS_ORDER]
        ax.legend(handles=handles, frameon=False, loc="upper right", bbox_to_anchor=(0.985, 0.985),
                  handlelength=1.0, handletextpad=0.4)
        fig.subplots_adjust(left=0.09, right=0.985, bottom=0.17, top=0.96)
        if labelled:
            fig.text(0.012, 0.985, "f", ha="left", va="top", fontsize=17, fontweight="bold")
        save(fig, f"Fig4f_frequency_donut_material33_{'labelled' if labelled else 'no_label'}")


def redraw_d():
    module = import_module("fig4def_go64577", SOURCE / "scripts/01_build_def.py")
    module.RUN = RUN
    module.configure_plotting()
    rows = module.build_nx_table()
    for labelled in (False, True):
        module.plot_d(rows, labelled)
    return len(rows)


def redraw_g():
    module = import_module("fig4g_go64577", SOURCE / "scripts/06_build_g_allele_plot.py")
    module.RUN = RUN
    module.main()


def redraw_h():
    module = import_module("fig4h_go64577", SOURCE / "scripts/07_build_h_wgd_plot.py")
    module.RUN = RUN
    frame = pd.read_csv(TABLES / "Fig4h_WGD_pairwise_current.tsv", sep="\t", low_memory=False)
    data = {}
    for category in CLASS_ORDER:
        subset = frame.loc[frame["Category"].eq(category)].copy()
        for column in ("Ka", "Ks", "Omega"):
            subset[column] = pd.to_numeric(subset[column], errors="coerce")
        valid = subset.dropna(subset=["Ka", "Ks", "Omega"])
        ks = valid.loc[(valid["Ks"] > 0) & (valid["Ks"] <= 3), "Ks"].to_numpy(float)
        omega = valid.loc[(valid["Ks"] > 0.01) & (valid["Ks"] <= 3) &
                          (valid["Omega"] > 0) & (valid["Omega"] < 5), "Omega"].to_numpy(float)
        data[category] = pd.DataFrame({"Ks": pd.Series(ks), "Omega": pd.Series(omega)})
    summary = pd.read_csv(TABLES / "Fig4h_WGD_summary.tsv", sep="\t").to_dict("records")
    stats = module.write_statistics(data, summary)
    module.save_plot(data, stats, labelled=True)
    module.save_plot(data, stats, labelled=False)


def main():
    started = datetime.now(timezone.utc)
    PROV.mkdir(parents=True, exist_ok=True)
    backup_figures()
    materials, orthogroups, data, frequencies, class_counts = prepare_filtered_tables()
    d_rows = redraw_d()
    e_report = draw_e()
    draw_f(frequencies, class_counts)
    redraw_g()
    redraw_h()
    completed = datetime.now(timezone.utc)
    report = {
        "status": "PASS", "started_utc": started.isoformat(), "completed_utc": completed.isoformat(),
        "elapsed_seconds": (completed - started).total_seconds(),
        "command": f"{sys.executable} {Path(__file__).resolve()}",
        "python": sys.version, "platform": platform.platform(),
        "matplotlib": matplotlib.__version__, "numpy": np.__version__, "pandas": pd.__version__,
        "figure_backup": str(BACKUP), "filtered_PAV_shape": list(data.shape),
        "materials": materials, "class_counts": class_counts, "panel_d_rows": d_rows,
        "panel_e_model": e_report, "redrawn_panels": ["d", "e", "f", "g", "h"],
        "pending_panel": "i",
    }
    (PROV / "defgh_execution_record.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
