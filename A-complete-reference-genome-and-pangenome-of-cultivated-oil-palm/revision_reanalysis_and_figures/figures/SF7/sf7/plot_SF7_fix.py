#!/usr/bin/env python3
"""Supplementary Fig. 7 redraw: FAD2 qPCR shown per primer assay.

Two input versions (select with --version):
  dct  (default) Source Data sheet qPCR_stage_summary: log2(FAD2/ACTIN) = -dCt, recomputed from
       raw Ct wells (mean Ct target - mean Ct ACTIN per stage x material x assay). Traceable to the
       instrument Ct grid; this is the version currently in the submission Source Data.
  lab  F2_SUPP_14_FAD2_qPCR_assay_level_tidy.tsv (2026-08-06): lab-supplied relative expression
       values from 数据整理.xlsx (hard-coded, calibrator not documented). This is the version the
       current Supplementary Fig. 7 image and legend q values were drawn from.
Panel a: one row per primer assay (technical measurements). Panel b: stage-level mean across the
three assays, FL versus each comparator, two-sided exact Wilcoxon signed-rank test paired by stage,
Benjamini-Hochberg within phase. q values are printed; no significance stars.
"""
import argparse
import pickle
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import wilcoxon

HERE = Path(__file__).resolve().parent
MM = 1 / 25.4
TEAL, RED, DARK = "#2A9D8F", "#D95F5F", "#242A30"
MAT = {
    "FL": dict(color=TEAL, marker="o", ls="-"),
    "TN": dict(color=RED, marker="o", ls="-"),
    "NS": dict(color="#8C8C8C", marker="s", ls="--"),
    "TK": dict(color="#4A4A4A", marker="^", ls=":"),
}
ORDER = ["FL", "TN", "NS", "TK"]

mpl.rcParams.update({
    "font.family": "Arial", "font.size": 7, "axes.labelsize": 7, "axes.titlesize": 7,
    "xtick.labelsize": 6, "ytick.labelsize": 6, "legend.fontsize": 6, "axes.linewidth": 0.6,
    "xtick.major.width": 0.6, "ytick.major.width": 0.6, "xtick.major.size": 2.5, "ytick.major.size": 2.5,
    "axes.unicode_minus": True, "mathtext.fontset": "custom", "mathtext.rm": "Arial",
    "mathtext.it": "Arial:italic", "mathtext.bf": "Arial:bold", "pdf.fonttype": 42, "svg.fonttype": "none",
})


def clean(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_axisbelow(True)


def bh(p):
    p = np.asarray(p, float)
    n = len(p)
    order = np.argsort(p)
    q = np.empty(n)
    prev = 1.0
    for rank in range(n, 0, -1):
        i = order[rank - 1]
        prev = min(prev, p[i] * n / rank)
        q[i] = prev
    return q


def load(version):
    """Return long table: stage_index(1-19), stage, phase, material, assay(1-3), value."""
    if version == "dct":
        # FIX 2026-09-24: corrected stage summary (50 d primer 3 F7/G7 swap restored to instrument CSV)
        s = pd.read_csv(HERE / "SD_SF7_stage_summary.tsv", sep="\t")
        rows = []
        for _, r in s.iterrows():
            for a in (1, 2, 3):
                rows.append(dict(stage_index=int(r["Stage index"]), stage=r["Stage"],
                                 phase="Development" if r["Phase"] == "Development" else "Postharvest",
                                 material=r["Material"], assay=a, value=r[f"Primer {a} log₂(target/ACTIN)"]))
        ylab = r"log$_2$(FAD2/ACTIN) ($-\Delta$Ct)"
        return pd.DataFrame(rows), ylab
    t = pd.read_csv(HERE / "F2_SUPP_14_FAD2_qPCR_assay_level_tidy.tsv", sep="\t")
    t = t.assign(stage_index=t.stage_index + 1, stage=t.timepoint,
                 phase=np.where(t.phase == "fruit_development", "Development", "Postharvest"),
                 material=t["sample"], assay=t.assay.str.extract(r"(\d)").astype(int)[0],
                 value=t.log2_relative_expression)
    return t[["stage_index", "stage", "phase", "material", "assay", "value"]], r"log$_2$ relative FAD2 RNA"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--version", choices=["dct", "lab"], default="dct")
    args = ap.parse_args()
    long, ylab = load(args.version)
    tag = "" if args.version == "dct" else "_ALT_lab_relative_expression"
    stages = long[["stage_index", "stage"]].drop_duplicates().sort_values("stage_index")
    stage_lab = [s.replace("d", " d").replace("h", " h") for s in stages.stage]

    # stage-level mean across assays + tests
    mean = long.groupby(["stage_index", "stage", "phase", "material"], as_index=False).value.mean()
    tests = []
    for ph in ("Development", "Postharvest"):
        piv = mean[mean.phase == ph].pivot(index="stage_index", columns="material", values="value")
        ps, diffs = [], []
        for m in ("TN", "NS", "TK"):
            ps.append(wilcoxon(piv.FL, piv[m], alternative="two-sided", method="exact").pvalue)
            diffs.append(np.median(piv.FL - piv[m]))
        for m, p, q, dm in zip(("TN", "NS", "TK"), ps, bh(ps), diffs):
            tests.append(dict(phase=ph, contrast=f"FL - {m}", n_paired_stages=len(piv),
                              median_paired_difference=dm, P=p, q_BH_within_phase=q))
    tests = pd.DataFrame(tests)

    fig = plt.figure(figsize=(180 * MM, 175 * MM))
    gsa = fig.add_gridspec(3, 1, left=0.085, right=0.86, top=0.965, bottom=0.44, hspace=0.35)
    gsb = fig.add_gridspec(1, 2, left=0.085, right=0.86, top=0.33, bottom=0.055, wspace=0.12)
    lo, hi = long.value.min(), long.value.max()
    pad = 0.08 * (hi - lo)
    axes_a = []
    for i, a in enumerate((1, 2, 3)):
        ax = fig.add_subplot(gsa[i, 0])
        sub = long[long.assay == a]
        for m in ORDER:
            s = sub[sub.material == m].sort_values("stage_index")
            st = MAT[m]
            ax.plot(s.stage_index, s.value, color=st["color"], ls=st["ls"], lw=1.0, marker=st["marker"],
                    ms=2.6, label=m, zorder=3 if m in ("FL", "TN") else 2)
        ax.axvline(13.5, color="#B9C0C7", lw=0.6, ls="--", zorder=0)
        ax.set_xlim(0.5, 19.5)
        ax.set_ylim(lo - pad, hi + pad)
        ax.set_xticks(stages.stage_index, stage_lab if i == 2 else [""] * len(stage_lab),
                      rotation=45 if i == 2 else 0, ha="right" if i == 2 else "center")
        ax.text(0.005, 0.97, f"Primer assay {a}", transform=ax.transAxes, ha="left", va="top", fontsize=6.5,
                fontweight="bold", color=DARK)
        if i == 0:
            ax.text(7, 1.02, "Development", transform=ax.get_xaxis_transform(), ha="center", va="bottom",
                    fontsize=6, color="#555555")
            ax.text(16.5, 1.02, "Postharvest", transform=ax.get_xaxis_transform(), ha="center",
                    va="bottom", fontsize=6, color="#555555")
            ax.legend(frameon=False, loc="upper left", bbox_to_anchor=(1.01, 1.0), handlelength=2.2)
        if i == 1:
            ax.set_ylabel(ylab)
        if i == 2:
            ax.set_xlabel("Sampling stage")
        ax.grid(axis="y", color="#EEEEEE", lw=0.5)
        clean(ax)
        axes_a.append(ax)

    rng = np.random.default_rng(7)
    axes_b = []
    for j, ph in enumerate(("Development", "Postharvest")):
        ax = fig.add_subplot(gsb[0, j], sharey=axes_b[0] if axes_b else None)
        sub = mean[mean.phase == ph]
        vals = [sub[sub.material == m].value.values for m in ORDER]
        bp = ax.boxplot(vals, positions=range(1, 5), widths=0.55, showfliers=False, patch_artist=True,
                        medianprops=dict(color=DARK, lw=1.0), whiskerprops=dict(lw=0.6),
                        capprops=dict(lw=0.6), boxprops=dict(lw=0.6))
        for patch, m in zip(bp["boxes"], ORDER):
            patch.set_facecolor(MAT[m]["color"])
            patch.set_alpha(0.25)
        for k, (m, v) in enumerate(zip(ORDER, vals), start=1):
            ax.scatter(k + rng.uniform(-0.12, 0.12, len(v)), v, s=7, color=MAT[m]["color"], lw=0, zorder=3)
        ax.axhline(0, color="#B9C0C7", lw=0.5, zorder=0)
        ax.set_xticks(range(1, 5), ORDER)
        ax.set_title(f"{ph} ({sub.stage_index.nunique()} stages)", pad=3)
        clean(ax)
        axes_b.append(ax)
    axes_b[0].set_ylabel("Stage-level mean across assays\n" + ylab)
    plt.setp(axes_b[1].get_yticklabels(), visible=False)
    ymax = mean.value.max()
    ymin = mean.value.min()
    step = 0.11 * (ymax - ymin)
    for ax, ph in zip(axes_b, ("Development", "Postharvest")):
        t = tests[tests.phase == ph].reset_index(drop=True)
        for k, (_, r) in enumerate(t.iterrows()):
            x2 = 2 + k
            y = ymax + step * (0.6 + k)
            ax.plot([1, 1, x2, x2], [y - step * 0.18, y, y, y - step * 0.18], color="#555555", lw=0.6)
            ax.text((1 + x2) / 2, y + step * 0.05, f"q = {r.q_BH_within_phase:.3g}", ha="center", va="bottom",
                    fontsize=5.5)
        ax.set_ylim(ymin - step * 0.5, ymax + step * 3.9)

    for ax, s in [(axes_a[0], "a"), (axes_b[0], "b")]:
        bb = ax.get_position()
        fig.text(bb.x0 - 0.07, bb.y1 + 0.012, s, fontsize=9, fontweight="bold", va="bottom")

    out = HERE / "out"
    out.mkdir(exist_ok=True)
    for ext in ("pdf", "png"):
        fig.savefig(out / f"Supplementary_Fig_07{tag}.{ext}", dpi=600, facecolor="white")
    long.sort_values(["assay", "stage_index", "material"]).to_csv(
        out / f"SF7a_assay_level_values{tag}.tsv", sep="\t", index=False)
    mean.to_csv(out / f"SF7b_stage_level_means{tag}.tsv", sep="\t", index=False)
    tests.to_csv(out / f"SF7b_paired_stage_tests{tag}.tsv", sep="\t", index=False)
    print(tests.round(4).to_string())


if __name__ == "__main__":
    main()
