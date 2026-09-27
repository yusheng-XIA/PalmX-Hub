#!/usr/bin/env python3
"""Build two publication-oriented FL/TN ASE figures panel-by-panel.

The layouts deliberately follow the two user-supplied examples while all
statistics come from the accepted graph-based FL/TN ASE workflow.  FL never
receives parental cis/trans labels; both varieties use replicated ASE calls.
"""

from __future__ import annotations

import bisect
import gzip
import math
import re
from pathlib import Path

import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize, TwoSlopeNorm
from matplotlib.gridspec import GridSpecFromSubplotSpec
from matplotlib.lines import Line2D
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, PathPatch, Rectangle
from matplotlib.path import Path as MplPath
from matplotlib.transforms import Bbox
import numpy as np
import pandas as pd
from scipy import stats


OUT = Path(__file__).resolve().parent
ANALYSIS = Path("${ANALYSIS_DIR}")
FIG3 = ANALYSIS / "22_answer_reviews/00_ms/03_V3/03_figure3"
ASE_RUN = FIG3 / "01_ASE/00_shared/runs/RUN-ASE-HETEROSIS-DOWNSTREAM-V2-001/output"
ASE_ROOT = FIG3 / "01_ASE"
CURRENT = ANALYSIS / "22_answer_reviews/00_ms/03_V3/02_figure/_runs/RUN-ASE-FIGURES-CURRENT-V3-001/output"
EXACT = ANALYSIS / "22_answer_reviews/00_ms/03_V3/02_figure/_runs/RUN-ASE-FIGURES-EXACT-EXAMPLE-V3-002/output"
TE_ROOT = ANALYSIS / "22_answer_reviews/03_TE_re"

CLASS_ORDER = ["HapDom", "Sub", "NoDiff", "NoASE"]
CLASS_LABEL = {
    "HapDom": "Stable bias",
    "Sub": "Bias switching",
    "NoDiff": "Stage-limited",
    "NoASE": "No robust ASE",
}
CLASS_COLOR = {
    "HapDom": "#E76F51",
    "Sub": "#5AB4AC",
    "NoDiff": "#F2C14E",
    "NoASE": "#B9BDC2",
}
TN_A, TN_B = "#3B82B4", "#E07A5F"
FL_A, FL_B = "#009E73", "#CC79A7"
STAGES = ["0d", "15d", "35d", "50d", "65d", "80d", "95d", "110d", "125d",
          "140d", "155d", "170d", "185d", "12h", "24h", "36h", "48h", "60h", "72h"]
PHASES = ["Early", "Mid", "Late", "Postharvest"]
PHASE_LABELS = ["Early\n0–50 d", "Mid\n65–110 d", "Late\n125–185 d", "Postharvest\n12–72 h"]


def set_style() -> None:
    mpl.rcParams.update({
        "font.family": "DejaVu Sans",
        "font.size": 7.2,
        "axes.titlesize": 8.6,
        "axes.labelsize": 8,
        "xtick.labelsize": 6.6,
        "ytick.labelsize": 6.6,
        "legend.fontsize": 6.5,
        "axes.linewidth": 0.7,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "svg.fonttype": "none",
        "savefig.facecolor": "white",
    })


def save_figure(fig: plt.Figure, stem: str) -> None:
    for ext in ("pdf", "svg", "png"):
        kw = {"bbox_inches": "tight", "facecolor": "white"}
        if ext == "png":
            kw["dpi"] = 600
        fig.savefig(OUT / f"{stem}.{ext}", **kw)
    plt.close(fig)


def save_panel_groups(fig: plt.Figure, groups: dict[str, list[plt.Axes]]) -> None:
    """Export selected axes as standalone, editable panels without re-rendering."""
    for stem, axes in groups.items():
        keep = set(axes)
        visibility = {axis: axis.get_visible() for axis in fig.axes}
        for axis in fig.axes:
            axis.set_visible(axis in keep)
        fig.canvas.draw()
        renderer = fig.canvas.get_renderer()
        boxes = [axis.get_tightbbox(renderer) for axis in axes if axis.get_visible()]
        if not boxes:
            for axis, visible in visibility.items():
                axis.set_visible(visible)
            continue
        bbox = Bbox.union(boxes).transformed(fig.dpi_scale_trans.inverted()).padded(0.10)
        for ext in ("pdf", "svg", "png"):
            kw = {"bbox_inches": bbox, "facecolor": "white"}
            if ext == "png":
                kw["dpi"] = 600
            fig.savefig(OUT / f"{stem}.{ext}", **kw)
        for axis, visible in visibility.items():
            axis.set_visible(visible)
    fig.canvas.draw()


def panel(ax: plt.Axes, letter: str, x: float = -0.10, y: float = 1.06) -> None:
    ax.text(x, y, letter, transform=ax.transAxes, fontsize=14, fontweight="bold",
            ha="left", va="bottom")


def clean(ax: plt.Axes, grid: str | None = None) -> None:
    ax.spines[["top", "right"]].set_visible(False)
    if grid:
        ax.grid(axis=grid, color="#D9DEE2", lw=0.45, alpha=0.75)


def read_gene_gff(path: Path) -> pd.DataFrame:
    rows = []
    with path.open(errors="replace") as handle:
        for line in handle:
            if not line or line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "gene":
                continue
            match = re.search(r"(?:^|;)ID=([^;]+)", f[8])
            if not match:
                continue
            chrom = f[0].split("__")[-1]
            rows.append((match.group(1), chrom, int(f[3]), int(f[4]), f[6]))
    return pd.DataFrame(rows, columns=["gene_id", "chrom", "start", "end", "strand"])


def read_ltr_gff(path: Path) -> dict[str, tuple[list[int], list[int]]]:
    raw: dict[str, list[tuple[int, int]]] = {}
    with path.open(errors="replace") as handle:
        for line in handle:
            if not line or line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 5:
                continue
            raw.setdefault(f[0].split("__")[-1], []).append((int(f[3]), int(f[4])))
    out = {}
    for chrom, values in raw.items():
        values.sort()
        starts = [x[0] for x in values]
        prefix_max = []
        current = -1
        for _, end in values:
            current = max(current, end)
            prefix_max.append(current)
        out[chrom] = (starts, prefix_max)
    return out


def any_overlap(index: dict[str, tuple[list[int], list[int]]], chrom: str,
                lo: int, hi: int) -> bool:
    if chrom not in index or hi < 1:
        return False
    lo = max(1, lo)
    starts, prefix_max = index[chrom]
    i = bisect.bisect_right(starts, hi) - 1
    return i >= 0 and prefix_max[i] >= lo


def oriented_segment(row, region: str, i: int, n: int) -> tuple[int, int]:
    start, end, strand = int(row.start), int(row.end), row.strand
    if region == "body":
        length = max(1, end - start + 1)
        a = start + math.floor(i * length / n)
        b = start + math.ceil((i + 1) * length / n) - 1
        if strand == "-":
            a2 = end - math.ceil((i + 1) * length / n) + 1
            b2 = end - math.floor(i * length / n)
            a, b = a2, b2
        return min(a, b), max(a, b)
    width = 5000 / n
    a0 = math.floor(i * width)
    b0 = math.ceil((i + 1) * width) - 1
    if region == "upstream":
        if strand == "+":
            return start - 5000 + a0, start - 5000 + b0
        return end + 1 + (5000 - 1 - b0), end + 1 + (5000 - 1 - a0)
    if strand == "+":
        return end + 1 + a0, end + 1 + b0
    return start - 5000 + (5000 - 1 - b0), start - 5000 + (5000 - 1 - a0)


def build_te_profiles(classes: pd.DataFrame, bridge: pd.DataFrame) -> pd.DataFrame:
    configs = [
        ("FL", "Africa hap2", "gene_africa",
         ASE_ROOT / "00_shared/runs/RUN-FL-VG-RNAIDX-V2-001/attempt_1/annotation/fl_ase.chr_only.gff3",
         TE_ROOT / "Africa_hap2/intact_LTR.gff3"),
        ("FL", "American hap1", "gene_american",
         ASE_ROOT / "02_seedless/ase/annotation/American_hap1.graph.prepared.gff3",
         TE_ROOT / "American_hap1/intact_LTR.gff3"),
        ("TN", "Dura-like", "gene_dura",
         ASE_ROOT / "00_shared/runs/RUN-TN-VG-RNAIDX-V2-001/attempt_2/annotation/tn_ase.chr_only.gff3",
         TE_ROOT / "dura/intact_LTR.gff3"),
        ("TN", "Pisifera-like", "gene_pisifera",
         ASE_ROOT / "01_bk/ase/annotation/EG_pisifera.graph.prepared.gff3",
         TE_ROOT / "pisifera/intact_LTR.gff3"),
    ]
    result = []
    for analysis, haplotype, gene_col, gff, te_gff in configs:
        base = classes[classes.analysis == analysis][["gene_id", "overall_class"]].copy()
        if gene_col in ("gene_africa", "gene_dura"):
            mapping = base.rename(columns={"gene_id": gene_col})
        else:
            left = "gene_africa" if analysis == "FL" else "gene_dura"
            mapping = base.merge(bridge[[left, gene_col]], left_on="gene_id", right_on=left,
                                 how="inner")[[gene_col, "overall_class"]]
        genes = read_gene_gff(gff).merge(mapping, left_on="gene_id", right_on=gene_col,
                                         how="inner")
        te_index = read_ltr_gff(te_gff)
        specs = [("Upstream 5 kb", "upstream", 36), ("Gene body", "body", 72),
                 ("Downstream 5 kb", "downstream", 36)]
        offset = 0
        for region_label, region_key, n_bins in specs:
            for cls in CLASS_ORDER:
                subset = genes[genes.overall_class == cls]
                if subset.empty:
                    continue
                hit = np.zeros(n_bins, dtype=int)
                for row in subset.itertuples(index=False):
                    for i in range(n_bins):
                        lo, hi = oriented_segment(row, region_key, i, n_bins)
                        hit[i] += int(any_overlap(te_index, row.chrom, lo, hi))
                for i, count in enumerate(hit):
                    result.append({
                        "analysis": analysis, "haplotype": haplotype,
                        "overall_class": cls, "region": region_label,
                        "bin_within_region": i + 1, "plot_bin": offset + i,
                        "genes": len(subset), "genes_overlapping_intact_LTR": int(count),
                        "occupancy_percent": 100 * count / len(subset),
                    })
            offset += n_bins
    result = pd.DataFrame(result)
    result.to_csv(OUT / "source_C_intact_LTR_profiles.tsv", sep="\t", index=False)
    return result


def flow_patch(ax, x0, x1, y0a, y0b, y1a, y1b, color, alpha=0.38):
    c = 0.42 * (x1 - x0)
    verts = [(x0, y0a), (x0 + c, y0a), (x1 - c, y1a), (x1, y1a),
             (x1, y1b), (x1 - c, y1b), (x0 + c, y0b), (x0, y0b), (x0, y0a)]
    codes = [MplPath.MOVETO, MplPath.CURVE4, MplPath.CURVE4, MplPath.CURVE4,
             MplPath.LINETO, MplPath.CURVE4, MplPath.CURVE4, MplPath.CURVE4,
             MplPath.CLOSEPOLY]
    ax.add_patch(PathPatch(MplPath(verts, codes), facecolor=color,
                           edgecolor="none", alpha=alpha, zorder=1))


def draw_two_column_alluvial(ax: plt.Axes, counts: pd.DataFrame) -> None:
    total = counts.orthogroups.sum()
    gap, barw = 0.012, 0.10
    bounds = {}
    for side, col in [(0, "overall_class_FL"), (1, "overall_class_TN")]:
        sums = counts.groupby(col).orthogroups.sum()
        usable = 1 - gap * (len(CLASS_ORDER) - 1)
        y = 1.0
        for cat in CLASS_ORDER:
            h = usable * sums.get(cat, 0) / total
            bounds[(side, cat)] = (y - h, y)
            y -= h + gap
    lc = {c: bounds[(0, c)][0] for c in CLASS_ORDER}
    rc = {c: bounds[(1, c)][0] for c in CLASS_ORDER}
    for source in CLASS_ORDER:
        for target in CLASS_ORDER:
            q = counts[(counts.overall_class_FL == source) &
                       (counts.overall_class_TN == target)]
            k = int(q.orthogroups.iloc[0]) if len(q) else 0
            if not k:
                continue
            h = (1 - gap * 3) * k / total
            flow_patch(ax, barw / 2, 1 - barw / 2,
                       lc[source], lc[source] + h, rc[target], rc[target] + h,
                       CLASS_COLOR[source])
            lc[source] += h
            rc[target] += h
    for side, title in [(0, "FL"), (1, "TN")]:
        for cat in CLASS_ORDER:
            y0, y1 = bounds[(side, cat)]
            ax.add_patch(Rectangle((side - barw / 2, y0), barw, y1-y0,
                                   facecolor=CLASS_COLOR[cat], edgecolor="white",
                                   lw=0.7, zorder=4))
            if y1 - y0 > 0.055:
                ax.text(side + (0.075 if side == 0 else -0.075), (y0+y1)/2,
                        f"{100*(y1-y0)/(1-gap*3):.1f}%", va="center",
                        ha="left" if side == 0 else "right", fontsize=5.7)
        ax.text(side, 1.035, title, ha="center", va="bottom", fontweight="bold")
    ax.set_xlim(-0.30, 1.30); ax.set_ylim(-0.01, 1.08); ax.axis("off")


def kde(ax: plt.Axes, values, color, label, xmax: float, ls: str = "-") -> None:
    v = pd.to_numeric(pd.Series(values), errors="coerce").to_numpy(float)
    v = v[np.isfinite(v) & (v >= 0) & (v <= xmax)]
    if len(v) < 20 or np.unique(v).size < 3:
        return
    x = np.linspace(0, xmax, 320)
    y = stats.gaussian_kde(v)(x)
    ax.plot(x, y, color=color, lw=1.15, ls=ls, label=label)


def build_candidate_source(gene_stage: pd.DataFrame) -> pd.DataFrame:
    candidates = pd.DataFrame([
        ("FL", "evm.TU.chr09B.1430", "MADS34"),
        ("FL", "evm.TU.chr08B.1689", "ARF2"),
        ("TN", "evm.TU.chr10.635", "DGAT"),
        ("TN", "evm.TU.chr07.37", "ACCase"),
    ], columns=["analysis", "gene_id", "display_label"])
    out = gene_stage.merge(candidates, on=["analysis", "gene_id"], how="inner")
    out["allele_A_fraction"] = out.ref_fraction
    out["allele_B_fraction"] = 1 - out.ref_fraction
    out.to_csv(OUT / "source_H_candidate_trajectories.tsv", sep="\t", index=False)
    return out


def build_figure_one(mirror: pd.DataFrame, classes: pd.DataFrame,
                     te: pd.DataFrame, overlap_counts: pd.DataFrame,
                     kaks: pd.DataFrame, snp: pd.DataFrame,
                     candidate: pd.DataFrame) -> None:
    fig = plt.figure(figsize=(15.7, 11.2))
    gs = fig.add_gridspec(3, 16, height_ratios=[1.02, 1.05, 1.0],
                          hspace=0.40, wspace=0.80)
    panel_groups: dict[str, list[plt.Axes]] = {}

    # A: allele-pair framework.
    ax = fig.add_subplot(gs[0, 0:5]); ax.axis("off"); panel(ax, "A", -0.05, 1.02)
    panel_groups["01A_allele_framework"] = [ax]
    ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    ax.text(0.02, 0.95, "Graph-based allele framework", fontweight="bold", fontsize=9)
    rows = [(0.66, "FL", "Africa hap2", "American hap1", FL_A, FL_B, "17,616 eligible genes"),
            (0.25, "TN", "Dura/TK-like", "Pisifera/NS-like", TN_A, TN_B, "11,409 eligible genes")]
    for y, name, a, b, ca, cb, nlab in rows:
        ax.text(0.02, y+0.11, name, fontsize=10, fontweight="bold")
        for x, label, color in [(0.17, a, ca), (0.65, b, cb)]:
            box = FancyBboxPatch((x, y), 0.25, 0.17, boxstyle="round,pad=0.012,rounding_size=0.025",
                                 facecolor=color, edgecolor="#30343B", lw=0.8, alpha=0.90)
            ax.add_patch(box); ax.text(x+0.125, y+0.085, label, color="white",
                                      ha="center", va="center", fontweight="bold", fontsize=7.5)
        ax.add_patch(FancyArrowPatch((0.43, y+0.085), (0.64, y+0.085), arrowstyle="<->",
                                     mutation_scale=10, lw=1.0, color="#4C5661"))
        ax.text(0.535, y+0.145, "1:1 allele pair", ha="center", fontsize=6.2)
        ax.text(0.535, y-0.045, nlab, ha="center", color="#5F6770", fontsize=6.2)
    ax.text(0.02, 0.02, "Replicate-aware beta-binomial ASE; BH FDR < 0.05; |log2(A/B)| ≥ 0.5",
            fontsize=6.1, color="#555D66")

    # B: overall ASE-class proportions.
    ax = fig.add_subplot(gs[0, 5:8]); panel(ax, "B", -0.16, 1.02)
    panel_groups["01B_ASE_class_proportions"] = [ax]
    summary = (classes.groupby(["analysis", "overall_class"]).size()
               .rename("genes").reset_index())
    summary["percentage"] = summary.groupby("analysis").genes.transform(lambda x: 100*x/x.sum())
    summary.to_csv(OUT / "source_B_ASE_class_proportions.tsv", sep="\t", index=False)
    left = np.zeros(2)
    for cat in CLASS_ORDER:
        vals = [float(summary[(summary.analysis == a) & (summary.overall_class == cat)].percentage.iloc[0])
                for a in ["FL", "TN"]]
        ax.barh([1, 0], vals, left=left, color=CLASS_COLOR[cat], height=0.48,
                edgecolor="white", lw=0.5, label=CLASS_LABEL[cat])
        for yi, l, v in zip([1, 0], left, vals):
            if v > 8: ax.text(l+v/2, yi, f"{v:.1f}%", ha="center", va="center", fontsize=5.8)
        left += vals
    ax.set_yticks([1, 0], ["FL", "TN"]); ax.set_xlim(0, 100); ax.set_xlabel("Eligible genes (%)")
    ax.set_title("Temporal ASE classes", loc="left", fontweight="bold")
    ax.legend(frameon=False, loc="lower center", bbox_to_anchor=(0.5, -0.58), ncol=2)
    clean(ax, "x")

    # C: four-haplotype intact-LTR profiles.
    cgrid = GridSpecFromSubplotSpec(4, 1, subplot_spec=gs[0, 8:16], hspace=0.10)
    hap_order = [("FL", "Africa hap2"), ("FL", "American hap1"),
                 ("TN", "Dura-like"), ("TN", "Pisifera-like")]
    caxes = []
    for i, (analysis, hap) in enumerate(hap_order):
        cax = fig.add_subplot(cgrid[i, 0]); caxes.append(cax)
        part = te[(te.analysis == analysis) & (te.haplotype == hap)]
        for cat in CLASS_ORDER:
            z = part[part.overall_class == cat].sort_values("plot_bin")
            cax.plot(z.plot_bin, z.occupancy_percent, color=CLASS_COLOR[cat], lw=1.0,
                     ls="--" if cat == "NoASE" else "-", label=CLASS_LABEL[cat])
        cax.axvline(35.5, color="#444", lw=0.6, ls="--"); cax.axvline(107.5, color="#444", lw=0.6, ls="--")
        cax.text(0.01, 0.82, f"{analysis} · {hap}", transform=cax.transAxes, fontsize=6.2,
                 bbox=dict(facecolor="white", edgecolor="none", alpha=0.75, pad=1))
        cax.set_xlim(0, 143); cax.set_ylim(bottom=0); cax.grid(axis="y", color="#E0E4E7", lw=0.35)
        cax.spines[["top", "right"]].set_visible(False)
        if i < 3: cax.set_xticklabels([])
        else:
            cax.set_xticks([0, 35.5, 71.5, 107.5, 143], ["−5 kb", "TSS", "Gene body", "TES", "+5 kb"])
    caxes[0].set_title("Intact-LTR occupancy around ASE genes", loc="left", fontweight="bold")
    caxes[1].set_ylabel("Genes overlapping intact LTR (%)")
    caxes[0].legend(frameon=False, ncol=4, loc="upper right", bbox_to_anchor=(1.0, 1.80))
    panel(caxes[0], "C", -0.10, 1.10)
    panel_groups["01C_intact_LTR_profiles"] = caxes

    # D: stage-resolved mirror bars.
    dgrid = GridSpecFromSubplotSpec(2, 1, subplot_spec=gs[1, 0:8], hspace=0.15)
    daxes = []
    for i, (analysis, ca, cb, la, lb) in enumerate([
        ("TN", TN_A, TN_B, "Dura-like", "Pisifera-like"),
        ("FL", FL_A, FL_B, "Africa hap2", "American hap1")]):
        dax = fig.add_subplot(dgrid[i, 0]); daxes.append(dax)
        z = mirror[(mirror.analysis == analysis) & mirror.ase_call.isin(["Allele_A_biased", "Allele_B_biased"])]
        p = z.pivot(index="stage", columns="ase_call", values="percentage").reindex(STAGES)
        x = np.arange(19)
        dax.bar(x, p.Allele_A_biased, color=ca, width=0.75, label=la, edgecolor="white", lw=0.3)
        dax.bar(x, -p.Allele_B_biased, color=cb, width=0.75, label=lb, edgecolor="white", lw=0.3)
        dax.axhline(0, color="#303030", lw=0.7); dax.axvspan(12.5, 18.5, color="#FFF5E8", zorder=-3)
        dax.set_ylim(-36, 36); dax.set_yticks([-30,-15,0,15,30], ["30","15","0","15","30"])
        dax.text(0.01, 0.88, analysis, transform=dax.transAxes, fontweight="bold")
        dax.legend(frameon=False, ncol=2, loc="upper right")
        clean(dax, "y")
        if i == 0: dax.set_xticklabels([])
        else: dax.set_xticks(x, STAGES, rotation=45, ha="right")
    daxes[0].set_title("Stage-resolved directional ASE", loc="left", fontweight="bold")
    daxes[0].set_ylabel("Robust ASE genes (%)"); daxes[1].set_ylabel("Robust ASE genes (%)")
    panel(daxes[0], "D", -0.08, 1.08)
    panel_groups["01D_stage_mirror_ASE"] = daxes
    mirror.to_csv(OUT / "source_D_stage_mirror.tsv", sep="\t", index=False)

    # E: cross-variety alluvial.
    ax = fig.add_subplot(gs[1, 8:16]); panel(ax, "E", -0.08, 1.04)
    panel_groups["01E_FL_TN_alluvial"] = [ax]
    draw_two_column_alluvial(ax, overlap_counts)
    ax.set_title("Ortholog ASE-class conservation and rewiring", loc="left", fontweight="bold")
    handles = [Rectangle((0,0),1,1,color=CLASS_COLOR[c]) for c in CLASS_ORDER]
    ax.legend(handles, [CLASS_LABEL[c] for c in CLASS_ORDER], frameon=False, ncol=4,
              loc="lower center", bbox_to_anchor=(0.5, -0.12))
    ax.text(0.5, -0.03, f"Reviewed 1:1 orthogroups, n = {overlap_counts.orthogroups.sum():,}",
            transform=ax.transAxes, ha="center", fontsize=6.2, color="#555D66")
    overlap_counts.to_csv(OUT / "source_E_FL_TN_alluvial.tsv", sep="\t", index=False)

    # F: Ka/Ks density with Ks inset.
    ax = fig.add_subplot(gs[2, 0:5]); panel(ax, "F", -0.12, 1.05)
    panel_groups["01F_KaKs"] = [ax]
    for analysis, ls in [("FL", "-"), ("TN", "--")]:
        for cat in ["HapDom", "NoDiff", "Sub"]:
            z = kaks[(kaks.analysis == analysis) & (kaks.ASE_type == cat)]
            kde(ax, z.Ka_Ks, CLASS_COLOR[cat], f"{cat} · {analysis}", 3.0, ls)
    ax.axvline(1, color="#333", lw=0.7, ls=":"); ax.set_xlim(0,3)
    ax.set_xlabel("Ka/Ks ratio"); ax.set_ylabel("Density"); clean(ax)
    ax.set_title("Coding-sequence evolution", loc="left", fontweight="bold")
    inset = ax.inset_axes([0.54, 0.48, 0.43, 0.46])
    for analysis, ls in [("FL", "-"), ("TN", "--")]:
        for cat in ["HapDom", "NoDiff", "Sub"]:
            z = kaks[(kaks.analysis == analysis) & (kaks.ASE_type == cat)]
            kde(inset, z.Ks, CLASS_COLOR[cat], "", 0.10, ls)
    inset.set_xlim(0,0.10); inset.set_xlabel("Ks", fontsize=6); inset.set_ylabel("Density", fontsize=6)
    inset.tick_params(labelsize=5.5); inset.spines[["top","right"]].set_visible(False)
    legend_handles = [Line2D([0],[0],color=CLASS_COLOR[c],lw=1.4,label=CLASS_LABEL[c])
                      for c in ["HapDom","NoDiff","Sub"]]
    legend_handles += [Line2D([0],[0],color="#30343B",lw=1.2,ls="-",label="FL"),
                       Line2D([0],[0],color="#30343B",lw=1.2,ls="--",label="TN")]
    ax.legend(handles=legend_handles, frameon=False, ncol=2, fontsize=5.2,
              loc="lower right", bbox_to_anchor=(1.0, 0.02))
    kaks.to_csv(OUT / "source_F_KaKs_by_ASE_class.tsv", sep="\t", index=False)

    # G: SNP density around genes, both varieties.
    ggrid = GridSpecFromSubplotSpec(1, 2, subplot_spec=gs[2, 5:10], wspace=0.20)
    gaxes = []
    snp_summary_rows = []
    for j, analysis in enumerate(["FL", "TN"]):
        gax = fig.add_subplot(ggrid[0, j]); gaxes.append(gax)
        part = snp[snp.analysis == analysis]
        data, pos, colors = [], [], []
        x = 0
        for region_col, region in [("upstream2kb_per_kb", "Up"), ("gene_body_per_kb", "Gene"),
                                   ("downstream2kb_per_kb", "Down")]:
            groups = []
            for cat in ["HapDom", "NoDiff", "Sub"]:
                v = np.log10(1 + part.loc[part.overall_class == cat, region_col].dropna().to_numpy())
                data.append(v); pos.append(x); colors.append(CLASS_COLOR[cat]); groups.append(v); x += 1
                q = np.quantile(v, [0.25,0.5,0.75])
                snp_summary_rows.append({"analysis":analysis,"region":region,"overall_class":cat,
                                         "genes":len(v),"q1":q[0],"median":q[1],"q3":q[2]})
            x += 0.65
        bp = gax.boxplot(data, positions=pos, widths=0.72, showfliers=False, patch_artist=True,
                         medianprops={"color":"#222","lw":0.8},
                         whiskerprops={"color":"#555","lw":0.55}, capprops={"color":"#555","lw":0.55})
        for patch_, color in zip(bp["boxes"], colors): patch_.set_facecolor(color); patch_.set_alpha(0.72)
        centers = [1, 4.65, 8.3]
        gax.set_xticks(centers, ["Upstream\n2 kb", "Gene", "Downstream\n2 kb"])
        gax.set_title(analysis, fontweight="bold"); clean(gax, "y")
        if j == 0: gax.set_ylabel(r"log$_{10}$(1 + SNPs kb$^{-1}$)")
    gaxes[0].set_ylim(0, 2.55); gaxes[1].set_ylim(0, 2.55)
    gaxes[0].text(0, 1.08, "Diagnostic SNP density", transform=gaxes[0].transAxes,
                  fontweight="bold", fontsize=8.6)
    gaxes[1].legend(
        [Rectangle((0,0),1,1,color=CLASS_COLOR[c],alpha=0.72)
         for c in ["HapDom","NoDiff","Sub"]],
        [CLASS_LABEL[c] for c in ["HapDom","NoDiff","Sub"]],
        frameon=False, fontsize=5.3, loc="upper right")
    panel(gaxes[0], "G", -0.22, 1.08)
    panel_groups["01G_SNP_density"] = gaxes
    pd.DataFrame(snp_summary_rows).to_csv(OUT / "source_G_SNP_density_summary.tsv", sep="\t", index=False)

    # H: candidate allele-fraction trajectories.
    hgrid = GridSpecFromSubplotSpec(2, 2, subplot_spec=gs[2, 10:16], hspace=0.68, wspace=0.25)
    haxes = []
    for i, ((analysis, label), part_) in enumerate(candidate.groupby(["analysis", "display_label"], sort=False)):
        hax = fig.add_subplot(hgrid[i//2, i%2]); haxes.append(hax)
        z = part_.sort_values("stage_index"); x = z.stage_index.to_numpy()
        ca, cb = (FL_A, FL_B) if analysis == "FL" else (TN_A, TN_B)
        hax.plot(x, z.allele_A_fraction, "o-", color=ca, ms=2.2, lw=1.0, label="Allele A")
        hax.plot(x, z.allele_B_fraction, "o-", color=cb, ms=2.2, lw=1.0, label="Allele B")
        hax.axhline(0.5, color="#555", ls="--", lw=0.55); hax.axvspan(13.5,19.5,color="#FFF5E8",zorder=-3)
        hax.set_ylim(-0.02,1.02); hax.set_xlim(0.5,19.5); hax.set_title(f"{label} ({analysis})", fontsize=7.2, fontweight="bold")
        hax.set_xticks([1,5,9,13,16,19], [STAGES[k-1] for k in [1,5,9,13,16,19]], rotation=45, ha="right")
        clean(hax, "y")
        if i % 2 == 0: hax.set_ylabel("Allelic fraction")
    if haxes:
        haxes[0].legend(frameon=False, ncol=2, loc="lower center", bbox_to_anchor=(1.18, 1.13))
        haxes[0].text(0, 1.34, "Candidate trajectories", transform=haxes[0].transAxes,
                      fontweight="bold", fontsize=8.6)
        panel(haxes[0], "H", -0.28, 1.34)
        panel_groups["01H_candidate_trajectories"] = haxes

    fig.suptitle("Allele-specific expression landscapes of FL and TN oil palm",
                 fontsize=15, fontweight="bold", y=0.995)
    fig.text(0.995, 0.004,
             "A = Africa hap2 (FL) or Dura/TK-like (TN); B = American hap1 (FL) or Pisifera/NS-like (TN).",
             ha="right", fontsize=6.3, color="#555D66")
    save_panel_groups(fig, panel_groups)
    save_figure(fig, "01_FL_TN_ASE_landscape_reference1")


def chromosome_windows(snp: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for analysis in ["FL", "TN"]:
        part = snp[snp.analysis == analysis].copy()
        part["chrom_num"] = part.chrom.str.extract(r"chr(\d+)").astype(int)
        part["window_mb"] = ((part.start + part.end) / 2 // 1_000_000).astype(int)
        part["robust_any"] = part.overall_class != "NoASE"
        agg = (part.groupby(["chrom_num","window_mb"]).agg(
            eligible_genes=("gene_id","size"), ASE_genes=("robust_any","sum")).reset_index())
        agg["analysis"] = analysis
        rows.append(agg)
    out = pd.concat(rows, ignore_index=True)
    out.to_csv(OUT / "source2_A_ASE_genes_per_1Mb.tsv", sep="\t", index=False)
    return out


def build_upset(overlap: pd.DataFrame, phase_fl: pd.DataFrame,
                phase_tn: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    fl = phase_fl[["gene_id"] + PHASES].rename(columns={"gene_id":"gene_africa",
         **{p:f"FL_{p}" for p in PHASES}})
    tn = phase_tn[["gene_id"] + PHASES].rename(columns={"gene_id":"gene_dura",
         **{p:f"TN_{p}" for p in PHASES}})
    z = overlap.merge(fl, on="gene_africa", how="inner").merge(tn, on="gene_dura", how="inner")
    set_cols = [f"FL_{p}" for p in PHASES] + [f"TN_{p}" for p in PHASES]
    for c in set_cols: z[c] = z[c] != "NoASE"
    z["pattern"] = z[set_cols].astype(int).astype(str).agg("".join, axis=1)
    patterns = z.groupby("pattern").size().rename("orthogroups").reset_index().sort_values("orthogroups", ascending=False)
    patterns.to_csv(OUT / "source2_B_upset_intersections.tsv", sep="\t", index=False)
    set_sizes = pd.DataFrame({"set":set_cols,"orthogroups":[int(z[c].sum()) for c in set_cols]})
    set_sizes.to_csv(OUT / "source2_B_upset_set_sizes.tsv", sep="\t", index=False)
    return patterns, z, set_sizes


def build_trait_matrix(trait_family: pd.DataFrame) -> pd.DataFrame:
    d = trait_family.copy()
    family_score = d.groupby("family").robust_ASE_rows.sum().sort_values(ascending=False)
    keep_families = family_score.head(22).index
    d = d[d.family.isin(keep_families)].copy()
    module_score = (d.groupby(["family","trait_module"]).robust_ASE_rows.sum()
                    .reset_index().sort_values(["family","robust_ASE_rows"], ascending=[True,False]))
    module_choice = module_score.drop_duplicates("family")[["family","trait_module"]]
    d = d.merge(module_choice, on=["family","trait_module"], how="inner")
    rows = []
    for keys, q in d.groupby(["analysis","stage_group","trait_module","family"], observed=True):
        analysis, stage_group, module, family = keys
        tested = int(q.tested_rows.sum()); robust = int(q.robust_ASE_rows.sum())
        weight = q.tested_rows.clip(lower=1).to_numpy(float)
        median_ratio = float(np.average(q.median_log2_ratio.to_numpy(float), weights=weight))
        rows.append({"analysis":analysis,"stage_group":stage_group,"trait_module":module,
                     "family":family,"tested_rows":tested,"robust_ASE_rows":robust,
                     "median_log2_ratio":median_ratio,
                     "robust_ASE_percentage":100*robust/tested if tested else np.nan})
    d = pd.DataFrame(rows)
    abbrev = {"Oil biosynthesis & storage":"Oil", "De-novo / saturated FA":"Sat. FA",
              "Unsaturated FA":"Unsat. FA", "TAG assembly & oil body":"TAG",
              "Lipid oxidation / antioxidant":"Oxidation", "Shell / cell wall / lignin":"Shell"}
    d["row_label"] = d.trait_module.map(abbrev).fillna(d.trait_module) + " · " + d.family.astype(str)
    order = (d.groupby("row_label").robust_ASE_rows.sum().sort_values(ascending=False).index)
    d["row_order"] = pd.Categorical(d.row_label, order, ordered=True)
    d = d.sort_values("row_order")
    d.to_csv(OUT / "source2_D_trait_family_ASE.tsv", sep="\t", index=False)
    return d


def draw_chromosome_panel(ax: plt.Axes, windows: pd.DataFrame) -> None:
    max_mb = int(windows.window_mb.max()) + 1
    vmax = max(1, float(windows.ASE_genes.quantile(0.98)))
    fl_cmap = LinearSegmentedColormap.from_list("fl", ["#F4FAF7", "#009E73"])
    tn_cmap = LinearSegmentedColormap.from_list("tn", ["#FFF7F2", "#E07A5F"])
    for chrom in range(1,17):
        y = 16 - chrom
        ax.text(-5.5, y, f"Chr{chrom}", ha="right", va="center", fontsize=6.3)
        for analysis, dy, cmap in [("FL",0.14,fl_cmap),("TN",-0.14,tn_cmap)]:
            z = windows[(windows.analysis==analysis)&(windows.chrom_num==chrom)]
            lookup = dict(zip(z.window_mb, z.ASE_genes))
            length = int(z.window_mb.max())+1 if len(z) else 0
            for mb in range(length):
                ax.add_patch(Rectangle((mb, y+dy-0.11), 1.0, 0.21,
                                       facecolor=cmap(min(1, lookup.get(mb,0)/vmax)),
                                       edgecolor="none"))
        ax.plot([0,max_mb],[y,y],color="#E6E9EB",lw=0.35,zorder=-2)
    ax.set_xlim(-7,max_mb); ax.set_ylim(-0.7,15.7)
    ax.set_xlabel("Chromosome position (Mb)"); ax.set_yticks([]); clean(ax)
    ax.set_title("ASE genes within 1-Mb windows", loc="left", fontweight="bold")
    ax.legend([Rectangle((0,0),1,1,color=FL_A),Rectangle((0,0),1,1,color=TN_B)],
              ["FL", "TN"], frameon=False, ncol=2, loc="upper right")
    ax.text(0.99, 0.955, "Darker = more ASE genes", transform=ax.transAxes,
            ha="right", va="top", fontsize=5.6, color="#5F6770")


def draw_upset(parent_spec, patterns: pd.DataFrame, set_sizes: pd.DataFrame,
               fig: plt.Figure) -> list[plt.Axes]:
    ug = GridSpecFromSubplotSpec(2, 2, subplot_spec=parent_spec,
                                 height_ratios=[0.56,0.44], width_ratios=[0.30,0.70],
                                 hspace=0.04, wspace=0.05)
    blank = fig.add_subplot(ug[0,0]); blank.axis("off")
    top = fig.add_subplot(ug[0,1]); mat = fig.add_subplot(ug[1,1], sharex=top)
    size_ax = fig.add_subplot(ug[1,0], sharey=mat)
    p = patterns.head(15).reset_index(drop=True)
    x = np.arange(len(p))
    top.bar(x, p.orthogroups, color="#3C4048", width=0.72)
    top.set_ylabel("Orthogroups"); top.set_xticks([]); clean(top, "y")
    top.set_title("Top phase × variety intersections", loc="left", fontweight="bold")
    labels = [f"FL {q}" for q in ["Early","Mid","Late","Post"]] + [f"TN {q}" for q in ["Early","Mid","Late","Post"]]
    for j, pattern_ in enumerate(p.pattern):
        bits = [int(v) for v in pattern_]
        for i, bit in enumerate(bits):
            mat.scatter(j, i, s=11, color="#30343B" if bit else "#D9DDE0", zorder=3)
        on = [i for i,b in enumerate(bits) if b]
        if len(on)>1: mat.plot([j,j],[min(on),max(on)],color="#30343B",lw=0.75,zorder=2)
    mat.set_yticks(range(8)); mat.set_yticklabels([]); mat.invert_yaxis(); mat.set_xticks([])
    mat.spines[:].set_visible(False); mat.grid(axis="x",color="#ECEEEF",lw=0.35)
    sizes = set_sizes.set_index("set").reindex([f"FL_{p}" for p in PHASES] + [f"TN_{p}" for p in PHASES])
    size_ax.barh(range(8), sizes.orthogroups, color=[FL_A]*4+[TN_B]*4, height=0.56)
    size_ax.set_yticks(range(8), labels); size_ax.invert_xaxis()
    size_ax.set_xlabel("Set size", fontsize=6.2); size_ax.tick_params(axis="x", labelsize=5.5)
    size_ax.spines[["top","right","left"]].set_visible(False)
    return [top, mat, size_ax]


def draw_volcano(ax: plt.Axes, d: pd.DataFrame, analysis: str,
                 labels: dict[str,str]) -> None:
    z = d[(d.analysis==analysis)&(d.stage=="95d")&d.eligible].copy()
    z["mlog10q"] = -np.log10(z.padj.clip(lower=1e-300)).clip(upper=80)
    bg = ~z.robust_ase; sig = z.robust_ase
    ax.scatter(z.loc[bg,"log2_allele_ratio"], z.loc[bg,"mlog10q"], s=2.0,
               color="#CED3D8", alpha=0.35, rasterized=True)
    ax.scatter(z.loc[sig,"log2_allele_ratio"], z.loc[sig,"mlog10q"], s=2.5,
               color="#D73027", alpha=0.28, rasterized=True)
    ax.axvline(-0.5,color="#555",ls="--",lw=0.55); ax.axvline(0.5,color="#555",ls="--",lw=0.55)
    ax.axhline(-math.log10(0.05),color="#555",ls=":",lw=0.55)
    candidates = z[z.gene_id.isin(labels)].copy()
    fallback = z[sig].nlargest(4,"mlog10q")
    q = pd.concat([candidates,fallback]).drop_duplicates("gene_id").head(6)
    ax.scatter(q.log2_allele_ratio, q.mlog10q, s=12, color="#B2182B",
               edgecolor="white", linewidth=0.35, zorder=5)
    for row in q.itertuples():
        name = labels.get(row.gene_id, row.gene_id.split(".")[-1])
        ax.annotate(name,(row.log2_allele_ratio,row.mlog10q),xytext=(2,2),textcoords="offset points",
                    fontsize=5.2,color="#20242A")
    ax.set_xlim(-11,11); ax.set_ylim(0,82); ax.set_title(analysis, fontweight="bold")
    ax.set_xlabel("log2(allele A / allele B)"); clean(ax, "y")
    if analysis=="FL": ax.set_ylabel("−log10(FDR)")


def build_figure_two(windows: pd.DataFrame, patterns: pd.DataFrame,
                     volcano: pd.DataFrame, trait: pd.DataFrame) -> None:
    fig = plt.figure(figsize=(15.5, 9.3))
    gs = fig.add_gridspec(2, 3, width_ratios=[1.08,0.92,0.92], height_ratios=[0.94,1.06],
                          wspace=0.34, hspace=0.34)
    panel_groups: dict[str, list[plt.Axes]] = {}
    ax_a = fig.add_subplot(gs[:,0]); draw_chromosome_panel(ax_a, windows); panel(ax_a,"A",-0.10,1.02)
    panel_groups["02A_chromosome_1Mb"] = [ax_a]
    set_sizes = pd.read_csv(OUT / "source2_B_upset_set_sizes.tsv", sep="\t")
    upset_axes = draw_upset(gs[0,1], patterns, set_sizes, fig); panel(upset_axes[0],"B",-0.50,1.04)
    panel_groups["02B_phase_upset"] = upset_axes

    cg = GridSpecFromSubplotSpec(1,2,subplot_spec=gs[0,2],wspace=0.16)
    labels = {"evm.TU.chr09B.1430":"MADS34","evm.TU.chr08B.1689":"ARF2",
              "evm.TU.chr10.635":"DGAT","evm.TU.chr07.37":"ACCase"}
    caxes=[]
    for i,analysis in enumerate(["FL","TN"]):
        ax=fig.add_subplot(cg[0,i]); caxes.append(ax); draw_volcano(ax,volcano,analysis,labels)
    caxes[0].text(0,1.10,"Mid-development ASE (95 d)",transform=caxes[0].transAxes,
                  fontweight="bold",fontsize=8.6); panel(caxes[0],"C",-0.24,1.10)
    panel_groups["02C_volcano_95d"] = caxes
    volcano[(volcano.stage=="95d")&volcano.eligible].to_csv(
        OUT/"source2_C_volcano_95d.tsv",sep="\t",index=False)

    ax = fig.add_subplot(gs[1,1:3]); panel(ax,"D",-0.08,1.03)
    phase_order=["Days 0–65","Days 80–140","Days 155–185","Hours 12–72"]
    col_order=[(a,p) for a in ["FL","TN"] for p in phase_order]
    rows=trait[["row_label","row_order"]].drop_duplicates().sort_values("row_order").row_label.tolist()
    ymap={r:i for i,r in enumerate(rows)}; xmap={k:i for i,k in enumerate(col_order)}
    vals=[]
    for r in trait.itertuples():
        vals.append(abs(float(r.median_log2_ratio)))
    lim=max(1.0,float(np.quantile(vals,0.95))) if vals else 1.0
    for r in trait.itertuples():
        key=r.row_label; x=xmap[(r.analysis,r.stage_group)]; y=ymap[key]
        ax.scatter(x,y,s=5+0.55*r.robust_ASE_percentage,c=[r.median_log2_ratio],
                   cmap="RdBu_r",norm=TwoSlopeNorm(vmin=-lim,vcenter=0,vmax=lim),
                   edgecolor="#525960",linewidth=0.25)
    ax.set_xlim(-0.6,7.6); ax.set_ylim(-0.7,len(rows)-0.3); ax.invert_yaxis()
    ax.set_xticks(range(8),["D0–65","D80–140","D155–185","H12–72"]*2,rotation=35,ha="right")
    ax.set_yticks(range(len(rows)),rows); ax.tick_params(axis="y",labelsize=5.2)
    ax.axvline(3.5,color="#30343B",lw=0.8)
    ax.text(0.25,1.015,"FL",transform=ax.transAxes,ha="center",fontweight="bold")
    ax.text(0.75,1.015,"TN",transform=ax.transAxes,ha="center",fontweight="bold")
    ax.set_title("Oil, fatty-acid, oxidation and shell-related gene families",loc="left",fontweight="bold",pad=18)
    ax.grid(color="#E2E6E9",lw=0.4); ax.spines[["top","right"]].set_visible(False)
    sm=mpl.cm.ScalarMappable(norm=TwoSlopeNorm(vmin=-lim,vcenter=0,vmax=lim),cmap="RdBu_r")
    cb=fig.colorbar(sm,ax=ax,fraction=0.025,pad=0.02); cb.set_label("Median log2(A/B)")
    panel_groups["02D_trait_family_matrix"] = [ax, cb.ax]
    for s,pct in [(18,25),(32,50),(60,100)]: ax.scatter([],[],s=s,facecolor="white",edgecolor="#525960",label=f"{pct}%")
    ax.legend(title="Robust ASE",frameon=False,ncol=3,loc="lower right",bbox_to_anchor=(1.0,-0.18))

    fig.suptitle("Genomic distribution, persistence and functional context of FL/TN ASE",
                 fontsize=15,fontweight="bold",y=0.995)
    fig.text(0.995,0.004,"All ASE panels use the accepted 19-stage × 3-replicate graph-based analysis.",
             ha="right",fontsize=6.3,color="#555D66")
    save_panel_groups(fig, panel_groups)
    save_figure(fig,"02_FL_TN_ASE_genome_function_reference2")


def write_readme() -> None:
    text = """# FL/TN allele-specific expression figures\n\nThis flat directory contains two publication-oriented figure sets rebuilt panel-by-panel from the accepted current graph-based ASE analysis. No legacy expression values were used as current replicates.\n\n## Allele orientation\n\n- FL: allele A = Africa hap2; allele B = American hap1.\n- TN: allele A = Dura/TK-like; allele B = Pisifera/NS-like.\n\n## Confirmatory ASE definition\n\nEach variety has 19 stages and three biological replicates per stage. Calls require at least two qualifying replicates, pooled informative depth >=30, absolute log2 allelic ratio >=0.5, BH FDR <0.05, and the accepted mapping-bias sensitivity criterion.\n\n## Panel mapping\n\nFigure 01 follows reference image 1: A framework; B ASE-class composition; C intact-LTR profile; D stage mirror bars; E FL-to-TN alluvial; F Ka/Ks and Ks; G diagnostic SNP density; H allele-fraction trajectories.\n\nFigure 02 follows reference image 2: A 1-Mb chromosome distribution; B UpSet intersections; C 95-d volcano plots; D trait-family ASE matrix.\n\nTE panel C reports **intact LTR occupancy**, using coordinate-matched EDTA intact-LTR annotations for all four haplotypes. It is not labelled as total TE content. Ka/Ks uses one-to-one codon alignments and KaKs_Calculator GMYN. SNP density uses strand-aware 2-kb upstream, gene-body and 2-kb downstream intervals from the current diagnostic biallelic VCFs.\n\nEvery panel is separately exported as `01A`–`01H` and `02A`–`02D`; all files remain in this single flat directory. PDF and SVG are editable vectors, PNG is 600 dpi, and every quantitative panel has a matching source TSV.\n"""
    (OUT / "README.md").write_text(text, encoding="utf-8")


def main() -> None:
    set_style()
    mirror = pd.read_csv(CURRENT / "ASE_mirror_counts_current.tsv", sep="\t")
    classes = pd.read_csv(CURRENT / "gene_overall_class_current.tsv", sep="\t")
    bridge = pd.read_csv(ASE_RUN / "reviewed_four_genome_orthogroup_bridge.tsv", sep="\t")
    overlap_gene = pd.read_csv(CURRENT / "ASE_overlap_one_to_one_current.tsv", sep="\t")
    overlap_counts = pd.read_csv(CURRENT / "ASE_overlap_alluvial_counts_current.tsv", sep="\t")
    kaks = pd.read_csv(EXACT / "KaKs_by_current_ASE_class.tsv", sep="\t")
    kaks = kaks[(pd.to_numeric(kaks.Ka_Ks, errors="coerce") >= 0) &
                (pd.to_numeric(kaks.Ka_Ks, errors="coerce") <= 3) &
                (pd.to_numeric(kaks.Ks, errors="coerce") >= 0) &
                (pd.to_numeric(kaks.Ks, errors="coerce") <= 0.10)].copy()
    snp = pd.concat([
        pd.read_csv(CURRENT / "FL_diagnostic_SNP_density_current.tsv.gz", sep="\t"),
        pd.read_csv(CURRENT / "TN_diagnostic_SNP_density_current.tsv.gz", sep="\t")], ignore_index=True)
    gene_stage = pd.read_csv(ASE_RUN / "gene_stage_ASE.tsv.gz", sep="\t")
    trait_family = pd.read_csv(ASE_RUN / "trait_family_ASE_summary.tsv", sep="\t")

    te_source = OUT / "source_C_intact_LTR_profiles.tsv"
    te = pd.read_csv(te_source, sep="\t") if te_source.exists() else build_te_profiles(classes, bridge)
    candidate = build_candidate_source(gene_stage)
    build_figure_one(mirror, classes, te, overlap_counts, kaks, snp, candidate)

    windows = chromosome_windows(snp)
    phase_fl = pd.read_csv(CURRENT / "FL_four_phase_gene_classes_current.tsv", sep="\t")
    phase_tn = pd.read_csv(CURRENT / "TN_four_phase_gene_classes_current.tsv", sep="\t")
    patterns, _, _ = build_upset(overlap_gene, phase_fl, phase_tn)
    trait = build_trait_matrix(trait_family)
    build_figure_two(windows, patterns, gene_stage, trait)
    write_readme()


if __name__ == "__main__":
    main()
