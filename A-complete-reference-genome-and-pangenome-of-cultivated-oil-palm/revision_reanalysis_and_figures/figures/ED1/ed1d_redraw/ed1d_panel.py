#!/usr/bin/env python3
"""ED1d redraw: HiFi and ONT read alignments around the two FL terminal-extension junctions.

Data (cluster work dir ${CLUSTER_WORK}/ed1d_redraw, copied to ./data):
  reads   = user 01_seedless_results/12_evaluate/{02_american,01_africa}/02_remapping/{hifi,ont}_remap.sorted.bam
            (minimap2 2.26 map-hifi / map-ont, --secondary=no, against American_hap1.fa / Africa_hap2.fa, which are
            identical in these windows to the final FL-Hap1 / FL-Hap2 chromosomes; seq/refcheck.tsv)
  tables  = s03_tables.py (pysam): per-alignment table (primary + supplementary; secondary/unmapped dropped) and
            100-bp mean depth (all alignments; MAPQ >= 20)
  telomere= (CCCTAAA)n / (TTTAGGG)n density from the final assembly sequence (seq/*.fa)

Usage:
  python3 ed1d_panel.py            -> ED1d_panel.pdf/.png (standalone, same size as panel d in ED1) + Source Data TSVs
  from ed1d_panel import draw_d    -> draw into an existing figure (used by plot_ED1_newd.py)
"""
import re
import sys
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
DATA = HERE / "data"
sys.path.insert(0, str(HERE.parents[1] / "fix/beautify/common") if (HERE.parents[1] / "fix").exists()
                else str(HERE.parent / "beautify/common"))
from palA import tint  # noqa: E402

MM = 1 / 25.4
MAPQ_MIN = 20
FLANK = 1000          # an alignment spans the junction when it extends >= 1 kb on both sides
CLIP_MIN = 500        # soft/hard clip >= 500 bp is marked
N_ROWS = 12           # read rows shown per platform

# ED1c colours (HiFi red, ONT green, filled gap/added sequence brown)
COL = {"HiFi": "#F28080", "ONT": "#4DBFAD"}
COL_DARK = {"HiFi": "#B8474A", "ONT": "#23806F"}      # minus strand: same hue, darker
C_ADDED, C_DRAFT, C_TELO, C_LOWQ = "#C07A45", "#9A9A9A", "#1F1F1F", "#D0D0D0"

LOCI = [
    # junctions from s06_junction.sh (terminal 150 kb of the Verkko draft aligned to the final ends, minimap2 asm5):
    # draft hap1_10:15,426-150,000 -> chr12A:75,516-210,114 (+); draft hap2_29:129,253,371-129,399,336 -> chr07B:
    # 128,856,911-129,402,880 (+). Sequence outside these blocks is new in the final assembly.
    dict(key="chr12A_left", chrom="chr12A", hap="FL-Hap1", start=1, end=130000, junction=75515,
         added=(1, 75515), side="left", fasta="chr12A_head.fa", fa_off=1, unit="kb", ticks=[1] + list(range(20000, 130001, 20000)),
         title="FL-Hap1 chr12A, left-end extension (75.5 kb)"),
    dict(key="chr07B_right", chrom="chr07B", hap="FL-Hap2", start=129352881, end=129452880, junction=129402880,
         added=(129402881, 130504653), side="right", fasta="chr07B_tail.fa", fa_off=128400000, unit="Mb",
         ticks=range(129360000, 129460001, 20000), title="FL-Hap2 chr07B, right-end extension (1.10 Mb)"),
]


# ------------------------------------------------------------------ data
def read_fa(p):
    return "".join(l.strip() for l in open(p) if not l.startswith(">")).upper()


def telomere_bins(loc, step=200, min_motif=10):
    """(start, end) runs where >= min_motif telomere heptamers occur per `step` bp (plant TTTAGGG)."""
    s = read_fa(DATA / loc["fasta"])
    off = loc["fa_off"]
    a, b = loc["start"] - off, loc["end"] - off + 1
    runs = []
    for i in range(a, b, step):
        w = s[i:i + step]
        n = len(re.findall("CCCTAAA", w)) + len(re.findall("TTTAGGG", w))
        if n >= min_motif:
            p0, p1 = off + i, off + i + len(w) - 1
            if runs and p0 - runs[-1][1] <= step + 1:
                runs[-1][1] = p1
            else:
                runs.append([p0, p1])
    # chromosome-end telomere outside the window (for the annotation)
    end_runs = []
    if loc["side"] == "right":
        for i in range(len(s) - 20000, len(s), step):
            w = s[i:i + step]
            if len(re.findall("CCCTAAA", w)) + len(re.findall("TTTAGGG", w)) >= min_motif:
                end_runs.append(off + i)
    return runs, (min(end_runs), off + len(s) - 1) if end_runs else None


def load(loc, pl):
    r = pd.read_csv(DATA / f"tables/{loc['key']}.{pl.lower()}.reads.tsv", sep="\t")
    d = pd.read_csv(DATA / f"tables/{loc['key']}.{pl.lower()}.depth100.tsv", sep="\t")
    J = loc["junction"]
    r["span"] = (r.ref_start_1b <= J - FLANK + 1) & (r.ref_end_1b >= J + FLANK)
    return r, d


def span_counts(r):
    sp = r[r.span]
    return dict(n=len(sp), n_q=int((sp.mapq >= MAPQ_MIN).sum()), n_reads=sp.read.nunique(),
                n_reads_q=sp[sp.mapq >= MAPQ_MIN].read.nunique())


def pack(r, loc, rows=N_ROWS, gap_frac=0.004, seed=7):
    """Greedy row packing: spanning alignments first (MAPQ-high first), then a random sample ordered by start."""
    width = loc["end"] - loc["start"] + 1
    gap = gap_frac * width
    rng = np.random.default_rng(seed)
    vis = r[(r.ref_end_1b >= loc["start"]) & (r.ref_start_1b <= loc["end"])].copy()
    vis["x0"] = vis.ref_start_1b.clip(lower=loc["start"])
    vis["x1"] = vis.ref_end_1b.clip(upper=loc["end"])
    first = vis[vis.span].sort_values(["mapq", "ref_start_1b"], ascending=[False, True])
    rest = vis[~vis.span]
    rest = rest.iloc[rng.permutation(len(rest))]
    occ = [[] for _ in range(rows)]          # occupied intervals per row
    placed, tried = [], 0
    for part in (first, rest):
        for a in part.itertuples(index=False):
            tried += 1
            if tried > 20000 and part is rest:
                break
            for k in range(rows):
                if all(a.x1 + gap < b0 or a.x0 > b1 + gap for b0, b1 in occ[k]):
                    occ[k].append((a.x0, a.x1))
                    placed.append((k, a))
                    break
    return placed, len(vis)


# ------------------------------------------------------------------ drawing
def _fmt_x(loc):
    if loc["unit"] == "kb":
        return mpl.ticker.FuncFormatter(lambda v, _: f"{round(v / 1e3):.0f}")
    return mpl.ticker.FuncFormatter(lambda v, _: f"{v / 1e6:.2f}")


def draw_locus(fig, box, loc, W, H, log=None):
    """box = (x_mm, y_top_mm, w_mm, h_mm) of the track area (title drawn by caller)."""
    x, y, w, h = box
    lab_w = 9.0            # left label column
    ann_w = 12.5           # right annotation column (spanning counts)
    tw = w - lab_w - ann_w
    heights = [("struct", 1.4), ("gap", 0.5), ("cov_HiFi", 2.8), ("reads_HiFi", 3.9), ("gap", 0.6),
               ("cov_ONT", 2.8), ("reads_ONT", 3.9)]
    fixed = sum(v for _, v in heights)
    axis_h = h - fixed      # remaining for x-axis ticks + label
    J = loc["junction"]
    X0, X1 = loc["start"] - 0.5, loc["end"] + 0.5
    axes, ymid = {}, {}
    yy = y
    for name, hh in heights:
        if name != "gap":
            axes[name] = fig.add_axes([(x + lab_w) / W, 1 - (yy + hh) / H, tw / W, hh / H])
            axes[name].set_xlim(X0, X1)
            ymid[name] = yy + hh / 2
        yy += hh
    stats = {}

    # structure track
    ax = axes["struct"]
    ax.set_ylim(0, 1); ax.axis("off")
    a0, a1 = loc["added"]
    d0, d1 = (J + 1, loc["end"]) if loc["side"] == "left" else (loc["start"], J)
    ax.add_patch(mpl.patches.Rectangle((d0, 0.15), d1 - d0 + 1, 0.7, color=C_DRAFT, lw=0))
    ax.add_patch(mpl.patches.Rectangle((max(a0, loc["start"]), 0.15), min(a1, loc["end"]) - max(a0, loc["start"]) + 1,
                                       0.7, color=C_ADDED, lw=0))
    runs, end_telo = telomere_bins(loc)
    for p0, p1 in runs:
        ax.add_patch(mpl.patches.Rectangle((p0, 0.0), max(p1 - p0 + 1, (X1 - X0) * 0.006), 1.0, color=C_TELO, lw=0))
    stats["telomere_runs_in_window"] = ";".join(f"{p0}-{p1}" for p0, p1 in runs) or "none"
    stats["chromosome_end_telomere"] = f"{end_telo[0]}-{end_telo[1]}" if end_telo else ""
    if loc["side"] == "right":
        ax.annotate("", xy=(X1 + (X1 - X0) * 0.035, 0.5), xytext=(X1, 0.5), annotation_clip=False,
                    arrowprops=dict(arrowstyle="-|>", lw=0.5, color=C_ADDED, mutation_scale=4, shrinkA=0, shrinkB=0))
        fig.text((x + lab_w + tw + 1.8) / W, 1 - (y + 0.75) / H,
                 f"to telomere\n({end_telo[0] / 1e6:.2f} Mb)" if end_telo else "to chr. end",
                 fontsize=5, va="center", ha="left", linespacing=0.95)

    # label column
    def lab(name, text, color="black", bold=False):
        a = axes[name]
        a.text(-0.012, 0.5, text, transform=a.transAxes, ha="right", va="center", fontsize=5, color=color,
               fontweight="bold" if bold else "normal", linespacing=0.95)
    lab("struct", "Assembly")

    for pl in ("HiFi", "ONT"):
        r, d = load(loc, pl)
        sc = span_counts(r)
        stats[pl] = sc
        # coverage (log10 depth + 1 so zero depth plots at the baseline)
        ax = axes[f"cov_{pl}"]
        xs = np.repeat(np.r_[d.bin_start_1b.values, d.bin_end_1b.values[-1] + 1], 2)[1:-1]
        ya = np.repeat(d.mean_depth_all.values, 2)
        yq = np.repeat(d.mean_depth_mapq20.values, 2)
        ax.fill_between(xs, 1, ya + 1, color=tint(COL[pl], 0.35), lw=0)
        ax.fill_between(xs, 1, yq + 1, color=COL[pl], lw=0)
        ax.set_yscale("log"); ax.set_ylim(1, 2e4)
        ax.set_yticks([1, 1001], ["0", "10³"]); ax.minorticks_off()
        ax.tick_params(axis="y", labelsize=5, length=1.2, width=0.4, pad=0.8)
        ax.set_xticks([]); ax.spines[["top", "right", "bottom"]].set_visible(False)
        ax.spines["left"].set_linewidth(0.4)
        ax.axhline(1, color="#606060", lw=0.3)
        ax.text(-0.075, 0.5, pl, transform=ax.transAxes, ha="right", va="center", fontsize=5.5,
                color=COL_DARK[pl], fontweight="bold")
        # reads
        ax = axes[f"reads_{pl}"]
        placed, n_vis = pack(r, loc)
        stats[pl]["alignments_in_window"] = n_vis
        stats[pl]["alignments_shown"] = len(placed)
        rows = N_ROWS
        ax.set_ylim(rows - 0.4, -0.6); ax.axis("off")
        segs, cols = [], []
        clipx, clipy = [], []
        for k, a in placed:
            c = C_LOWQ if a.mapq < MAPQ_MIN else (COL_DARK[pl] if a.strand == "-" else COL[pl])
            segs.append([(a.x0, k), (a.x1, k)]); cols.append(c)
            if a.left_softclip >= CLIP_MIN and a.ref_start_1b >= loc["start"]:
                clipx.append(a.ref_start_1b); clipy.append(k)
            if a.right_softclip >= CLIP_MIN and a.ref_end_1b <= loc["end"]:
                clipx.append(a.ref_end_1b); clipy.append(k)
        ax.add_collection(mpl.collections.LineCollection(segs, colors=cols, linewidths=0.75, capstyle="butt"))
        ax.scatter(clipx, clipy, marker="|", s=2.2, lw=0.35, color="black", zorder=3)
        lab(f"reads_{pl}", "Reads")
        axes[f"reads_{pl}"].text(-0.012, 0.5, "", transform=axes[f"reads_{pl}"].transAxes)
        # spanning-read annotation (right column)
        fig.text((x + lab_w + tw + 1.2) / W, 1 - ymid[f"reads_{pl}"] / H,
                 f"{sc['n_reads']} spanning\n({sc['n_reads_q']} MAPQ ≥ {MAPQ_MIN})", fontsize=5, va="center",
                 ha="left", linespacing=0.95, color="black")

    # junction line across all tracks
    for name, a in axes.items():
        a.axvline(J + 0.5, color="#303030", lw=0.4, ls=(0, (2, 1.4)), zorder=5)
    # x axis on the bottom track
    axx = fig.add_axes([(x + lab_w) / W, 1 - (y + fixed + 0.01) / H, tw / W, 0.01 / H])
    axx.set_xlim(X0, X1); axx.set_yticks([])
    for s in ("left", "right", "top"):
        axx.spines[s].set_visible(False)
    axx.spines["bottom"].set_linewidth(0.4)
    axx.set_xticks([t for t in loc["ticks"] if X0 <= t <= X1]); axx.set_xlim(X0, X1)
    axx.xaxis.set_major_formatter(_fmt_x(loc))
    axx.tick_params(axis="x", labelsize=5, length=1.5, width=0.4, pad=0.8)
    axx.set_xlabel(f"Position on {loc['chrom']} ({loc['unit']})", fontsize=5, labelpad=0.8)
    stats["axis_h_mm"] = round(axis_h, 2)
    if log is not None:
        log.append((loc, stats))
    return stats


def legend_row(fig, x, y_mm, W, H):
    """One-line key (drawn by caller at the panel-title baseline)."""
    items = [("line2", (COL["HiFi"], COL["ONT"]), "+ strand"), ("line2", (COL_DARK["HiFi"], COL_DARK["ONT"]), "− strand"),
             ("line", C_LOWQ, f"MAPQ < {MAPQ_MIN}"), ("tick", "black", "clip ≥ 0.5 kb"),
             ("box", C_DRAFT, "draft-derived"), ("box", C_ADDED, "added in finishing"), ("box", C_TELO, "(TTTAGGG)n"),
             ("dash", "#303030", "junction")]
    ax = fig.add_axes([0, 0, 1, 1], facecolor="none"); ax.set_xlim(0, W); ax.set_ylim(H, 0); ax.axis("off")
    xx = x
    for kind, c, t in items:
        if kind == "line2":
            ax.plot([xx, xx + 1.2], [y_mm, y_mm], color=c[0], lw=0.75, solid_capstyle="butt")
            ax.plot([xx + 1.2, xx + 2.4], [y_mm, y_mm], color=c[1], lw=0.75, solid_capstyle="butt")
        elif kind == "line":
            ax.plot([xx, xx + 2.4], [y_mm, y_mm], color=c, lw=0.75, solid_capstyle="butt")
        elif kind == "tick":
            ax.plot([xx + 1.2, xx + 1.2], [y_mm - 0.5, y_mm + 0.5], color=c, lw=0.35)
        elif kind == "box":
            ax.add_patch(mpl.patches.Rectangle((xx, y_mm - 0.45), 2.4, 0.9, color=c, lw=0))
        else:
            ax.plot([xx + 1.2, xx + 1.2], [y_mm - 0.7, y_mm + 0.7], color=c, lw=0.4, ls=(0, (2, 1.4)))
        tt = ax.text(xx + 2.9, y_mm, t, fontsize=5, va="center", ha="left")
        bb = tt.get_window_extent(renderer=fig.canvas.get_renderer())
        xx += 2.9 + bb.width / fig.dpi * 25.4 + 2.0
    return xx


def draw_d(fig, W, H, d_top, x0=3, gap=4, letter=True):
    """Panel d as laid out in ED1: two loci side by side; returns per-locus stats."""
    dw = (W - 6 - gap) / 2
    log = []
    for k, loc in enumerate(LOCI):
        bx = x0 + k * (dw + gap)
        fig.text(bx / W, 1 - (d_top + 2.6) / H, loc["title"], fontsize=6, va="bottom")
        draw_locus(fig, (bx, d_top + 5.4, dw, H - d_top - 5.4 - 0.2), loc, W, H, log)
    legend_row(fig, x0 + 7.5, d_top + 4.1, W, H)
    return log


def write_source_data(log, out):
    rows, cov = [], []
    for loc, st in log:
        for pl in ("HiFi", "ONT"):
            s = st[pl]
            rows.append(dict(locus=loc["key"], chromosome=loc["chrom"], haplotype=loc["hap"], platform=pl,
                             junction=f"{loc['chrom']}:{loc['junction']}|{loc['junction'] + 1}",
                             added_sequence=f"{loc['chrom']}:{loc['added'][0]:,}-{loc['added'][1]:,}",
                             window=f"{loc['chrom']}:{loc['start']:,}-{loc['end']:,}",
                             spanning_reads_all=s["n_reads"], spanning_reads_MAPQ20=s["n_reads_q"],
                             spanning_alignments_all=s["n"], spanning_alignments_MAPQ20=s["n_q"],
                             alignments_in_window=s["alignments_in_window"], alignments_drawn=s["alignments_shown"],
                             telomere_runs_in_window=st["telomere_runs_in_window"],
                             chromosome_end_telomere=st["chromosome_end_telomere"]))
            d = pd.read_csv(DATA / f"tables/{loc['key']}.{pl.lower()}.depth100.tsv", sep="\t")
            d.insert(0, "platform", pl); d.insert(0, "chromosome", loc["chrom"]); d.insert(0, "locus", loc["key"])
            cov.append(d)
    note = (f"# spanning = primary or supplementary alignment (secondary excluded) covering >= {FLANK} bp on both "
            f"sides of the junction; reads counted once; MAPQ threshold {MAPQ_MIN}\n")
    with open(out / "SourceData_ED1d_junction_spanning_reads.tsv", "w") as fh:
        fh.write(note); pd.DataFrame(rows).to_csv(fh, sep="\t", index=False)
    pd.concat(cov).to_csv(out / "SourceData_ED1d_coverage_100bp.tsv", sep="\t", index=False)
    return pd.DataFrame(rows)


if __name__ == "__main__":
    mpl.rcParams.update({"font.family": "Arial", "font.size": 5, "axes.linewidth": 0.4, "pdf.fonttype": 42,
                         "axes.unicode_minus": True})
    W, H = 180, 31
    fig = plt.figure(figsize=(W * MM, H * MM))
    log = draw_d(fig, W, H, d_top=1.0)
    fig.savefig(HERE / "ED1d_panel.pdf", dpi=600, facecolor="white")
    fig.savefig(HERE / "ED1d_panel.png", dpi=600, facecolor="white")
    print(write_source_data(log, HERE).to_string())
