#!${DATA_DIR}/miniconda3/bin/python
"""Quantify and plot the complementary value of graph-pangenome SV-GWAS.

The figure deliberately avoids claiming that SV-GWAS detects more significant
traits than SNP-GWAS. It tests the defensible claim that graph-derived SVs span
additional variant classes and add physically distinct association intervals.
"""
from __future__ import annotations

import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd


RUN = Path(__file__).resolve().parents[1]
OUT = RUN / "figures/graph_pangenome_value"
TABLES = RUN / "tables"
SV_META = Path("${ANALYSIS_DIR}/14_pan_genome/06_Minigraph/Pangenie/02_sv_combined/步骤四_过滤分类统计/sv_type_stats/sv_qc.per_sv.tsv")
TYPE_COLORS = {"DEL": "#0072B2", "INS": "#D55E00", "MNV/COMPLEX": "#7B3294"}
CLASS_COLORS = {"both": "#7B3294", "SNP_only": "#0072B2", "SV_only": "#D55E00", "neither": "#B8B8B8"}


def interval_gap(a_start: int, a_end: int, b_start: int, b_end: int) -> int:
    return max(0, b_start - a_end, a_start - b_end)


def setup_style() -> None:
    cjk = Path("/usr/share/fonts/google-noto-cjk/NotoSansCJK-Regular.ttc")
    if cjk.is_file():
        font_manager.fontManager.addfont(str(cjk))
        family = font_manager.FontProperties(fname=str(cjk)).get_name()
    else:
        family = "DejaVu Sans"
    plt.rcParams.update({
        "font.family": "sans-serif", "font.sans-serif": [family, "DejaVu Sans"],
        "axes.unicode_minus": False, "pdf.fonttype": 42, "ps.fonttype": 42,
        "font.size": 9, "axes.titlesize": 10.5, "axes.labelsize": 9.5,
    })


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    snp = pd.read_csv(RUN / "snp/summary.tsv", sep="\t")
    sv = pd.read_csv(RUN / "sv/summary.tsv", sep="\t")
    comp = snp[["category", "trait", "n_positive", "n_tests", "top_p", "bonferroni"]].merge(
        sv[["trait", "n_tests", "top_p", "bonferroni"]], on="trait", suffixes=("_snp", "_sv"), validate="one_to_one")
    comp["snp_support_log10"] = np.log10(comp.bonferroni_snp / comp.top_p_snp)
    comp["sv_support_log10"] = np.log10(comp.bonferroni_sv / comp.top_p_sv)
    comp["snp_significant"] = comp.top_p_snp < comp.bonferroni_snp
    comp["sv_significant"] = comp.top_p_sv < comp.bonferroni_sv
    comp["discovery_class"] = np.select(
        [comp.snp_significant & comp.sv_significant, comp.snp_significant, comp.sv_significant],
        ["both", "SNP_only", "SV_only"], default="neither")
    comp.to_csv(TABLES / "graph_pangenome_trait_level_comparison.tsv", sep="\t", index=False)

    loci = pd.read_csv(TABLES / "all_trait_gwas_loci.tsv", sep="\t")
    gw = loci[loci.signal_level.eq("genomewide_bonferroni")].copy()
    sv_loci = gw[gw.modality.eq("SV")].copy()
    snp_loci = gw[gw.modality.eq("SNP")]
    nearest_dist, nearest_id = [], []
    for r in sv_loci.itertuples(index=False):
        q = snp_loci[(snp_loci.trait == r.trait) & (snp_loci.chrom == r.chrom)]
        if q.empty:
            nearest_dist.append(np.inf); nearest_id.append("")
            continue
        candidates = [(interval_gap(int(r.start), int(r.end), int(x.start), int(x.end)), x.locus_id) for x in q.itertuples(index=False)]
        d, locus_id = min(candidates)
        nearest_dist.append(d); nearest_id.append(locus_id)
    sv_loci["nearest_significant_snp_distance_bp"] = nearest_dist
    sv_loci["nearest_significant_snp_locus"] = nearest_id
    sv_loci["snp_proximity_class"] = np.select(
        [sv_loci.nearest_significant_snp_distance_bp.eq(0), sv_loci.nearest_significant_snp_distance_bp.le(250_000)],
        ["overlapping SNP interval", "within 250 kb"], default=">250 kb / no same-chromosome SNP")
    sv_loci.to_csv(TABLES / "graph_pangenome_significant_sv_window_novelty.tsv", sep="\t", index=False)

    known = pd.read_csv(TABLES / "trait_relevant_known_gene_overlaps.tsv", sep="\t")
    known = known[known.signal_level.eq("genomewide_bonferroni")].copy()
    snp_pairs = set(map(tuple, known[known.modality.eq("SNP")][["trait", "gene_id"]].drop_duplicates().to_numpy()))
    sv_pairs = set(map(tuple, known[known.modality.eq("SV")][["trait", "gene_id"]].drop_duplicates().to_numpy()))
    known["modality_pair_class"] = known.apply(
        lambda r: "shared" if (r.trait, r.gene_id) in snp_pairs & sv_pairs else ("SV_only" if r.modality == "SV" else "SNP_only"), axis=1)
    sv_unique_known = known[known.apply(lambda r: r.modality == "SV" and (r.trait, r.gene_id) in sv_pairs - snp_pairs, axis=1)].copy()
    sv_unique_known.to_csv(TABLES / "graph_pangenome_sv_only_known_gene_support.tsv", sep="\t", index=False)

    meta = pd.read_csv(SV_META, sep="\t")
    meta["event_span_bp"] = meta[["ref_len", "alt_len"]].max(axis=1).clip(lower=1)
    type_stats = meta.groupby("svtype", sort=False).agg(
        n_events=("id", "size"), median_event_span_bp=("event_span_bp", "median"),
        p95_event_span_bp=("event_span_bp", lambda x: x.quantile(.95))).reset_index()
    type_stats["fraction_of_all_sv"] = type_stats.n_events / len(meta)
    type_stats.to_csv(TABLES / "graph_pangenome_sv_event_spectrum.tsv", sep="\t", index=False)

    proximity_order = ["overlapping SNP interval", "within 250 kb", ">250 kb / no same-chromosome SNP"]
    proximity = sv_loci.snp_proximity_class.value_counts().reindex(proximity_order, fill_value=0)
    metrics = pd.DataFrame([
        ("phenotypes_tested", len(comp), "traits"),
        ("snp_significant_traits", int(comp.snp_significant.sum()), "traits"),
        ("sv_significant_traits", int(comp.sv_significant.sum()), "traits"),
        ("graph_sv_events", len(meta), "events"),
        ("graph_complex_events", int((meta.svtype == "MNV/COMPLEX").sum()), "events"),
        ("graph_events_span_gt50bp", int((meta.event_span_bp > 50).sum()), "events"),
        ("graph_events_span_gt1kb", int((meta.event_span_bp > 1000).sum()), "events"),
        ("significant_sv_reporting_windows", len(sv_loci), "250-kb-clustered windows"),
        ("significant_sv_windows_gt250kb_from_snp", int((sv_loci.nearest_significant_snp_distance_bp > 250_000).sum()), "windows"),
        ("sv_trait_gene_pairs", len(sv_pairs), "unique trait-gene pairs"),
        ("sv_only_trait_gene_pairs", len(sv_pairs - snp_pairs), "unique trait-gene pairs"),
    ], columns=["metric", "value", "unit"])
    metrics.to_csv(TABLES / "graph_pangenome_value_metrics.tsv", sep="\t", index=False)

    setup_style()
    fig = plt.figure(figsize=(13.4, 10.2))
    gs = fig.add_gridspec(2, 2, width_ratios=[1.02, .98], height_ratios=[1, 1], hspace=.36, wspace=.28)

    # A: distribution of graph-derived event spans.
    ax = fig.add_subplot(gs[0, 0])
    bins = np.logspace(0, math.log10(meta.event_span_bp.max()), 48)
    for typ in ["DEL", "INS", "MNV/COMPLEX"]:
        vals = meta.loc[meta.svtype.eq(typ), "event_span_bp"]
        ax.hist(vals, bins=bins, histtype="step", lw=1.8, color=TYPE_COLORS[typ], label=f"{typ}: {len(vals):,}")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.axvline(50, color="#555", ls="--", lw=.8)
    ax.axvline(1000, color="#555", ls=":", lw=.8)
    ax.set_xlabel("图泛基因组变异事件跨度（bp，对数尺度）")
    ax.set_ylabel("事件数（对数尺度）")
    ax.set_title("A  图泛比对覆盖 SNP 难以表示的多尺度变异", loc="left", fontweight="bold")
    ax.legend(frameon=False, fontsize=8, loc="upper right")
    ax.text(.02, .04,
            f"总计 {len(meta):,} 个SV事件\n{(meta.event_span_bp > 50).mean():.1%} >50 bp；{(meta.event_span_bp > 1000).mean():.1%} >1 kb",
            transform=ax.transAxes, va="bottom", fontsize=9, bbox=dict(boxstyle="round,pad=.3", fc="white", ec="#CCCCCC", alpha=.92))
    ax.spines[["top", "right"]].set_visible(False)

    # B: fair trait-level comparison normalized by each modality threshold.
    ax = fig.add_subplot(gs[0, 1])
    for cls in ["neither", "SNP_only", "both", "SV_only"]:
        q = comp[comp.discovery_class.eq(cls)]
        ax.scatter(q.snp_support_log10, q.sv_support_log10, s=np.where(q.discovery_class.eq("both"), 46, 30),
                   color=CLASS_COLORS[cls], alpha=.82, edgecolor="white", linewidth=.45, label=f"{cls}: {len(q)}")
    ax.axhline(0, color="#444", lw=.8); ax.axvline(0, color="#444", lw=.8)
    lim = max(abs(comp[["snp_support_log10", "sv_support_log10"]].to_numpy()).max(), 1) * 1.08
    ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim)
    ax.set_xlabel(r"SNP关联强度  $\log_{10}(P_{Bonf}/P_{top})$")
    ax.set_ylabel(r"SV关联强度  $\log_{10}(P_{Bonf}/P_{top})$")
    ax.set_title("B  两类GWAS在表型层面相互验证，而非SV数量取胜", loc="left", fontweight="bold")
    ax.legend(frameon=False, fontsize=7.7, loc="lower right")
    for r in comp[comp.discovery_class.eq("both")].sort_values("sv_support_log10", ascending=False).head(4).itertuples():
        ax.annotate(r.trait.replace("_", " "), (r.snp_support_log10, r.sv_support_log10), xytext=(4, 4), textcoords="offset points", fontsize=6.8)
    ax.spines[["top", "right"]].set_visible(False); ax.grid(color="#EAEAEA", lw=.45, zorder=0)

    # C: significant SV windows relative to significant SNP windows.
    ax = fig.add_subplot(gs[1, 0])
    prox_colors = ["#5E81AC", "#EBCB8B", "#D55E00"]
    left = 0
    for label, color in zip(proximity_order, prox_colors):
        value = int(proximity[label]); ax.barh([.18], [value], left=left, height=.34, color=color, edgecolor="white", label=label)
        ax.text(left + value / 2, .18, f"{value}\n({value/len(sv_loci):.1%})", ha="center", va="center", fontsize=9, fontweight="bold")
        left += value
    ax.set_xlim(0, len(sv_loci)); ax.set_ylim(-.72, .78); ax.set_yticks([])
    ax.set_xlabel("全基因组显著SV报告区间数（同一表型内比较）")
    ax.set_title("C  40.9%的显著SV区间位于显著SNP邻域之外", loc="left", fontweight="bold")
    legend_labels = ["与显著SNP区间重叠", "距显著SNP≤250 kb", "距显著SNP>250 kb/同染色体无显著SNP"]
    handles = [Line2D([0], [0], color=c, lw=7, label=l) for c, l in zip(prox_colors, legend_labels)]
    ax.legend(handles=handles, frameon=False, fontsize=7.5, loc="lower left", bbox_to_anchor=(0, -.02), ncol=1)
    ax.text(.01, .92, f"已知基因支持：SV关联到 {len(sv_pairs)} 个表型–基因组合，其中 {len(sv_pairs-snp_pairs)} 个未被显著SNP窗口覆盖",
            transform=ax.transAxes, fontsize=8.6)
    ax.spines[["top", "right", "left"]].set_visible(False)

    # D: a concrete SV-added locus for nut length.
    ax = fig.add_subplot(gs[1, 1])
    focus = gw[(gw.trait.eq("Nut_length_mm")) & (gw.chrom.eq("chr01B")) & (gw.start <= 12_000_000)].copy()
    ypos = {"SNP": 1.0, "SV": 0.0}
    for mod, color in [("SNP", "#0072B2"), ("SV", "#D55E00")]:
        q = focus[focus.modality.eq(mod)]
        for r in q.itertuples(index=False):
            ax.plot([r.start / 1e6, r.end / 1e6], [ypos[mod], ypos[mod]], color=color, lw=5.2, solid_capstyle="butt", alpha=.82)
            ax.scatter(r.lead_pos / 1e6, ypos[mod], s=18 + 7 * max(0, -math.log10(r.lead_p) - 7), color=color, edgecolor="white", lw=.45, zorder=3)
    cinv_start, cinv_end = 9_280_465 / 1e6, 9_287_063 / 1e6
    ax.axvspan(cinv_start, cinv_end, color="#009E73", alpha=.28, lw=0)
    ax.annotate("中性转化酶 CINV1\nSV区间直接覆盖\n最近显著SNP约2.4 Mb外",
                ((cinv_start + cinv_end) / 2, .05), xytext=(8.05, .48), textcoords="data", ha="center", fontsize=8.3,
                arrowprops=dict(arrowstyle="->", color="#166B55", lw=1.0),
                bbox=dict(boxstyle="round,pad=.25", fc="white", ec="#85B8A8"))
    ax.set_xlim(.5, 12); ax.set_ylim(-.55, 1.55); ax.set_yticks([0, 1], ["SV", "SNP"])
    ax.set_xlabel("chr01B位置（Mb）")
    ax.set_title("D  坚果长：SV补充SNP峰之外的功能候选区间", loc="left", fontweight="bold")
    ax.grid(axis="x", color="#EAEAEA", lw=.45); ax.spines[["top", "right", "left"]].set_visible(False)
    ax.legend(handles=[Line2D([0], [0], color="#0072B2", lw=5, label="显著SNP报告区间"),
                       Line2D([0], [0], color="#D55E00", lw=5, label="显著SV报告区间"),
                       Line2D([0], [0], color="#009E73", lw=8, alpha=.35, label="CINV1基因")],
              frameon=False, fontsize=7.7, loc="upper right")

    fig.suptitle("图泛基因组SV扩展SNP-GWAS可见的遗传关联空间", fontsize=16, fontweight="bold", y=.985)
    fig.text(.5, .945, "证据来自互补变异尺度、独立关联区间与功能候选增益；不等同于‘SV显著表型数多于SNP’", ha="center", fontsize=10.2, color="#444")
    fig.text(.5, .012,
             "注：表型已排除精确0值及非零分布上下各2.5%极端值；SNP与SV均为EMMAX校正。报告区间按250 kb聚类，物理距离不等同于正式LD独立性。",
             ha="center", fontsize=8.3, color="#444")
    fig.subplots_adjust(top=.90, bottom=.10, left=.075, right=.98)
    fig.savefig(OUT / "graph_pangenome_sv_complements_snp_gwas.png", dpi=450, bbox_inches="tight")
    fig.savefig(OUT / "graph_pangenome_sv_complements_snp_gwas.pdf", bbox_inches="tight")
    plt.close(fig)

    report = RUN / "reports/GRAPH_PANGENOME_VALUE_INTERPRETATION.md"
    report.write_text(
        "# 图泛基因组对GWAS的增益\n\n"
        f"- 图泛比对获得 {len(meta):,} 个SV事件，其中 {(meta.svtype == 'MNV/COMPLEX').sum():,} 个为MNV/COMPLEX；"
        f"{(meta.event_span_bp > 50).mean():.1%} 的事件跨度>50 bp，{(meta.event_span_bp > 1000).mean():.1%}>1 kb。\n"
        f"- 60个表型中SNP-GWAS有{int(comp.snp_significant.sum())}个、SV-GWAS有{int(comp.sv_significant.sum())}个达到全基因组Bonferroni显著；"
        f"SV显著的{int(comp.sv_significant.sum())}个表型均与SNP结果重合，因此不支持‘SV发现更多显著表型’。\n"
        f"- 115个显著SV报告区间中，{int((sv_loci.nearest_significant_snp_distance_bp > 250_000).sum())}个（{(sv_loci.nearest_significant_snp_distance_bp > 250_000).mean():.1%}）"
        "距同一表型的显著SNP区间>250 kb，支持SV提供额外关联区间。\n"
        f"- SV显著区间关联到{len(sv_pairs)}个有文献支持的表型–基因组合，其中{len(sv_pairs-snp_pairs)}个未被显著SNP窗口覆盖；"
        "该组合为坚果长–中性转化酶CINV1。\n\n"
        "推荐表述：‘Graph-pangenome-derived SVs complement SNP-GWAS by interrogating multi-base and complex variation and by revealing significant association intervals outside SNP-defined peaks.’\n"
        "不推荐表述：‘SV-GWAS detects more significant traits than SNP-GWAS.’\n",
        encoding="utf-8")
    print(metrics.to_string(index=False))


if __name__ == "__main__":
    main()
