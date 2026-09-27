#!${DATA_DIR}/miniconda3/bin/python
"""Aggregate all completed SV EMMAX scans and draw Manhattan/QQ panels."""
from __future__ import annotations

import csv
import heapq
import math
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import chi2

RUN = Path(__file__).resolve().parents[1]
META = Path("${ANALYSIS_DIR}/14_pan_genome/06_Minigraph/Pangenie/02_sv_combined/步骤四_过滤分类统计/sv_type_stats/sv_qc.per_sv.tsv")
CHROMS = [f"chr{i:02d}B" for i in range(1, 17)]


def sample_n(category: str, trait: str) -> int:
    audit = pd.read_csv(RUN / "manifests/phenotype_filter_audit.tsv", sep="\t")
    hit = audit[(audit.category == category) & (audit.trait == trait)]
    return int(hit.iloc[0].n_retained)


def plot_trait(category: str, trait: str, ps: Path, meta: pd.DataFrame) -> dict:
    out = RUN / "sv/results" / category / trait
    out.mkdir(parents=True, exist_ok=True)
    expected = len(meta)
    threshold = .05 / expected
    stride = max(1, expected // 100_000)
    plot_stride = max(1, expected // 80_000)
    plot_points, qq, sig_rows, top = [], [], [], []
    minp, top_id, n = 1.0, "", 0
    meta_id = meta["id"].to_numpy(); meta_chrom = meta["chrom"].to_numpy(); meta_pos = meta["pos"].to_numpy()
    meta_type = meta["svtype"].to_numpy(); meta_len = meta["svlen"].to_numpy()
    with ps.open() as fh:
        for i, line in enumerate(fh):
            x = line.split()
            if len(x) < 4:
                continue
            sv_id, beta, se, pv = x[0], float(x[-3]), float(x[-2]), max(float(x[-1]), np.nextafter(0, 1))
            if i >= expected or sv_id != meta_id[i]:
                raise RuntimeError(f"{trait}: SV order mismatch at row {i + 1}: {sv_id}")
            chrom, pos, svtype, svlen = meta_chrom[i], int(meta_pos[i]), meta_type[i], meta_len[i]
            if pv < minp:
                minp, top_id = pv, sv_id
            item = (-pv, sv_id, chrom, pos, beta, se, svtype, svlen)
            if len(top) < 1000:
                heapq.heappush(top, item)
            elif pv < -top[0][0]:
                heapq.heapreplace(top, item)
            if i % plot_stride == 0 or pv < 1e-5:
                plot_points.append((chrom, pos, pv, svtype))
            if i % stride == 0:
                qq.append(pv)
            if pv < threshold:
                sig_rows.append((sv_id, chrom, pos, svtype, svlen, beta, se, pv, threshold))
            n += 1
    if n != expected:
        raise RuntimeError(f"{trait}: expected {expected} SV rows, observed {n}")
    cols = ["SV", "CHR", "BP", "SVTYPE", "SVLEN", "BETA", "SE", "P", "BONFERRONI"]
    pd.DataFrame(sig_rows, columns=cols).to_csv(out / "significant_svs.tsv", sep="\t", index=False)
    top_rows = sorted([(-r[0], *r[1:]) for r in top])
    pd.DataFrame(top_rows, columns=["P", "SV", "CHR", "BP", "BETA", "SE", "SVTYPE", "SVLEN"]).to_csv(out / "top1000_svs.tsv", sep="\t", index=False)

    chrom_max = meta.groupby("chrom").pos.max().to_dict()
    offsets, cursor, gap = {}, 0, 2_000_000
    for chrom in CHROMS:
        offsets[chrom] = cursor
        cursor += int(chrom_max.get(chrom, 0)) + gap
    colors = {"DEL": "#D95F4A", "INS": "#3973B7", "MNV/COMPLEX": "#7B4F9D"}

    def draw_manhattan(ax):
        for svtype in ["DEL", "INS", "MNV/COMPLEX"]:
            pts = [(offsets[c] + bp, -math.log10(p)) for c, bp, p, t in plot_points if t == svtype]
            if pts:
                xx, yy = zip(*pts)
                ax.scatter(xx, yy, s=3.0 if svtype == "MNV/COMPLEX" else 2.2, color=colors[svtype], alpha=.68, lw=0, rasterized=True, label=svtype)
        ax.axhline(-math.log10(threshold), color="#B22222", ls="--", lw=.9, label="Bonferroni 0.05/N")
        centers = [offsets[c] + int(chrom_max.get(c, 0)) / 2 for c in CHROMS]
        ax.set_xticks(centers, CHROMS, rotation=45, ha="right", fontsize=7)
        ax.set_xlabel("Chromosome"); ax.set_ylabel("−log10(P)")
        ax.spines[["top", "right"]].set_visible(False); ax.grid(axis="y", color="#E6E6E6", lw=.5)

    q = np.sort(np.asarray(qq, dtype=float))
    obs = -np.log10(q); exp = -np.log10((np.arange(1, len(q) + 1) - .5) / len(q))
    lim = max(float(exp.max()), float(obs.max()))
    lam = float(np.median(chi2.isf(np.clip(q, np.nextafter(0, 1), 1), 1)) / chi2.ppf(.5, 1))

    fig, (ax, qax) = plt.subplots(1, 2, figsize=(15.2, 4.8), gridspec_kw={"width_ratios": [4.3, 1]})
    draw_manhattan(ax); ax.set_title(f"{trait}: zero- and 2.5% tail-excluded SV-GWAS", loc="left", fontsize=11)
    qax.scatter(exp, obs, s=5, color="#365F91", alpha=.55, lw=0, rasterized=True)
    qax.plot([0, lim], [0, lim], color="#B22222", lw=.8, ls="--")
    qax.set_xlabel("Expected −log10(P)"); qax.set_ylabel("Observed −log10(P)"); qax.set_title(f"QQ; λGC={lam:.3f}", fontsize=10)
    qax.spines[["top", "right"]].set_visible(False); fig.tight_layout()
    fig.savefig(out / "manhattan_qq.png", dpi=420); fig.savefig(out / "manhattan_qq.pdf"); plt.close(fig)

    mfig, maxis = plt.subplots(figsize=(14.2, 4.8)); draw_manhattan(maxis)
    maxis.set_title(f"{trait}: zero- and 2.5% tail-excluded SV-GWAS", loc="left", fontsize=11)
    maxis.legend(frameon=False, fontsize=7, ncol=4); mfig.tight_layout()
    mfig.savefig(out / "manhattan.png", dpi=420); mfig.savefig(out / "manhattan.pdf"); plt.close(mfig)
    qfig, qaxis = plt.subplots(figsize=(5.4, 5.0)); qaxis.scatter(exp, obs, s=6, color="#365F91", alpha=.55, lw=0, rasterized=True)
    qaxis.plot([0, lim], [0, lim], color="#B22222", lw=.8, ls="--"); qaxis.set_xlabel("Expected −log10(P)"); qaxis.set_ylabel("Observed −log10(P)")
    qaxis.set_title(f"{trait}: SV-GWAS QQ; λGC={lam:.3f}", fontsize=10); qaxis.spines[["top", "right"]].set_visible(False); qfig.tight_layout()
    qfig.savefig(out / "qq.png", dpi=420); qfig.savefig(out / "qq.pdf"); plt.close(qfig)
    result = {"category": category, "trait": trait, "n_positive": sample_n(category, trait), "n_tests": n,
            "top_sv": top_id, "top_p": minp, "bonferroni": threshold, "n_significant": len(sig_rows), "lambda_gc": lam, "status": "complete"}
    pd.DataFrame([result]).to_csv(out / "summary_row.tsv", sep="\t", index=False)
    return result


def main(task_id=None):
    meta = pd.read_csv(META, sep="\t", usecols=["chrom", "pos", "id", "svtype", "svlen"])
    tasks = pd.read_csv(RUN / "manifests/sv_tasks.tsv", sep="\t")
    selected = tasks.iloc[[task_id - 1]] if task_id is not None else tasks
    rows = []
    for r in selected.itertuples(index=False):
        ps = Path(str(r.out_prefix) + ".ps")
        if not ps.is_file() or ps.stat().st_size == 0:
            raise FileNotFoundError(ps)
        rows.append(plot_trait(r.category, r.trait, ps, meta))
        print(f"sv_plot_done trait={r.trait}", flush=True)
    if task_id is None:
        pd.DataFrame(rows).to_csv(RUN / "sv/summary.tsv", sep="\t", index=False)
    print(f"sv_traits={len(rows)} expected_svs={len(meta)}")


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else None)
