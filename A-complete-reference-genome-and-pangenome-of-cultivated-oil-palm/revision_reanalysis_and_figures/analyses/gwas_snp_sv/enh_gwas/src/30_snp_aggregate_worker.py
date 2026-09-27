#!${DATA_DIR}/miniconda3/bin/python
"""Aggregate and plot one zero- and phenotype-tail-excluded SNP-GWAS trait."""
from __future__ import annotations

import csv
import fcntl
import heapq
import math
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import chi2


RUN = Path(__file__).resolve().parents[1]
TARGET_START, TARGET_END = 15_600_000, 16_200_000


def main(task_id: int):
    with (RUN / "manifests/sv_tasks.tsv").open() as fh:
        traits = list(csv.DictReader(fh, delimiter="\t"))
    trait_row = traits[task_id - 1]
    category, trait = trait_row["category"], trait_row["trait"]
    with (RUN / "manifests/snp_tasks.tsv").open() as fh:
        rows = [r for r in csv.DictReader(fh, delimiter="\t") if r["category"] == category and r["trait"] == trait]
    rows.sort(key=lambda r: r["chrom"])
    if len(rows) != 16:
        raise RuntimeError(f"{trait}: expected 16 chromosome tasks, found {len(rows)}")
    expected_total = sum(int(r["expected_rows"]) for r in rows)
    threshold = 0.05 / expected_total
    sample_stride = max(1, expected_total // 100_000)
    plot_stride = max(1, expected_total // 80_000)
    top_heap = []
    bin_best = {}
    sampled = []
    qq = []
    chrom_max = {}
    complete_chr = 0
    n = 0
    minp, top_marker = 1.0, ""
    outdir = RUN / "snp/results" / category / trait
    outdir.mkdir(parents=True, exist_ok=True)
    lock_handle = (outdir / "aggregate.flock").open("w")
    fcntl.flock(lock_handle, fcntl.LOCK_EX)
    ready = outdir / "summary_row.tsv"
    if ready.is_file() and (outdir / "manhattan.png").is_file() and (outdir / "qq.png").is_file():
        print(f"existing complete aggregate trait={trait}")
        return
    region_path = outdir / "chr01B_15p6_16p2_snps.tsv"
    sig_path = outdir / "significant_snps.tsv"
    with region_path.open("w") as region, sig_path.open("w") as sig:
        region.write("SNP\tCHR\tBP\tBETA\tSE\tP\n")
        sig.write("SNP\tCHR\tBP\tBETA\tSE\tP\tBONFERRONI\n")
        for row in rows:
            chrom = row["chrom"]
            marker_index = RUN / "manifests" / f"snp_markers_chr{chrom}.tsv"
            ps = Path(row["out_prefix"] + ".ps")
            expected = int(row["expected_rows"])
            if not ps.is_file() or ps.stat().st_size == 0:
                continue
            with ps.open() as pf:
                observed = sum(1 for _ in pf)
            if observed != expected:
                raise RuntimeError(f"{trait} chr{chrom}: expected {expected}, observed {observed}")
            complete_chr += 1
            if not marker_index.is_file():
                raise FileNotFoundError(marker_index)
            with marker_index.open() as tf, ps.open() as pf:
                for tl, pl in zip(tf, pf):
                    t, p = tl.split(), pl.split()
                    marker, bp = t[0], int(t[1])
                    beta, se, pv = float(p[-3]), float(p[-2]), float(p[-1])
                    pv = max(pv, np.nextafter(0, 1))
                    chrom_max[chrom] = max(chrom_max.get(chrom, 0), bp)
                    if pv < minp:
                        minp, top_marker = pv, marker
                    item = (-pv, marker, chrom, bp, beta, se)
                    if len(top_heap) < 1000:
                        heapq.heappush(top_heap, item)
                    elif pv < -top_heap[0][0]:
                        heapq.heapreplace(top_heap, item)
                    key = (chrom, bp // 50_000)
                    if key not in bin_best or pv < bin_best[key][0]:
                        bin_best[key] = (pv, marker, bp)
                    if n % plot_stride == 0:
                        sampled.append((chrom, bp, pv))
                    if n % sample_stride == 0:
                        qq.append(pv)
                    if chrom == "01B" and TARGET_START <= bp <= TARGET_END:
                        region.write(f"{marker}\tchr01B\t{bp}\t{beta}\t{se}\t{pv}\n")
                    if pv < threshold:
                        sig.write(f"{marker}\tchr{chrom}\t{bp}\t{beta}\t{se}\t{pv}\t{threshold}\n")
                    n += 1
    if complete_chr != 16 or n != expected_total:
        raise RuntimeError(f"{trait}: incomplete aggregate chromosomes={complete_chr}/16 rows={n}/{expected_total}")
    status = "complete"
    threshold = 0.05 / n if n else math.nan
    top_rows = sorted([(-x[0], x[1], x[2], x[3], x[4], x[5]) for x in top_heap])
    with (outdir / "top1000_snps.tsv").open("w") as fh:
        fh.write("P\tSNP\tCHR\tBP\tBETA\tSE\tBONFERRONI_SIGNIFICANT\n")
        for pv, marker, chrom, bp, beta, se in top_rows:
            fh.write(f"{pv}\t{marker}\tchr{chrom}\t{bp}\t{beta}\t{se}\t{pv < threshold}\n")
    # Add the strongest point per 50-kb bin so sharp signals are retained.
    points = set(sampled)
    points.update((chrom, bp, pv) for (chrom, _), (pv, _, bp) in bin_best.items())
    chroms = [f"{i:02d}B" for i in range(1, 17)]
    gap = 2_000_000
    offsets, cursor = {}, 0
    for chrom in chroms:
        offsets[chrom] = cursor
        cursor += chrom_max.get(chrom, 0) + gap
    fig, (ax, qax) = plt.subplots(1, 2, figsize=(15.2, 4.8), gridspec_kw={"width_ratios": [4.3, 1]})
    colors = ["#365F91", "#D9903D"]
    for idx, chrom in enumerate(chroms):
        vals = [(offsets[chrom] + bp, -math.log10(pv)) for c, bp, pv in points if c == chrom]
        if vals:
            x, y = zip(*vals); ax.scatter(x, y, s=2.2, color=colors[idx % 2], alpha=.70, rasterized=True, lw=0)
    if n:
        ax.axhline(-math.log10(threshold), color="#B22222", ls="--", lw=.9, label="Bonferroni 0.05/N")
    centers = [offsets[c] + chrom_max.get(c, 0) / 2 for c in chroms]
    ax.set_xticks(centers, [f"chr{c}" for c in chroms], rotation=45, ha="right", fontsize=7)
    ax.set_ylabel("−log10(P)"); ax.set_xlabel("Chromosome")
    ax.set_title(f"{trait}: zero- and 2.5% tail-excluded SNP-GWAS (all 16 chromosomes)", loc="left", fontsize=11)
    ax.spines[["top", "right"]].set_visible(False); ax.grid(axis="y", color="#E6E6E6", lw=.5)
    q = np.sort(np.asarray(qq, dtype=float))
    if len(q):
        obs = -np.log10(q); exp = -np.log10((np.arange(1, len(q) + 1) - .5) / len(q))
        qax.scatter(exp, obs, s=5, color="#365F91", alpha=.55, lw=0, rasterized=True)
        lim = max(float(exp.max()), float(obs.max())); qax.plot([0, lim], [0, lim], color="#B22222", lw=.8, ls="--")
        chisq = chi2.isf(np.clip(q, np.nextafter(0, 1), 1), 1)
        lambda_gc = float(np.median(chisq) / chi2.ppf(.5, 1))
    else:
        lambda_gc = math.nan
    qax.set_xlabel("Expected −log10(P)"); qax.set_ylabel("Observed −log10(P)")
    qax.set_title(f"QQ; λGC={lambda_gc:.3f}", fontsize=10)
    qax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(outdir / "zero_excluded_manhattan_qq.png", dpi=420)
    fig.savefig(outdir / "zero_excluded_manhattan_qq.pdf")
    plt.close(fig)

    # Separate publication panels in addition to the combined diagnostic.
    mfig, maxis = plt.subplots(figsize=(14.2, 4.8))
    for idx, chrom in enumerate(chroms):
        vals = [(offsets[chrom] + bp, -math.log10(pv)) for c, bp, pv in points if c == chrom]
        if vals:
            xx, yy = zip(*vals); maxis.scatter(xx, yy, s=2.2, color=colors[idx % 2], alpha=.70, rasterized=True, lw=0)
    maxis.axhline(-math.log10(threshold), color="#B22222", ls="--", lw=.9)
    maxis.set_xticks(centers, [f"chr{c}" for c in chroms], rotation=45, ha="right", fontsize=7)
    maxis.set_ylabel("−log10(P)"); maxis.set_xlabel("Chromosome")
    maxis.set_title(f"{trait}: zero- and 2.5% tail-excluded SNP-GWAS", loc="left", fontsize=11)
    maxis.spines[["top", "right"]].set_visible(False); maxis.grid(axis="y", color="#E6E6E6", lw=.5)
    mfig.tight_layout(); mfig.savefig(outdir / "manhattan.png", dpi=420); mfig.savefig(outdir / "manhattan.pdf"); plt.close(mfig)

    qfig, qaxis = plt.subplots(figsize=(5.4, 5.0))
    if len(q):
        qaxis.scatter(exp, obs, s=6, color="#365F91", alpha=.55, lw=0, rasterized=True)
        qaxis.plot([0, lim], [0, lim], color="#B22222", lw=.8, ls="--")
    qaxis.set_xlabel("Expected −log10(P)"); qaxis.set_ylabel("Observed −log10(P)")
    qaxis.set_title(f"{trait}: SNP-GWAS QQ; λGC={lambda_gc:.3f}", fontsize=10)
    qaxis.spines[["top", "right"]].set_visible(False); qfig.tight_layout()
    qfig.savefig(outdir / "qq.png", dpi=420); qfig.savefig(outdir / "qq.pdf"); plt.close(qfig)
    pheno = pd_read_summary(category, trait)
    with (outdir / "summary_row.tsv").open("w") as fh:
        fh.write("category\ttrait\tn_positive\tn_chr\tn_tests\ttop_marker\ttop_p\tbonferroni\tlambda_gc\tstatus\n")
        fh.write(f"{category}\t{trait}\t{pheno}\t{complete_chr}\t{n}\t{top_marker}\t{minp}\t{threshold}\t{lambda_gc}\t{status}\n")
    print(f"trait={trait} n={n} n_chr={complete_chr} top={top_marker} p={minp} lambda={lambda_gc} status={status}")


def pd_read_summary(category, trait):
    with (RUN / "manifests/phenotype_filter_audit.tsv").open() as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["category"] == category and r["trait"] == trait:
                return int(r["n_retained"])
    return 0


if __name__ == "__main__":
    main(int(sys.argv[1]))
