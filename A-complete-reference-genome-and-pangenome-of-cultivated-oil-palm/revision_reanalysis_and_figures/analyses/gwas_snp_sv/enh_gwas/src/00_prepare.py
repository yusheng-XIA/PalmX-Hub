#!${DATA_DIR}/miniconda3/bin/python
from __future__ import annotations

import csv
import hashlib
from pathlib import Path
import numpy as np

RUN = Path(__file__).resolve().parents[1]
GWAS = Path("${ANALYSIS_DIR}/05_GWAS/00_analysis/06_GWAS")
PREV = GWAS / "05_zero_excluded_sensitivity_20260721/manifests/snp_tasks.tsv"
CATEGORIES = ["cold_resistant", "growth", "photosynthesis", "quality", "yield"]
SV_ROOT = Path("${ANALYSIS_DIR}/14_pan_genome/06_Minigraph/Pangenie/02_sv_combined")
SV_TPED = SV_ROOT / "步骤七_SV_GWAS/kinship/sv_assoc.tped"
SV_KIN = SV_ROOT / "步骤七_SV_GWAS/kinship/sv.kinf"
SV_COV = SV_ROOT / "步骤五_群体结构/pca/sv_covariates_5PC_with_intercept.txt"
LOW_Q, HIGH_Q, MIN_N = 0.025, 0.975, 50


def digest_head_tail(path: Path, block: int = 1024 * 1024) -> str:
    h = hashlib.sha256(); size = path.stat().st_size
    with path.open("rb") as fh:
        h.update(fh.read(block))
        if size > block:
            fh.seek(max(0, size - block)); h.update(fh.read(block))
    return h.hexdigest()


def parse_value(x: str):
    if x in {"NA", "", "-9"}:
        return None
    try:
        return float(x)
    except ValueError:
        return None


def main() -> None:
    if RUN.is_symlink():
        raise RuntimeError("run directory cannot be a symlink")
    for p in [PREV, SV_TPED, SV_KIN, SV_COV]:
        if not p.is_file() or p.stat().st_size == 0:
            raise FileNotFoundError(p)

    chrom_info = {}
    with PREV.open() as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            chrom_info.setdefault(row["chrom"], {
                "tped_prefix": row["tped_prefix"], "kinship": row["kinship"],
                "expected_rows": int(row["expected_rows"]),
            })
    if len(chrom_info) != 16:
        raise RuntimeError(f"expected 16 chromosome inputs, found {len(chrom_info)}")

    audit, eligible, identity = [], [], []
    for category in CATEGORIES:
        for source in sorted((GWAS / category / "01_phenotype").glob("*.txt")):
            if "with_header" in source.name or ".backup" in source.name:
                continue
            raw = []
            with source.open(encoding="utf-8", errors="replace") as fh:
                for lineno, line in enumerate(fh, 1):
                    f = line.split()
                    if len(f) < 3:
                        raise RuntimeError(f"invalid phenotype line {source}:{lineno}")
                    raw.append((f[0], f[1], parse_value(f[2])))
            positive = np.asarray([v for _, _, v in raw if v is not None and v != 0], dtype=float)
            lo = float(np.quantile(positive, LOW_Q, method="linear")) if len(positive) else np.nan
            hi = float(np.quantile(positive, HIGH_Q, method="linear")) if len(positive) else np.nan
            n_missing = sum(v is None for _, _, v in raw)
            n_zero = sum(v == 0 for _, _, v in raw if v is not None)
            n_low = sum(v is not None and v != 0 and v < lo for _, _, v in raw)
            n_high = sum(v is not None and v != 0 and v > hi for _, _, v in raw)
            retained = [(fid, iid, v if v is not None and v != 0 and lo <= v <= hi else None) for fid, iid, v in raw]
            n_keep = sum(v is not None for _, _, v in retained)
            outdir = RUN / "phenotypes" / category; outdir.mkdir(parents=True, exist_ok=True)
            target = outdir / source.name
            with target.open("w") as fh:
                for fid, iid, v in retained:
                    fh.write(f"{fid}\t{iid}\t{'NA' if v is None else format(v, '.15g')}\n")
            if [(a, b) for a, b, _ in raw] != [(a, b) for a, b, _ in retained]:
                raise RuntimeError(f"sample order mismatch: {source.stem}")
            status = "eligible" if n_keep >= MIN_N else "skip_lt50_after_zero_tail_filter"
            audit.append([category, source.stem, source, target, len(raw), n_missing, n_zero,
                          len(positive), lo, hi, n_low, n_high, n_keep, status])
            identity.append([source, source.stat().st_size, source.stat().st_mtime_ns, digest_head_tail(source)])
            if status == "eligible":
                eligible.append((category, source.stem, target))

    manifests = RUN / "manifests"; manifests.mkdir(parents=True, exist_ok=True)
    with (manifests / "phenotype_filter_audit.tsv").open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["category", "trait", "source", "filtered", "n_total", "n_missing_original",
                    "n_zero_excluded", "n_nonzero_before_tail_filter", "q025", "q975",
                    "n_low_tail_excluded", "n_high_tail_excluded", "n_retained", "status"])
        w.writerows(audit)

    snp_tasks = []; task_id = 1
    for category, trait, pheno in eligible:
        for chrom in sorted(chrom_info):
            info = chrom_info[chrom]
            out = RUN / "snp/results" / category / trait / f"emmax_chr{chrom}"
            snp_tasks.append([task_id, category, trait, chrom, pheno, info["tped_prefix"],
                              info["kinship"], out, info["expected_rows"]]); task_id += 1
    with (manifests / "snp_tasks.tsv").open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["task_id", "category", "trait", "chrom", "phenotype", "tped_prefix", "kinship", "out_prefix", "expected_rows"])
        w.writerows(snp_tasks)
    with (manifests / "sv_tasks.tsv").open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["task_id", "category", "trait", "phenotype", "out_prefix"])
        for i, (category, trait, pheno) in enumerate(eligible, 1):
            w.writerow([i, category, trait, pheno, RUN / "sv/results" / category / trait / trait])
    with (RUN / "provenance/input_identity.tsv").open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["path", "size_bytes", "mtime_ns", "sha256_head_tail"]); w.writerows(identity)
        for p in [PREV, SV_TPED, SV_KIN, SV_COV]:
            w.writerow([p, p.stat().st_size, p.stat().st_mtime_ns, digest_head_tail(p)])
    print(f"traits_total={len(audit)} eligible={len(eligible)} snp_tasks={len(snp_tasks)} sv_tasks={len(eligible)}")


if __name__ == "__main__":
    main()
