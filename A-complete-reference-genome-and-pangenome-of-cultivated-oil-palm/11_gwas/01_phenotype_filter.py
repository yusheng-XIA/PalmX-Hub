#!/usr/bin/env python3
"""Phenotype filtering before GWAS.

For each of the 64 phenotypes (FID IID value, one file per trait under <in_dir>/<category>/):
exact zero values are excluded; among non-zero observations, values below the 2.5th or above the 97.5th
percentile of the non-zero distribution are excluded; traits with < 50 observations after filtering are
not analysed (60 traits retained). Excluded values are written as NA so that the sample order is unchanged.

usage: 01_phenotype_filter.py <in_dir> <out_dir>
"""
import csv
import sys
from pathlib import Path

import numpy as np

LOW_Q, HIGH_Q, MIN_N = 0.025, 0.975, 50


def value(x):
    try:
        return None if x in {"NA", "", "-9"} else float(x)
    except ValueError:
        return None


def main(in_dir, out_dir):
    in_dir, out_dir = Path(in_dir), Path(out_dir)
    audit, tasks = [], []
    for src in sorted(in_dir.glob("*/*.txt")):
        category, trait = src.parent.name, src.stem
        raw = [(f[0], f[1], value(f[2])) for f in (l.split() for l in open(src)) if len(f) >= 3]
        nz = np.array([v for *_, v in raw if v is not None and v != 0], float)
        lo, hi = (np.quantile(nz, LOW_Q), np.quantile(nz, HIGH_Q)) if len(nz) else (np.nan, np.nan)
        kept = [(a, b, v if v is not None and v != 0 and lo <= v <= hi else None) for a, b, v in raw]
        n_keep = sum(v is not None for *_, v in kept)
        status = "eligible" if n_keep >= MIN_N else "skip_lt50"
        audit.append([category, trait, len(raw), sum(v == 0 for *_, v in raw), len(nz), lo, hi, n_keep, status])
        if status == "eligible":
            dst = out_dir / "phenotypes" / category / src.name
            dst.parent.mkdir(parents=True, exist_ok=True)
            with open(dst, "w") as fh:
                for a, b, v in kept:
                    fh.write(f"{a}\t{b}\t{'NA' if v is None else format(v, '.15g')}\n")
            tasks.append([len(tasks) + 1, category, trait, dst])
    (out_dir / "manifests").mkdir(parents=True, exist_ok=True)
    with open(out_dir / "manifests/phenotype_filter_audit.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["category", "trait", "n_total", "n_zero", "n_nonzero", "q025", "q975", "n_retained", "status"])
        w.writerows(audit)
    with open(out_dir / "manifests/sv_tasks.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["task_id", "category", "trait", "phenotype"])
        w.writerows(tasks)
    print(f"traits={len(audit)} retained={len(tasks)}")


if __name__ == "__main__":
    main(*sys.argv[1:3])
