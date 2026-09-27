#!/usr/bin/env python3
"""Parallel, deterministic equivalent of legacy 30_extract_highconf_records.py."""
import argparse
import concurrent.futures as cf
import glob
import os
import re
import shutil
from pathlib import Path

SUPP_RE = re.compile(r"(?:^|;)SUPP=([0-9]+)(?:;|$)")
SUPP_VEC_RE = re.compile(r"(?:^|;)SUPP_VEC=([01]+)(?:;|$)")
SVTYPE_RE = re.compile(r"(?:^|;)SVTYPE=([^;]+)(?:;|$)")
SVLEN_RE = re.compile(r"(?:^|;)SVLEN=(-?\d+)(?:;|$)")
END_RE = re.compile(r"(?:^|;)END=([0-9]+)(?:;|$)")
ORDER_DEFAULT = ["syri", "svimasm", "cutesv"]
ASSEMBLY = {"syri", "svimasm"}
READ = {"cutesv"}
HEADER = ["Sample", "Chrom", "Start", "End", "SVTYPE", "SVLEN_bp",
          "SVLEN_raw", "SUPP", "Caller_Combo", "Size_Bin"]


def size_bin(value):
    if value is None:
        return "NA"
    value = abs(value)
    if value < 1000:
        return "50bp-1kb"
    if value < 10000:
        return "1-10kb"
    if value < 100000:
        return "10-100kb"
    return ">100kb"


def first_int(regex, text):
    match = regex.search(text)
    return int(match.group(1)) if match else None


def process_one(task):
    vcf, merged_dir, part_dir = task
    sample = Path(vcf).name.replace(".3caller.vcf", "")
    order_file = Path(merged_dir) / f"{sample}.callers_order.txt"
    if order_file.exists():
        order = [line.strip() for line in order_file.open() if line.strip()]
    else:
        order = list(ORDER_DEFAULT)
    if not order:
        order = list(ORDER_DEFAULT)
    part = Path(part_dir) / f"{sample}.tsv"
    written = 0
    with open(vcf, encoding="utf-8") as source, part.open("w", encoding="utf-8") as out:
        for line in source:
            if line.startswith("#"):
                continue
            columns = line.rstrip("\n").split("\t")
            if len(columns) < 8:
                continue
            info = columns[7]
            vec_match = SUPP_VEC_RE.search(info)
            if not vec_match:
                continue
            vector = vec_match.group(1)
            callers = [order[i] for i, bit in enumerate(vector)
                       if i < len(order) and bit == "1"]
            caller_set = set(callers)
            if not (caller_set & READ) or not (caller_set & ASSEMBLY):
                continue
            type_match = SVTYPE_RE.search(info)
            svtype = type_match.group(1).strip().upper() if type_match else "NA"
            if svtype in {"BND", "TRANS"}:
                svtype = "TRA"
            start = int(columns[1])
            svlen_raw = first_int(SVLEN_RE, info)
            svlen_bp = abs(svlen_raw) if svlen_raw is not None else None
            end = first_int(END_RE, info)
            if end is None:
                end = start + (svlen_bp - 1) if (svlen_bp and svtype in {"DEL", "DUP", "INV"}) else start
            supp = first_int(SUPP_RE, info) or len(callers)
            combo = "+".join(caller for caller in ORDER_DEFAULT if caller in caller_set)
            row = [sample, columns[0], str(start), str(end), svtype,
                   str(svlen_bp) if svlen_bp is not None else "NA",
                   str(svlen_raw) if svlen_raw is not None else "NA",
                   str(supp), combo, size_bin(svlen_bp)]
            out.write("\t".join(row) + "\n")
            written += 1
    return sample, str(part), written


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--merged-dir", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--workers", type=int, default=16)
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()
    out_path = Path(args.out)
    if out_path.exists() and not args.force:
        raise SystemExit(f"refuse to overwrite {out_path}; use --force")
    out_path.parent.mkdir(parents=True, exist_ok=True)
    vcfs = sorted(glob.glob(os.path.join(args.merged_dir, "*.3caller.vcf")))
    if not vcfs:
        raise SystemExit("no merged VCFs found")
    part_dir = out_path.parent / "highconf_parallel_parts"
    part_dir.mkdir(parents=True, exist_ok=True)
    tasks = [(vcf, args.merged_dir, str(part_dir)) for vcf in vcfs]
    results = {}
    with cf.ProcessPoolExecutor(max_workers=args.workers) as pool:
        for sample, part, count in pool.map(process_one, tasks):
            results[sample] = (part, count)
            print(f"[extract-parallel] {sample}: {count}", flush=True)
    tmp = Path(str(out_path) + ".tmp")
    total = 0
    with tmp.open("w", encoding="utf-8") as out:
        out.write("\t".join(HEADER) + "\n")
        for vcf in vcfs:
            sample = Path(vcf).name.replace(".3caller.vcf", "")
            part, count = results[sample]
            with open(part, encoding="utf-8") as source:
                shutil.copyfileobj(source, out)
            total += count
    os.replace(tmp, out_path)
    print(f"[extract-parallel] wrote {out_path} total_records={total} samples={len(vcfs)}")


if __name__ == "__main__":
    main()
