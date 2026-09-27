#!/usr/bin/env python3
"""Extract SyRI large-rearrangement records (INV / DUP / TRA) as an independent,
SyRI single-evidence SV set that runs in PARALLEL to the Tier-1 consensus set.

Rationale: read-based CuteSV is poor at large inversions, duplications and
translocations, so the Tier-1 "CuteSV INTERSECT assembly" intersection collapses
these three classes to almost nothing (population: INV 268 / DUP 102 / TRA 122
clusters). SyRI, being synteny-based, is the appropriate caller for exactly these
large rearrangements. We therefore report INV/DUP/TRA from SyRI as a separate
parallel layer, so the manuscript can use:
  * DEL / INS  -> Tier-1 read+assembly consensus (high-confidence)
  * INV/DUP/TRA -> this SyRI single-evidence set (more complete large rearr.)

This does NOT modify Tier-1 (scripts/30/31) nor the Tier-2/Tier-3 layers
(scripts/40/41); INVDP/INVTR/CPG/CPL stay in Tier-2, HDR in Tier-3.

Input: per-sample normalized SyRI VCFs under results/norm/*.syri.norm.vcf
(coordinates already in Africa_hap2 reference space, chr01B..chr16B). Each
sample has both a .vcf and a .vcf.gz; we keep one per sample (prefer .vcf.gz) to
avoid the well-known double-counting trap.

Output columns match scripts/40 (and what scripts/41 consumes):
  Sample Chrom Start End SVTYPE SVLEN_bp SVLEN_raw SUPP Caller_Combo Size_Bin
SUPP is fixed to 1 and Caller_Combo to "syri" because these are single-caller
calls. INV/DUP/TRA all carry a clean reference span (POS..END) for clustering.
"""
import argparse
import glob
import gzip
import os
import re
import sys

SVTYPE_RE = re.compile(r"(?:^|;)SVTYPE=([^;]+)(?:;|$)")
SVLEN_RE = re.compile(r"(?:^|;)SVLEN=(-?\d+)(?:;|$)")
END_RE = re.compile(r"(?:^|;)END=([0-9]+)(?:;|$)")

# The three standard large rearrangement classes (user-requested scope).
# Complex inverted classes INVDP/INVTR and CPG/CPL remain in Tier-2.
KEEP_TYPES = ["INV", "DUP", "TRA"]

HEADER = [
    "Sample", "Chrom", "Start", "End", "SVTYPE",
    "SVLEN_bp", "SVLEN_raw", "SUPP", "Caller_Combo", "Size_Bin",
]


def size_bin(svlen_bp):
    if svlen_bp is None:
        return "NA"
    v = abs(svlen_bp)
    if v < 1000:
        return "50bp-1kb"
    if v < 10000:
        return "1-10kb"
    if v < 100000:
        return "10-100kb"
    return ">100kb"


def first_int(rx, s):
    m = rx.search(s)
    return int(m.group(1)) if m else None


def open_text(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8")
    return open(path, encoding="utf-8")


def resolve_sample_files(norm_dir):
    """Map sample -> one SyRI norm VCF path, preferring .vcf.gz over .vcf."""
    files = {}
    for path in sorted(glob.glob(os.path.join(norm_dir, "*.syri.norm.vcf*"))):
        if path.endswith(".tbi"):
            continue
        base = os.path.basename(path)
        sample = re.sub(r"\.syri\.norm\.vcf(\.gz)?$", "", base)
        if sample not in files or path.endswith(".gz"):
            files[sample] = path
    return files


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--norm-dir", required=True,
                    help="Directory with *.syri.norm.vcf[.gz].")
    ap.add_argument("--out-dir", required=True,
                    help="Output directory for the records TSV.")
    ap.add_argument("--types", default=",".join(KEEP_TYPES),
                    help="Comma-separated SVTYPEs to extract from the SyRI VCFs "
                         "(default INV,DUP,TRA). E.g. 'DEL,INS' for the SyRI indel comparison.")
    ap.add_argument("--out-name", default="syri_large_rearrangements.records.tsv",
                    help="Output TSV filename within --out-dir.")
    ap.add_argument("--min-svlen", type=int, default=50,
                    help="Minimum |SVLEN| in bp to keep for length-bearing types (default 50).")
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()

    type_order = [t.strip().upper() for t in args.types.split(",") if t.strip()]

    os.makedirs(args.out_dir, exist_ok=True)
    out_all = os.path.join(args.out_dir, args.out_name)
    if os.path.exists(out_all) and not args.force:
        sys.exit(f"[extract] refuse to overwrite {out_all} (use --force)")

    sample_files = resolve_sample_files(args.norm_dir)
    if not sample_files:
        sys.exit(f"[extract] no *.syri.norm.vcf[.gz] under {args.norm_dir}")

    keep_types = set(type_order)
    totals = {t: 0 for t in keep_types}

    tmp_all = out_all + ".tmp"
    with open(tmp_all, "w", encoding="utf-8") as fh_all:
        fh_all.write("\t".join(HEADER) + "\n")
        for sample in sorted(sample_files):
            vcf = sample_files[sample]
            counts = {t: 0 for t in keep_types}
            with open_text(vcf) as h:
                for line in h:
                    if line.startswith("#"):
                        continue
                    cols = line.rstrip("\n").split("\t")
                    if len(cols) < 8:
                        continue
                    info = cols[7]
                    svtype_m = SVTYPE_RE.search(info)
                    if not svtype_m:
                        continue
                    svtype = svtype_m.group(1).strip().upper()
                    if svtype not in keep_types:
                        continue
                    start = int(cols[1])
                    svlen_raw = first_int(SVLEN_RE, info)
                    svlen_bp = abs(svlen_raw) if svlen_raw is not None else None
                    end = first_int(END_RE, info)
                    if end is None:
                        end = start + (svlen_bp - 1) if svlen_bp else start
                    # length filter on length-bearing types (INV/DUP); TRA may
                    # keep records even with small spans but in practice carries
                    # a clean SVLEN; the >=min filter mirrors normalization.
                    if svlen_bp is not None and svlen_bp < args.min_svlen:
                        continue
                    row = [
                        sample, cols[0], str(start), str(end), svtype,
                        str(svlen_bp) if svlen_bp is not None else "NA",
                        str(svlen_raw) if svlen_raw is not None else "NA",
                        "1", "syri", size_bin(svlen_bp),
                    ]
                    fh_all.write("\t".join(row) + "\n")
                    counts[svtype] += 1
                    totals[svtype] += 1
            print(f"[extract] {sample}: "
                  + " ".join(f"{t}={counts[t]}" for t in type_order), file=sys.stderr)

    os.replace(tmp_all, out_all)

    grand = sum(totals.values())
    print(f"[extract] samples={len(sample_files)}")
    print(f"[extract] {'/'.join(type_order)} records = {grand}")
    for t in type_order:
        print(f"           {t}: {totals[t]}")
    print(f"[extract] wrote: {out_all}")


if __name__ == "__main__":
    main()
