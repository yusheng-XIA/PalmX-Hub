#!/usr/bin/env python3
"""Ka, Ks and Ka/Ks of one-to-one allelic gene pairs between the two haplotypes of a hybrid.

Pairs: FL, one-to-one OrthoFinder orthogroups between FL-Hap2 and FL-Hap1; TN, the 20,952 curated gene pairs.
Coding sequences are translated, aligned at the protein level (Biopython global alignment, BLOSUM62,
gap opening 10, gap extension 0.5) and back-translated to codon alignments. Pairs are discarded when a CDS is
shorter than 100 codons or has ambiguous bases or internal stops, or when amino-acid identity < 50% or the
alignment gap fraction > 50%. Ka/Ks is estimated with KaKs_Calculator (gamma-MYN, -m GMYN).

usage: 06_allelic_kaks.py --pairs pairs.tsv --cds-a hapA_cds.fa --cds-b hapB_cds.fa --prefix FL
  pairs.tsv: columns pair_group, gene_a, gene_b (gene_a on allele A: FL-Hap2 or dura/TK-like haplotype)
"""
import argparse
import re
import subprocess

import pandas as pd
from Bio import SeqIO, pairwise2
from Bio.Align import substitution_matrices
from Bio.Seq import Seq

KAKS_COLUMNS = ["Sequence", "Method", "Ka", "Ks", "Ka/Ks", "P-Value(Fisher)", "Length", "S-Sites", "N-Sites",
                "Fold-Sites", "Substitutions", "S-Substitutions", "N-Substitutions", "Fold-S-Substitutions",
                "Fold-N-Substitutions", "Divergence-Time", "Substitution-Rate-Ratio", "GC", "ML-Score", "AICc",
                "Akaike-Weight", "Model"]


def read_cds(path):
    return {r.id.split()[0]: str(r.seq).upper().replace("U", "T") for r in SeqIO.parse(path, "fasta")}


def clean_cds(seq, min_aa):
    seq = re.sub(r"\s+", "", seq)
    seq = re.sub(r"[^ACGTN]", "N", seq)
    seq = seq[: len(seq) - len(seq) % 3]
    if len(seq) < min_aa * 3:
        return None, None, "short_cds"
    if "N" in seq:
        return None, None, "ambiguous_N"
    aa = str(Seq(seq).translate(table=1))
    if aa.endswith("*"):
        aa, seq = aa[:-1], seq[:-3]
    if "*" in aa:
        return None, None, "internal_stop"
    if len(aa) < min_aa:
        return None, None, "short_protein"
    return seq, aa, "ok"


def codon_align(cds1, aa1, cds2, aa2, matrix):
    if len(aa1) == len(aa2):
        ident = sum(a == b for a, b in zip(aa1, aa2)) / len(aa1)
        return cds1, cds2, ident, 0.0
    aln = pairwise2.align.globalds(aa1, aa2, matrix, -10, -0.5, one_alignment_only=True)[0]
    c1 = [cds1[i:i + 3] for i in range(0, len(cds1), 3)]
    c2 = [cds2[i:i + 3] for i in range(0, len(cds2), 3)]
    i1 = i2 = matches = compared = gaps = 0
    n1, n2 = [], []
    for a, b in zip(aln.seqA, aln.seqB):
        n1.append("---" if a == "-" else c1[i1]); i1 += a != "-"
        n2.append("---" if b == "-" else c2[i2]); i2 += b != "-"
        if a == "-" or b == "-":
            gaps += 1
        else:
            compared += 1; matches += a == b
    return "".join(n1), "".join(n2), matches / compared if compared else 0.0, gaps / len(aln.seqA)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pairs", required=True)
    ap.add_argument("--cds-a", required=True)
    ap.add_argument("--cds-b", required=True)
    ap.add_argument("--prefix", required=True)
    ap.add_argument("--kaks-bin", default="KaKs_Calculator")
    ap.add_argument("--min-aa", type=int, default=100)
    ap.add_argument("--min-identity", type=float, default=0.50)
    ap.add_argument("--max-gap", type=float, default=0.50)
    a = ap.parse_args()

    pairs = pd.read_csv(a.pairs, sep="\t", dtype=str)
    ca, cb = read_cds(a.cds_a), read_cds(a.cds_b)
    matrix = substitution_matrices.load("BLOSUM62")
    qc, axt = [], f"{a.prefix}.codon_alignment.axt"
    with open(axt, "w") as out:
        for r in pairs.itertuples(index=False):
            pid = f"{r.pair_group}__{r.gene_a}__{r.gene_b}"
            sa, sb = ca.get(r.gene_a), cb.get(r.gene_b)
            if sa is None or sb is None:
                qc.append((pid, r.gene_a, r.gene_b, "fail", "missing_cds", None, None)); continue
            na, aa_a, why_a = clean_cds(sa, a.min_aa)
            nb, aa_b, why_b = clean_cds(sb, a.min_aa)
            if na is None or nb is None:
                qc.append((pid, r.gene_a, r.gene_b, "fail", f"cds_qc:{why_a}/{why_b}", None, None)); continue
            if abs(len(aa_a) - len(aa_b)) / max(len(aa_a), len(aa_b)) > a.max_gap:
                qc.append((pid, r.gene_a, r.gene_b, "fail", "length_gap_too_large", None, None)); continue
            x1, x2, ident, gap = codon_align(na, aa_a, nb, aa_b, matrix)
            status, note = "pass", "ok"
            if ident < a.min_identity:
                status, note = "fail", "low_identity"
            elif gap > a.max_gap:
                status, note = "fail", "high_gap_fraction"
            qc.append((pid, r.gene_a, r.gene_b, status, note, ident, gap))
            if status == "pass":
                out.write(f"{pid}\n{x1}\n{x2}\n\n")
    qc = pd.DataFrame(qc, columns=["pair_id", "gene_a", "gene_b", "status", "note", "aa_identity", "aa_gap_fraction"])
    raw = f"{a.prefix}.GMYN.raw.tsv"
    subprocess.run([a.kaks_bin, "-i", axt, "-o", raw, "-m", "GMYN"], check=True)
    with open(raw) as h:
        header = h.readline().startswith("Sequence")
    kaks = pd.read_csv(raw, sep="\t", header=0 if header else None)
    if not header:
        kaks.columns = KAKS_COLUMNS[: len(kaks.columns)]
    kaks.merge(qc, left_on="Sequence", right_on="pair_id", how="left").to_csv(
        f"{a.prefix}.KaKs_all.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
