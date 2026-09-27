#!/usr/bin/env python3
"""Genome-wide null for identical lipid-enzyme copy number in Cocos nucifera and three Elaeis genomes
(Supplementary Fig. 2a,g).

Focal genomes: Cocos_nucifera, American_hap1 (FL-Hap1), Dura (TK) and Pisifera (NS).
Background: all orthogroups of the 30-genome OrthoFinder run with >= 1 gene locus in each focal genome
(C. nucifera protein isoforms collapsed to gene loci). In each of 10,000 permutations every enzyme class is
replaced by an orthogroup drawn at random with the same C. nucifera copy number (1, 2, 3, 4, 5-6, 7-10, >10);
one-sided P = proportion of permutations with at least as many identical classes as observed (+1 pseudocount).

usage: 05_dosage_conservation_null.py Orthogroups.tsv enzyme_copy_number.tsv OUT_PREFIX
  enzyme_copy_number.tsv: Enzyme, Cocos_nucifera, American_hap1, Dura, Pisifera (curated gene-locus counts)
"""
import re
import sys

import numpy as np
import pandas as pd

FOC = ["Cocos_nucifera", "American_hap1", "Dura", "Pisifera"]
CN_BINS = [0.5, 1.5, 2.5, 3.5, 4.5, 6.5, 10.5, np.inf]
CN_LAB = ["1", "2", "3", "4", "5-6", "7-10", ">10"]
NPERM = 10_000
rng = np.random.default_rng(20260924)


def locus_counts(orthogroups_tsv):
    """Per-orthogroup gene-locus counts; Cocos isoform suffixes (.N) are collapsed."""
    rows = []
    with open(orthogroups_tsv) as f:
        species = f.readline().rstrip("\n").split("\t")[1:]
        for line in f:
            r = line.rstrip("\n").split("\t"); r += [""] * (len(species) + 1 - len(r))
            cnt = []
            for s, c in zip(species, r[1:]):
                genes = {x.split("|", 1)[-1] for x in c.split(", ")} if c else set()
                if s == "Cocos_nucifera":
                    genes = {re.sub(r"\.\d+$", "", g) for g in genes}
                cnt.append(len(genes))
            rows.append([r[0]] + cnt)
    return pd.DataFrame(rows, columns=["Orthogroup"] + species).set_index("Orthogroup")


def main(og_tsv, enzyme_tsv, out):
    og = locus_counts(og_tsv)
    U = og[og[FOC].gt(0).all(axis=1)].copy()
    U["cn_bin"] = pd.cut(U["Cocos_nucifera"], CN_BINS, labels=CN_LAB).astype(str)
    U["identical"] = U[FOC].nunique(axis=1).eq(1)

    enz = pd.read_csv(enzyme_tsv, sep="\t")
    enz["cn_bin"] = pd.cut(enz["Cocos_nucifera"], CN_BINS, labels=CN_LAB).astype(str)
    enz["identical"] = enz[FOC].nunique(axis=1).eq(1)
    obs = int(enz.identical.sum())

    pools = {b: g["identical"].to_numpy() for b, g in U.groupby("cn_bin")}
    null = np.zeros(NPERM, dtype=int)
    for b in enz.cn_bin:
        null += pools[b][rng.integers(0, len(pools[b]), NPERM)]
    p_upper = (1 + (null >= obs).sum()) / (NPERM + 1)
    p_lower = (1 + (null <= obs).sum()) / (NPERM + 1)

    U.groupby("cn_bin").agg(n_OG=("identical", "size"), frac_identical=("identical", "mean")).reindex(CN_LAB) \
        .to_csv(f"{out}.background_by_copy_number.tsv", sep="\t")
    pd.DataFrame({"n_identical": np.arange(len(enz) + 1),
                  "permutations": [int((null == k).sum()) for k in range(len(enz) + 1)]}) \
        .to_csv(f"{out}.null_distribution.tsv", sep="\t", index=False)
    summary = (f"background_orthogroups\t{len(U)}\nbackground_frac_identical\t{U.identical.mean():.4f}\n"
               f"observed_identical\t{obs}/{len(enz)}\nnull_mean\t{null.mean():.3f}\n"
               f"null_95pct\t{int(np.percentile(null, 2.5))}-{int(np.percentile(null, 97.5))}\n"
               f"P_upper\t{p_upper:.4g}\nP_lower\t{p_lower:.4g}\n")
    open(f"{out}.summary.tsv", "w").write(summary)
    print(summary)


if __name__ == "__main__":
    main(*sys.argv[1:4])
