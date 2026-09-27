#!/usr/bin/env python3
"""Association of the known SHELL mutations with fruit-form traits.

shMPOB : chr01B:3,259,207 (L29P in FL-Hap2 SHELL evm.model.chr01B.166; L28P in the original numbering)
shAVROS: chr01B:3,259,200 (K31N; K30N in the original numbering). The FL-Hap2 reference carries shAVROS (N31).
The number of mutant alleles per accession (shMPOB + shAVROS) is tested under the SNP-model EMMAX settings
(SNP kinship, intercept), alone and with the dosage of the lead SNP of the trait as an additional covariate.

usage: 07_shell_mutation_test.py --mutant chr01B:3259207=<base> --mutant chr01B:3259200=<base> \
           --trait Nut_weight_g --lead-snp chr01B:<pos>
"""
import argparse

import numpy as np
import pandas as pd

from config import BED, BIM, FAM, KIN_SNP, RUN
from emmaxpy import GLS, read_bed_rows, read_matrix, read_pheno


def dosage(bim, ids, site, allele=None):
    """Count of `allele` (default: bim A1) at chrom:pos."""
    chrom, pos = site.split(":")
    hit = bim[(bim.chr == chrom) & (bim.pos == int(pos))]
    if hit.empty:
        raise SystemExit(f"{site} not in the SNP set")
    i = int(hit.index[0])
    d = read_bed_rows(BED, len(ids), i, i + 1)[0].astype(float)          # count of A1
    if allele is not None:
        a1, a2 = hit.iloc[0][["a1", "a2"]]
        if allele == a2:
            d = 2 - d
        elif allele != a1:
            raise SystemExit(f"{site}: allele {allele} not in {a1}/{a2}")
    return d


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mutant", action="append", required=True, help="chrom:pos=mutant_base")
    ap.add_argument("--trait", required=True)
    ap.add_argument("--lead-snp")
    a = ap.parse_args()
    ids = [l.split()[1] for l in open(FAM)]
    bim = pd.read_csv(BIM, sep=r"\s+", header=None, names=["chr", "id", "cm", "pos", "a1", "a2"], dtype={"chr": str})
    tasks = pd.read_csv(RUN / "manifests/sv_tasks.tsv", sep="\t")
    cat = dict(zip(tasks.trait, tasks.category))[a.trait]
    y = read_pheno(RUN / f"phenotypes/{cat}/{a.trait}.txt", ids)
    k = np.isfinite(y)
    K = read_matrix(KIN_SNP)
    mut = sum(dosage(bim, ids, s.split("=")[0], s.split("=")[1]) for s in a.mutant)
    X = np.ones((len(ids), 1))
    beta, se, p = GLS(y[k], X[k], K[np.ix_(k, k)]).scan(mut[k][None, :])
    print(f"{a.trait}\tmutant_alleles\tbeta={beta[0]:.4g}\tse={se[0]:.4g}\tP={p[0]:.3g}")
    if a.lead_snp:
        lead = dosage(bim, ids, a.lead_snp)
        lead = np.where(np.isnan(lead), np.nanmean(lead[k]), lead)
        Xc = np.c_[X, lead]
        beta, se, p = GLS(y[k], Xc[k], K[np.ix_(k, k)]).scan(mut[k][None, :])
        print(f"{a.trait}\tmutant_alleles|lead_snp\tbeta={beta[0]:.4g}\tse={se[0]:.4g}\tP={p[0]:.3g}")


if __name__ == "__main__":
    main()
