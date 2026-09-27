"""In-sample minor-allele counts (MAC) for traits with SNP hits: per chromosome, MAC among the phenotyped accessions of
each trait for every marker (uint16 .npy), used to recompute lambda_GC and hit counts at in-sample MAF >= 0.05."""
import sys, numpy as np, pandas as pd
from config import *
from emmaxpy import read_bed_rows, read_pheno
chrom = sys.argv[1]
ids = [l.split()[1] for l in open(FAM)]; N = len(ids)
tasks = pd.read_csv(RUN / "manifests/sv_tasks.tsv", sep="\t")
TR = ["C12_0_Lauric_acid", "C14_0_Myristic_acid", "Flesh_thickness_mm", "C18_3n3_Alpha_linolenic_acid", "Stem_height_cm",
      "Nut_weight_g", "Nut_length_mm", "Shell_thickness_mm", "Shell_weight_g", "Nut_width_mm", "Petiole_width_cm", "Ca"]
keeps = {}
for cat, tr in zip(tasks.category, tasks.trait):
    if tr in TR: keeps[tr] = np.isfinite(read_pheno(RUN / f"phenotypes/{cat}/{tr}.txt", ids))
rng = pd.read_csv(W / "in/chr_ranges.tsv", sep="\t").set_index("chr").loc[chrom]
out = {tr: np.empty(rng.n, dtype=np.uint16) for tr in keeps}; nobs = {tr: np.empty(rng.n, dtype=np.uint16) for tr in keeps}
CH = 200000
for s in range(0, int(rng.n), CH):
    e = min(int(rng.n), s + CH)
    G = read_bed_rows(BED, N, int(rng.start) + s, int(rng.start) + e)
    for tr, k in keeps.items():
        g = G[:, k]; ok = ~np.isnan(g); ac = np.where(ok, g, 0).sum(1); n2 = 2 * ok.sum(1)
        out[tr][s:e] = np.minimum(ac, n2 - ac).astype(np.uint16); nobs[tr][s:e] = n2.astype(np.uint16)
for tr in keeps:
    d = W / "out/mac" / tr; d.mkdir(parents=True, exist_ok=True)
    np.save(d / f"{chrom}.mac.npy", out[tr]); np.save(d / f"{chrom}.n2.npy", nobs[tr])
print("done", chrom)
