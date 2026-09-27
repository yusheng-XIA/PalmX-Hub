#!/usr/bin/env python3
"""Re-derive favourable loci from the reported SV-GWAS (60 retained traits, zero/extreme-trimmed phenotypes,
EMMAX, Bonferroni 0.05/370,706) using the author rules (92_single_trait_sv_validation.py::choose_genotype,
build_ideal_haplotype_sv_catalog.py genotype-consistent core set, 09_mosaic_step9.py favour/allele rules)."""
import subprocess, sys, os
from collections import defaultdict, Counter
import numpy as np, pandas as pd
B = "${ANALYSIS_DIR}"
RUN = f"{B}/22_answer_reviews/00_ms/03_V3/05_figure/07_GWAS_zero_extreme_excluded_20260723"
SVVCF = f"{B}/14_pan_genome/06_Minigraph/Pangenie/02_sv_combined/步骤四_过滤分类统计/sv.qc.id.vcf.gz"
BUB = f"{B}/14_pan_genome/06_Minigraph/final_outputs/pangenome.clip.bub.vcf.gz"
TABIX = "${DATA_DIR2}/anaconda3/bin/tabix"
OUT = "${CLUSTER_WORK}/fav_rederive/out"; os.makedirs(OUT, exist_ok=True)
MIN_N = 5; BONF = 0.05 / 370706
TRAIT_DIRECTION = {"Kernel_oil_content": 1, "Kernel_length_mm": 1, "Kernel_width_mm": 1, "Kernel_weight_g": 1,
    "Flesh_thickness_mm": 1, "Cluster_weight_kg": 1, "Total_fruit_weight": 1, "Good_fruit_weight_kg": 1,
    "Poor_fruit_weight_kg": -1, "Poor_fruit_num": -1, "Shell_weight_g": -1, "Shell_thickness_mm": -1,
    "Kernel_moisture_content": -1, "Flesh_moisture_content": -1, "SOD_U_g": 1}
EXCLUDE = {"American_hap1"}; NAME_MAP = {"dura": "EG_dura", "pisifera": "EG_pisifera"}

pairs = pd.read_csv(f"{RUN}/tables/all_genomewide_significant_variant_trait_pairs.tsv", sep="\t")
sv = pairs[pairs.modality == "SV"].copy()
print("SV trait pairs:", len(sv), "unique SVs:", sv.variant_id.nunique(), "traits:", sv.trait.value_counts().to_dict())
assert (sv.p <= BONF).all()
sv.to_csv(f"{OUT}/reported_sv_pairs.tsv", sep="\t", index=False)

# genotypes + REF/ALT lengths from the SV VCF used for SV-GWAS
hdr = subprocess.run([TABIX, "-H", SVVCF], capture_output=True, text=True).stdout.splitlines()
samples = hdr[-1].split("\t")[9:]
geno, meta = {}, {}
ids = set(sv.variant_id)
for (c, p), _ in sv.groupby(["chrom", "pos"]):
    out = subprocess.run([TABIX, SVVCF, f"{c}:{p}-{p}"], capture_output=True, text=True).stdout
    for ln in out.splitlines():
        f = ln.split("\t")
        if f[2] not in ids: continue
        d = []
        for g in f[9:]:
            g = g.split(":")[0].replace("|", "/")
            if "." in g: d.append(np.nan); continue
            d.append(float(sum(int(x) != 0 for x in g.split("/"))))
        geno[f[2]] = np.array(d); meta[f[2]] = (f[0], int(f[1]), len(f[3]), len(f[4].split(",")[0]), f[4].count(",") + 1)
miss = ids - set(geno); print("SVs without genotype record:", len(miss))

def label(x):
    return "NA" if np.isnan(x) else ("0/0" if x <= .5 else "0/1" if x <= 1.5 else "1/1")

def pheno(trait, cat):
    rows = {}
    for ln in open(f"{RUN}/phenotypes/{cat}/{trait}.txt", encoding="utf-8-sig", errors="replace"):
        p = ln.strip().split()
        if len(p) < 3 or p[2] in {"", "NA", "nan", "NaN"}: continue
        try: rows[p[1]] = float(p[2])
        except ValueError: pass
    return pd.Series(rows)

rec_rows = []
for r in sv.itertuples():
    if r.variant_id not in geno: continue
    y = pheno(r.trait, r.category)
    df = pd.DataFrame({"IID": samples, "g": [label(x) for x in geno[r.variant_id]]}).set_index("IID").join(y.rename("t"), how="inner").dropna()
    df = df[df.g != "NA"]
    means, ns = {}, {}
    for gl in ("0/0", "0/1", "1/1"):
        v = df.t[df.g == gl]
        ns[gl] = len(v)
        if len(v) >= MIN_N: means[gl] = v.mean()
    if len(means) < 2: rec = "-"
    else:
        hi = max(means, key=means.get); lo = min(means, key=means.get)
        rec = (hi if TRAIT_DIRECTION[r.trait] > 0 else lo) if r.trait in TRAIT_DIRECTION else "nondirectional"
    rec_rows.append(dict(sv_id=r.variant_id, trait=r.trait, directional=r.trait in TRAIT_DIRECTION, recommended=rec,
                         n00=ns["0/0"], n01=ns["0/1"], n11=ns["1/1"], p=r.p))
R = pd.DataFrame(rec_rows); R.to_csv(f"{OUT}/single_trait_recommendation.tsv", sep="\t", index=False)
core = []
for s, g in R[R.directional].groupby("sv_id"):
    gs = sorted({x for x in g.recommended if x in {"0/0", "0/1", "1/1"}})
    cls = "consistent" if len(gs) == 1 else ("tradeoff" if len(gs) > 1 else "no_recommendation")
    c, p, rl, al, nalt = meta[s]
    core.append(dict(sv_id=s, chrom=c, pos=p, ref_len=rl, alt_len=al, n_alt=nalt, n_directional=len(g),
                     traits=";".join(sorted(g.trait)), cls=cls, rec=gs[0] if len(gs) == 1 else ";".join(gs)))
C = pd.DataFrame(core); C.to_csv(f"{OUT}/core_catalog.tsv", sep="\t", index=False)
print(C.cls.value_counts().to_dict()); print(C[C.cls == "consistent"].rec.value_counts().to_dict())
fav = C[(C.cls == "consistent") & C.rec.isin(["0/0", "1/1"]) & (C.ref_len != C.alt_len)].copy()
fav["Target"] = np.where(fav.rec == "0/0", "REF", "ALT")
print("favourable loci:", len(fav), fav.chrom.value_counts().to_dict())
# donor alleles from the bubble VCF (09_mosaic_step9.py::bub_alleles)
bh = subprocess.run([TABIX, "-H", BUB], capture_output=True, text=True).stdout.splitlines()
bs = bh[-1].split("\t")[9:]
al_rows = []
for r in fav.itertuples():
    out = subprocess.run([TABIX, BUB, f"{r.chrom}:{r.pos}-{r.pos}"], capture_output=True, text=True).stdout
    rec = None
    for ln in out.splitlines():
        f = ln.split("\t")
        if int(f[1]) == r.pos: rec = f; break
    if rec is None: continue
    alen = [len(rec[3])] + [len(a) for a in rec[4].split(",")]
    for name, g in zip(bs, rec[9:]):
        if name in EXCLUDE: continue
        a = g.split(":")[0].replace("|", "/").split("/")[0]
        if a in (".", ""): continue
        gi = int(a); L = alen[gi] if gi < len(alen) else -1
        al_rows.append(dict(sv_id=r.sv_id, donor=NAME_MAP.get(name, name), carries_ALT=int(L == r.alt_len)))
A = pd.DataFrame(al_rows); A.to_csv(f"{OUT}/donor_alleles.tsv", sep="\t", index=False)
fav[["sv_id", "chrom", "pos", "Target", "rec", "ref_len", "alt_len", "traits"]].rename(
    columns={"sv_id": "SV", "chrom": "Chrom", "pos": "Pos"}).to_csv(f"{OUT}/fav_loci_reported.tsv", sep="\t", index=False)
print("SVs with bubble record:", A.sv_id.nunique(), "of", len(fav))
