#!/usr/bin/env python3
"""Part 1: genome-wide vs trait-module ASE direction in FL and TN (gene-first, same metric as Fig. 3k)."""
import sys, numpy as np, pandas as pd
from pathlib import Path
from scipy import stats
sys.path.insert(0, str(Path(__file__).parent)); import tha
A = Path("${ANALYSIS_DIR}/22_answer_reviews")
FAT = A/"10_fatty_acid_pathway_reconstruction/tables/fa_pathway_gene_catalog_integrated.tsv"
EM = A/"09_translational_suppression_reexamine/7_转录组文件/04_GO/Africa_hap2/Africa_hap2.emapper.annotations"
F3 = A/"00_ms/03_V3/03_figure3"
BK = F3/"01_ASE/01_bk/ase/config/validated_one_to_one_gene_pairs.tsv"
SL = F3/"01_ASE/02_seedless/ase/config/validated_one_to_one_gene_pairs.tsv"
UNI = F3/"05_multiomics_integration/runs/RUN-MULTIOMICS-ALLELE-CNS-V4-20260723-001/outputs/stage3_existing_ASE_unification_attempt001/gene_stage_ASE_unified.tsv.gz"
FILT = A/"00_ms/05_MS/new_revision/RUN-R2-012-DIAGCONC-20260906-001/results/filtered_gene_stage_ASE.tsv.gz"
OUT = Path("${CLUSTER_WORK}/enh_A/out"); OUT.mkdir(exist_ok=True, parents=True)
SRC = sys.argv[1] if len(sys.argv) > 1 else "unified"
rng = np.random.default_rng(20260924); NPERM = 10000
WIN = ["Days 0–65", "Days 80–140", "Days 155–185", "Hours 12–72", "All stages"]

es = OUT/"empty_shell.tsv"; pd.DataFrame(columns=["Genome","Gene_ID","Identity","Coverage"]).to_csv(es, sep="\t", index=False)
cat = tha.build_catalog(FAT, EM, es)
bk = pd.read_csv(BK, sep="\t").rename(columns={"gene_a":"gene_dura"}); sl = pd.read_csv(SL, sep="\t").rename(columns={"gene_a":"gene_africa"})
bridge = bk[["orthogroup","gene_dura"]].merge(sl[["orthogroup","gene_africa"]], on="orthogroup")
mem_fl = cat[["gene_africa","trait_module","family","preferred_name"]].rename(columns={"gene_africa":"gene_id"}).assign(analysis="FL")
mem_tn = cat.merge(bridge, on="gene_africa")[["gene_dura","trait_module","family","preferred_name"]].rename(columns={"gene_dura":"gene_id"}).assign(analysis="TN")
mem = pd.concat([mem_fl, mem_tn]); mem["trait_module"] = mem.trait_module.astype(str)
mem.to_csv(OUT/"module_membership.tsv", sep="\t", index=False)

use = ["analysis","gene_id","stage","stage_group","eligible","robust_ase","ase_call","log2_allele_ratio","informative_fragments","ref_fragments","alt_fragments"]
ase = pd.read_csv(UNI if SRC == "unified" else FILT, sep="\t", usecols=use, low_memory=False)
for c in ["eligible","robust_ase"]: ase[c] = ase[c].astype(str).str.lower().eq("true")
e = ase[ase.eligible].copy()
e2 = e.copy(); e2["stage_group"] = "All stages"; e = pd.concat([e, e2])

# gene-level table per analysis x window
g = e.groupby(["analysis","stage_group","gene_id"]).agg(
        n_elig=("stage","size"), robust_any=("robust_ase","any"),
        med_all=("log2_allele_ratio","median"), mean_depth=("informative_fragments","mean"),
        refsum=("ref_fragments","sum"), altsum=("alt_fragments","sum")).reset_index()
rb = e[e.robust_ase].groupby(["analysis","stage_group","gene_id"]).log2_allele_ratio.median().rename("med_robust").reset_index()
g = g.merge(rb, how="left")
g.to_csv(OUT/f"gene_window_ASE_{SRC}.tsv.gz", sep="\t", index=False, compression="gzip")

def summ(x):
    r = x[x.robust_any]
    return dict(n_eligible=len(x), n_robust=len(r), pct_robust=100*len(r)/max(len(x),1),
                pct_B=100*(r.med_robust < 0).mean() if len(r) else np.nan,
                pct_A=100*(r.med_robust > 0).mean() if len(r) else np.nan,
                median_log2AB_robust=r.med_robust.median() if len(r) else np.nan,
                median_log2AB_alleligible=x.med_all.median())

rows = []; nulls = {}
for (an, w), x in g.groupby(["analysis","stage_group"]):
    base = summ(x); rows.append(dict(analysis=an, window=w, set="Genome-wide", **base,
                                     p_perm_lower=np.nan, p_perm_two=np.nan, p_perm_matched_lower=np.nan, p_wilcox=np.nan, p_fisher_B=np.nan))
    r = x[x.robust_any].reset_index(drop=True)
    vals = r.med_robust.values
    dec = pd.qcut(np.log10(r.mean_depth), 10, labels=False, duplicates="drop").values
    mods = mem[mem.analysis == an].groupby("trait_module").gene_id.apply(set)
    for m, gs in mods.items():
        xm = x[x.gene_id.isin(gs)]; s = summ(xm)
        inm = r.gene_id.isin(gs).values; k = inm.sum()
        if k < 3:
            rows.append(dict(analysis=an, window=w, set=m, **s)); continue
        obs = np.median(vals[inm])
        null = np.array([np.median(vals[rng.choice(len(vals), k, replace=False)]) for _ in range(NPERM)])
        # depth-decile matched null
        idx_by_dec = {d: np.where(dec == d)[0] for d in np.unique(dec)}
        need = pd.Series(dec[inm]).value_counts()
        nullm = np.empty(NPERM)
        for i in range(NPERM):
            pick = np.concatenate([rng.choice(idx_by_dec[d], c, replace=False) for d, c in need.items()])
            nullm[i] = np.median(vals[pick])
        c = np.median(null)
        rest = vals[~inm]
        fb = stats.fisher_exact([[int((vals[inm] < 0).sum()), int((vals[inm] > 0).sum())],
                                 [int((rest < 0).sum()), int((rest > 0).sum())]])
        rows.append(dict(analysis=an, window=w, set=m, **s,
                         p_perm_lower=(1 + (null <= obs).sum())/(NPERM + 1),
                         p_perm_two=(1 + (np.abs(null - c) >= abs(obs - c)).sum())/(NPERM + 1),
                         p_perm_matched_lower=(1 + (nullm <= obs).sum())/(NPERM + 1),
                         p_wilcox=stats.mannwhitneyu(vals[inm], rest).pvalue,
                         odds_B=fb[0], p_fisher_B=fb[1]))
res = pd.DataFrame(rows)
res.to_csv(OUT/f"part1_genome_vs_module_{SRC}.tsv", sep="\t", index=False, float_format="%.5g")
pd.set_option("display.width", 250); pd.set_option("display.max_columns", 30)
print(res[["analysis","window","set","n_eligible","n_robust","pct_B","pct_A","median_log2AB_robust","median_log2AB_alleligible","p_perm_lower","p_perm_matched_lower","p_wilcox","p_fisher_B"]].to_string())
