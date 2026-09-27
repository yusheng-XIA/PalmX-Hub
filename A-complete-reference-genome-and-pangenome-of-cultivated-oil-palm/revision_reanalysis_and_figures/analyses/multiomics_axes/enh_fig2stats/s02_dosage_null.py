#!/usr/bin/env python3
"""Genome-wide null for 'identical copy number in Cocos nucifera and the three Elaeis genomes' (SF2a, 17/24).

Inputs
  data/og_locus_counts.tsv      per-orthogroup gene-locus counts for the 30 species of the Fig. 2 OrthoFinder run
                                (Results_Feb05, v3.1.1); Cocos protein models collapsed to loci (suffix .N stripped;
                                the other palm proteomes have one model per locus). Built on ${COMPUTE_HOST} by s01_og_counts.sh.
  data/Orthogroups.GeneCount.tsv raw OrthoFinder protein-model counts (reported for comparison only)
  ../fig2_time/optionB/SF2_palA/SF2a_copy_number.tsv   24 curated enzyme classes (SF2a)
  ../../trace/SF/sf2/Fig2c_core_pathway_gene_assignments.tsv   703 curated loci -> orthogroup

Focal genomes: Cocos_nucifera, American_hap1 (= FL-Hap1), Dura (= TK), Pisifera (= NS).

Background universes
  U_any     orthogroups with >= 1 locus in at least one focal genome
  U_all4    orthogroups with >= 1 locus in each of the four focal genomes (all 24 SF2a classes satisfy this)
Matching (U_all4)
  copy number : Cocos locus count, bins 1,2,3,4,5-6,7-10,>10 (the SF2a colour bins)
  family size : total loci across the 30 species (for an enzyme class: summed over its linked orthogroups),
                bins by quintiles of U_all4
Permutation: 10,000 draws; for each of the 24 classes one orthogroup is drawn at random from its stratum;
statistic = number of the 24 draws with identical counts in the four focal genomes. P_upper = P(X >= 17),
P_lower = P(X <= 17) (each with +1 correction); two-sided = min(1, 2 * min).
"""
from pathlib import Path
import json
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
SP = HERE.parents[1]
OUT = HERE / "out"
OUT.mkdir(exist_ok=True)
FOC = ["Cocos_nucifera", "American_hap1", "Dura", "Pisifera"]
EL = ["American_hap1", "Dura", "Pisifera"]
NPERM = 10_000
rng = np.random.default_rng(20260924)

og = pd.read_csv(HERE / "data/og_locus_counts.tsv", sep="\t", index_col=0)
raw = pd.read_csv(HERE / "data/Orthogroups.GeneCount.tsv", sep="\t", index_col=0)
species = [c for c in og.columns]
assert len(species) == 30, len(species)
og["Total"] = og[species].sum(axis=1)
sf = pd.read_csv(SP / "fix/fig2_time/optionB/SF2_palA/SF2a_copy_number.tsv", sep="\t")
sf = sf.rename(columns={})
ga = pd.read_csv(SP / "trace/SF/sf2/Fig2c_core_pathway_gene_assignments.tsv", sep="\t")
ga = ga[ga.Accepted == True]
assert len(sf) == 24 and len(ga) == 703

def ident(df, cols):
    return df[cols].nunique(axis=1).eq(1)

CN_BINS = [0.5, 1.5, 2.5, 3.5, 4.5, 6.5, 10.5, np.inf]
CN_LAB = ["1", "2", "3", "4", "5-6", "7-10", ">10"]

res = {}
lines = []
def rep(k, v):
    res[k] = v
    lines.append(f"{k}\t{v}")

# ---- genome-wide rates ---------------------------------------------------------------------------
for name, tab in (("locus", og), ("protein_model", raw)):
    anyf = tab[FOC].gt(0).any(axis=1)
    all4 = tab[FOC].gt(0).all(axis=1)
    for uname, u in (("U_any", anyf), ("U_all4", all4)):
        sub = tab[u]
        rep(f"{name}.{uname}.n_OG", int(len(sub)))
        rep(f"{name}.{uname}.frac_identical_Cocos_3Elaeis", float(ident(sub, FOC).mean()))
        rep(f"{name}.{uname}.frac_identical_3Elaeis", float(ident(sub, EL).mean()))

U = og[og[FOC].gt(0).all(axis=1)].copy()
U["cn_bin"] = pd.cut(U["Cocos_nucifera"], CN_BINS, labels=CN_LAB).astype(str)
U["idC"] = ident(U, FOC)
U["idE"] = ident(U, EL)
q = np.quantile(U["Total"], [0.2, 0.4, 0.6, 0.8])
FS_BINS = [0] + list(q) + [np.inf]
U["fs_bin"] = pd.cut(U["Total"], FS_BINS, labels=[f"Q{i}" for i in range(1, 6)], include_lowest=True).astype(str)
rep("family_size_quintile_edges", [float(x) for x in q])

strat = (U.groupby("cn_bin").agg(n_OG=("idC", "size"), frac_identical_Cocos_3Elaeis=("idC", "mean"),
                                 frac_identical_3Elaeis=("idE", "mean")).reindex(CN_LAB).reset_index())

# ---- enzyme classes ------------------------------------------------------------------------------
sf["cn_bin"] = pd.cut(sf["Cocos_nucifera"], CN_BINS, labels=CN_LAB).astype(str)
ogs_by_enz = ga.dropna(subset=["Orthogroup"]).groupby("Assigned_enzyme").Orthogroup.apply(lambda s: sorted(set(s)))
sf["linked_OGs"] = sf.Enzyme.map(lambda e: ",".join(ogs_by_enz.get(e, [])))
sf["family_size_30sp"] = sf.Enzyme.map(lambda e: int(og.loc[ogs_by_enz.get(e, []), "Total"].sum()))
sf["fs_bin"] = pd.cut(sf["family_size_30sp"], FS_BINS, labels=[f"Q{i}" for i in range(1, 6)],
                      include_lowest=True).astype(str)
sf["idC"] = ident(sf, FOC)
sf["idE"] = ident(sf, EL)
obs = int(sf.idC.sum()); obsE = int(sf.idE.sum())
rep("observed_identical_Cocos_3Elaeis", obs)
rep("observed_identical_3Elaeis", obsE)

def perm(keys, pools, flag):
    draws = np.zeros((NPERM, len(keys)), dtype=bool)
    for j, k in enumerate(keys):
        pool = pools[k]
        draws[:, j] = pool[rng.integers(0, len(pool), NPERM)]
    return draws.sum(axis=1)

def summarise(tag, x, o):
    up = (1 + (x >= o).sum()) / (NPERM + 1)
    lo = (1 + (x <= o).sum()) / (NPERM + 1)
    rep(f"{tag}.null_mean", float(x.mean()))
    rep(f"{tag}.null_median", float(np.median(x)))
    rep(f"{tag}.null_95pct_interval", [int(np.percentile(x, 2.5)), int(np.percentile(x, 97.5))])
    rep(f"{tag}.null_mean_fraction", float(x.mean() / 24))
    rep(f"{tag}.P_upper(X>=obs)", float(up))
    rep(f"{tag}.P_lower(X<=obs)", float(lo))
    rep(f"{tag}.P_two_sided", float(min(1, 2 * min(up, lo))))

null_tab = {}
for flagcol, o, lab in (("idC", obs, "Cocos_3Elaeis"), ("idE", obsE, "3Elaeis")):
    # A: copy-number matched
    pools = {b: g[flagcol].to_numpy() for b, g in U.groupby("cn_bin")}
    xA = perm(sf.cn_bin.tolist(), pools, flagcol)
    summarise(f"A_copy_matched.{lab}", xA, o)
    # B: copy-number + family-size matched (fall back to copy-number stratum if a cell has < 20 OGs)
    poolsB, keysB = {}, []
    for _, r in sf.iterrows():
        cell = U[(U.cn_bin == r.cn_bin) & (U.fs_bin == r.fs_bin)]
        key = (r.cn_bin, r.fs_bin) if len(cell) >= 20 else (r.cn_bin, "any")
        if key not in poolsB:
            poolsB[key] = (cell if key[1] != "any" else U[U.cn_bin == r.cn_bin])[flagcol].to_numpy()
        keysB.append(key)
    xB = perm(keysB, poolsB, flagcol)
    summarise(f"B_copy_family_matched.{lab}", xB, o)
    null_tab[lab] = (xA, xB)
    if lab == "Cocos_3Elaeis":
        sf["pool_A_n"] = [len(pools[k]) for k in sf.cn_bin]
        sf["pool_A_frac_identical"] = [float(pools[k].mean()) for k in sf.cn_bin]
        sf["pool_B"] = ["/".join(k) for k in keysB]
        sf["pool_B_n"] = [len(poolsB[k]) for k in keysB]
        sf["pool_B_frac_identical"] = [float(poolsB[k].mean()) for k in keysB]

# ---- like-for-like: the linked orthogroups themselves (OG units) ---------------------------------
lip = sorted(set(ga.Orthogroup.dropna()))
L = og.loc[lip].copy()
L_all4 = L[L[FOC].gt(0).all(axis=1)]
rep("lipid_linked_OG.n", len(L)); rep("lipid_linked_OG.n_all4", len(L_all4))
rep("lipid_linked_OG.frac_identical_Cocos_3Elaeis(all4)", float(ident(L_all4, FOC).mean()))
rep("lipid_linked_OG.n_identical(all4)", int(ident(L_all4, FOC).sum()))
rep("lipid_linked_OG.frac_identical_3Elaeis(all4)", float(ident(L_all4, EL).mean()))
# copy-matched permutation for the OGs too
L2 = L_all4.copy(); L2["cn_bin"] = pd.cut(L2["Cocos_nucifera"], CN_BINS, labels=CN_LAB).astype(str)
pools = {b: g["idC"].to_numpy() for b, g in U.groupby("cn_bin")}
xL = np.zeros(NPERM, int)
for b in L2.cn_bin:
    xL += pools[b][rng.integers(0, len(pools[b]), NPERM)]
oL = int(ident(L_all4, FOC).sum())
up = (1 + (xL >= oL).sum()) / (NPERM + 1); lo = (1 + (xL <= oL).sum()) / (NPERM + 1)
rep("lipid_linked_OG.copy_matched_null_mean", float(xL.mean()))
rep("lipid_linked_OG.P_upper", float(up)); rep("lipid_linked_OG.P_lower", float(lo))

# ---- outputs ---------------------------------------------------------------------------------------
strat.to_csv(OUT / "dosage_background_by_copy_number.tsv", sep="\t", index=False)
cols = ["Pathway_section", "Enzyme", "N_linked_Orthogroups", "linked_OGs", "Cocos_nucifera", "American_hap1", "Dura",
        "Pisifera", "idC", "idE", "cn_bin", "family_size_30sp", "fs_bin", "pool_A_n", "pool_A_frac_identical",
        "pool_B", "pool_B_n", "pool_B_frac_identical"]
sf[cols].rename(columns={"idC": "identical_Cocos_and_3_Elaeis", "idE": "identical_3_Elaeis"}).to_csv(
    OUT / "dosage_enzyme_matching.tsv", sep="\t", index=False)
xA, xB = null_tab["Cocos_3Elaeis"]
dist = pd.DataFrame({"n_identical_of_24": np.arange(25)})
dist["A_copy_matched_permutations"] = [int((xA == k).sum()) for k in dist.n_identical_of_24]
dist["B_copy_family_matched_permutations"] = [int((xB == k).sum()) for k in dist.n_identical_of_24]
dist.to_csv(OUT / "dosage_null_distribution.tsv", sep="\t", index=False)
U.reset_index()[["Orthogroup"] + FOC + ["Total", "cn_bin", "fs_bin", "idC", "idE"]].rename(
    columns={"idC": "identical_Cocos_and_3_Elaeis", "idE": "identical_3_Elaeis"}).to_csv(
    OUT / "dosage_background_orthogroups_all4.tsv.gz", sep="\t", index=False)
(OUT / "dosage_null_summary.tsv").write_text("key\tvalue\n" + "\n".join(lines) + "\n")
print("\n".join(lines))
print(strat.to_string())
print(sf[["Enzyme", "Cocos_nucifera", "cn_bin", "family_size_30sp", "fs_bin", "pool_B", "pool_B_n", "pool_B_frac_identical"]].to_string())
