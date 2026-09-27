#!/usr/bin/env python3
"""Source Data tables for the new Supplementary Fig. 5d-i (OLE16a/OLE16b loci), from fix/idea_A/res and
fix/ase_genogroup/in/gene_stage_ASE_unified.tsv.gz.  Writes sd/*.tsv and prints the numbers quoted in the text."""
import gzip
from pathlib import Path

import numpy as np
import pandas as pd
from Bio import Phylo

HERE = Path(__file__).resolve().parent
A = HERE.parent
RES = A / "res"
OUT = HERE / "sd"
OUT.mkdir(exist_ok=True)
ST = ["0d", "15d", "35d", "50d", "65d", "80d", "95d", "110d", "125d", "140d", "155d", "170d", "185d",
      "12h", "24h", "36h", "48h", "60h", "72h"]
PH = ["185d", "12h", "24h", "36h", "48h", "60h", "72h"]
OLE16A = {"evm.TU.chr11B.1497", "evm.TU.chr11A.1059", "evm.TU.bk_hap1_chr9.1211", "evm.TU.bk_hap2_chr9.414"}
OLE16B = {"evm.TU.chr04B.697", "evm.TU.chr04A.610", "evm.TU.bk_hap1_chr3.2082"}

# ---- d: phylogeny (full IQ-TREE consensus tree; pruned tip set drawn) ---------------------------------
KEEP = {  # tip-name substring -> display label
    "chr11B.1497": "OLE16a FL-Hap2 (chr11B)", "chr11A.1059": "OLE16a FL-Hap1 (chr11A)",
    "bk_hap1_chr9.1211": "OLE16a TN-Hap1", "bk_hap2_chr9.414": "OLE16a TN-Hap2",
    "chr04B.697": "OLE16b FL-Hap2 (chr04B)", "chr04A.610": "OLE16b FL-Hap1 (chr04A)",
    "bk_hap1_chr3.2082": "OLE16b TN-Hap1",
    "Pda_LOC103704664": "Date palm LOC103704664", "Pda_LOC103700664": "Date palm LOC103700664",
    "Osa_OLE16": "Rice OLE16", "Zma_OLE16": "Maize OLE16", "Bse_OLE16": "Brome OLE16",
    "Ath_OLE1_": "Arabidopsis OLE1", "Ath_OLE3_": "Arabidopsis OLE3", "Pdu_OLE1": "Almond OLE1",
    "Sin_OLEL": "Sesame OLE-L", "Osa_OLE18": "Rice OLE18", "Zma_OLE18": "Maize OLE18",
    "Sin_OLEH1": "Sesame OLE-H1", "Ath_GRP17": "Arabidopsis GRP17", "Ath_GRP19": "Arabidopsis GRP19",
    "Pam_MRB53_021998": "Avocado MRB53_021998", "Pam_MRB53_006675": "Avocado MRB53_006675",
    "Pam_MRB53_024459": "Avocado MRB53_024459", "Llo_OLE": "Lily pollen oleosin",
    "Pel_OLE": "Pine oleosin", "Ath_At3g18570": "Arabidopsis At3g18570",
    "chr15B.53": "Oil palm chr15B.53", "chr06B.757": "Oil palm chr06B.757", "chr13B.967": "Oil palm chr13B.967",
}
t = Phylo.read(str(RES / "dom/oleosin_iq.contree"), "newick")
t.root_at_midpoint()
full_nwk = t.format("newick").strip()
rows = []
for c in t.get_terminals():
    k = next((k for k in KEEP if k in c.name), None)
    p = c.name.split("|")
    grp = p[1]
    if "OLE16a" in (KEEP.get(k) or ""):
        grp = "OLE16a"
    elif "OLE16b" in (KEEP.get(k) or ""):
        grp = "OLE16b"
    rows.append(dict(Tip=c.name, Category=grp, Shown_in_panel=("yes" if k else "no"), Display_label=KEEP.get(k, "")))
for c in t.get_terminals():
    if not any(k in c.name for k in KEEP):
        t.prune(c)
assert len(t.get_terminals()) == len(KEEP), (len(t.get_terminals()), len(KEEP))
for c in t.get_terminals():
    c.name = KEEP[next(k for k in KEEP if k in c.name)]
t.ladderize()
Phylo.write(t, str(HERE / "tree_pruned.nwk"), "newick")
d = pd.DataFrame(rows)
d = pd.concat([d, pd.DataFrame([dict(Tip="Full consensus tree (Newick, midpoint-rooted; UFBoot support)",
                                     Category="", Shown_in_panel="", Display_label=full_nwk)])])
d.to_csv(OUT / "SF5d_oleosin_phylogeny.tsv", sep="\t", index=False)

# ---- e: bulk RNA, OLE16a and OLE16b, per sample ----------------------------------------------------
L = pd.read_csv(RES / "rna/ld_gene_sample_long.tsv.gz", sep="\t")
e = L[L.gene.isin(["evm.TU.chr11B.1497", "evm.TU.chr04B.697"])].copy()
e["Locus"] = e.gene.map({"evm.TU.chr11B.1497": "OLE16a", "evm.TU.chr04B.697": "OLE16b"})
e["Stage index"] = e.stage.map({s: i + 1 for i, s in enumerate(ST)})
e = e.sort_values(["Locus", "genotype", "Stage index", "sample"])
e = e[["Locus", "gene", "genotype", "Stage index", "stage", "sample", "norm"]].rename(
    columns={"gene": "Gene (FL-Hap2)", "genotype": "Material", "stage": "Stage", "sample": "Library",
             "norm": "Normalized count (DESeq2 size factors, 114 libraries)"})
e.to_csv(OUT / "SF5e_OLE16_RNA.tsv", sep="\t", index=False)
em = L[L.gene.isin(["evm.TU.chr11B.1497"])].groupby(["genotype", "stage"]).norm.mean()

# ---- f: shared-peptide precursor quantities per sample ------------------------------------------
P = pd.read_csv(RES / "ld_precursor_rows.tsv.gz", sep="\t", low_memory=False)
shared = {"RVPGSEQLEQAR": "OLE16a", "VPGSEQLEQAR": "OLE16a", "RPPGFEQLEQAR": "OLE16b"}
f = P[P["Stripped.Sequence"].isin(shared) & (P.accepted.astype(str) == "True")].copy()
f["Locus"] = f["Stripped.Sequence"].map(shared)
k = f.key.str.split("|", expand=True)
f["Material"], f["Stage"], f["Replicate"] = k[0], k[1], k[2]
f = f[f.Stage.isin(ST[10:])]
f["Stage index"] = f.Stage.map({s: i + 1 for i, s in enumerate(ST)})
f = f.sort_values(["Locus", "Precursor.Id", "Material", "Stage index", "Replicate"])
f = f[["Locus", "Stripped.Sequence", "Precursor.Id", "Material", "Stage index", "Stage", "Replicate",
       "Precursor.Quantity", "Q.Value", "prots"]].rename(columns={"prots": "Mapped protein entries"})
f.to_csv(OUT / "SF5f_OLE16_shared_peptides.tsv", sep="\t", index=False)
Q = pd.read_csv(RES / "run_precursor_quantiles.tsv", sep="\t")
Q.to_csv(OUT / "SF5f_run_precursor_quantiles.tsv", sep="\t", index=False)
fm = f.groupby(["Precursor.Id", "Material", "Stage"])["Precursor.Quantity"].agg(["mean", "count"])

# ---- g: LD-coat protein families, 185 d-72 h -----------------------------------------------------
D = pd.read_csv(RES / "prot/ld_family_stage_means_1e6.tsv", sep="\t")
grows = []
for fam, lab in [("oleosin", "Oleosin"), ("REF", "LDAP (REF/SRPP)"), ("caleosin", "Caleosin")]:
    for m in ("FL", "TN"):
        for s in PH:
            grows.append(dict(Family=lab, Material=m, Stage=s,
                              **{"Summed directLFQ abundance (mean of 3 replicates)":
                                 D[D.family == fam][f"{m}{s}"].sum() * 1e6}))
g = pd.DataFrame(grows)
g.to_csv(OUT / "SF5g_LD_coat_families.tsv", sep="\t", index=False)
D.to_csv(OUT / "SF5g_LD_coat_protein_groups.tsv", sep="\t", index=False)

# ---- h: FL allele origin of OLE16a transcripts --------------------------------------------------
with gzip.open(A.parent / "ase_genogroup/in/gene_stage_ASE_unified.tsv.gz", "rt") as fh:
    U = pd.read_csv(fh, sep="\t")
h = U[(U.gene_id == "evm.TU.chr11B.1497") & (U.analysis == "FL")].copy()
h = h[h.informative_fragments.fillna(0) > 0]
h["FL-Hap1 fraction"] = h.alt_fragments / h.informative_fragments
h = h[["gene_id", "stage_index", "stage", "qualifying_replicates", "ref_fragments", "alt_fragments",
       "informative_fragments", "FL-Hap1 fraction", "robust_ase"]].rename(columns={
    "gene_id": "Gene", "stage_index": "Stage index", "stage": "Stage", "ref_fragments": "FL-Hap2 fragments",
    "alt_fragments": "FL-Hap1 fragments", "informative_fragments": "Allele-informative fragments"})
h.to_csv(OUT / "SF5h_OLE16a_FL_allele_origin.tsv", sep="\t", index=False)

# ---- i: snRNA, OLE16a-positive nuclei per cluster at 185 d ----------------------------------------
R = pd.read_csv(RES / "sn/ld_genes_by_cluster.tsv", sep="\t")
R = R[(R.gene == "evm.TU.chr11B.1497") & R.lib.isin(["FL_185", "TN_185"])].copy()
R["Cluster"] = "C" + R.cl.astype(str)
R = R[["lib", "Cluster", "n", "pct_pos"]].rename(columns={"lib": "Library", "n": "Nuclei",
                                                          "pct_pos": "% nuclei with >=1 OLE16a UMI"})
R.to_csv(OUT / "SF5i_OLE16a_snRNA_clusters.tsv", sep="\t", index=False)
S = pd.read_csv(RES / "sn/ole16a_positive_nuclei_summary.tsv", sep="\t")
S.to_csv(OUT / "SF5i_OLE16a_snRNA_libraries.tsv", sep="\t", index=False)

# ---- numbers for the text ------------------------------------------------------------------------
print("RNA OLE16a FL 155/170/185:", [round(em["FL", s]) for s in ("155d", "170d", "185d")],
      "FL <=140d max", round(max(em["FL", s] for s in ST[:10]), 1), "TN max", round(em["TN"].max(), 1))
for p in ("RVPGSEQLEQAR3", "VPGSEQLEQAR2", "RPPGFEQLEQAR3"):
    a, b = fm.loc[(p, "FL", "185d")], fm.loc[(p, "TN", "185d")]
    print(p, "185d FL", a["mean"], int(a["count"]), "TN", b["mean"], int(b["count"]), "FL/TN", a["mean"] / b["mean"])
gg = g.set_index(["Family", "Material", "Stage"]).iloc[:, 0]
print("oleosin:LDAP 185d FL", gg["Oleosin", "FL", "185d"] / gg["LDAP (REF/SRPP)", "FL", "185d"],
      "TN", gg["Oleosin", "TN", "185d"] / gg["LDAP (REF/SRPP)", "TN", "185d"])
print("OLE16a directLFQ 185d FL/TN", D.set_index("genes").loc["evm.TU.chr11B.1497", "FL185d"] /
      D.set_index("genes").loc["evm.TU.chr11B.1497", "TN185d"])
hh = h.set_index("Stage")
print("FL-Hap1 fraction:", hh["FL-Hap1 fraction"].round(4).to_dict())
print("min 170d-72h", hh.loc[ST[11:], "FL-Hap1 fraction"].min())
print(S[["lib", "n", "OLEa_pos"]].to_string())
