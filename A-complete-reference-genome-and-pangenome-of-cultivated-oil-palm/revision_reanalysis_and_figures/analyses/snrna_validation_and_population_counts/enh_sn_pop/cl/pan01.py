#!/usr/bin/env python3
"""Gene families (39-assembly PAV, Fig. 4e/f) and SVs (178,314 high-confidence catalogue, Fig. 5) absent from the
commercial-hybrid haplotypes TK, NS and TN (6 haplotypes) but present in other assemblies."""
import pandas as pd, numpy as np, collections
T = "${ANALYSIS_DIR}/22_answer_reviews/00_ms/03_V3/04_figure4/Fig4_d_i_pan39_material33_singletons_20260811"
SVF = "${ANALYSIS_DIR}/22_answer_reviews/00_ms/03_V3/05_figure/08_SV_hap39_Figure5abcd_20260808/results/final_catalog/oilpalm_hap38.final_highconfidence_sv.tsv"
O = "${CLUSTER_WORK}/enh_sn_pop/pan/"
pav = pd.read_csv(f"{T}/pan39_genome_redraw_20260812/tables/Orthogroups.PAV.genome39.GO_singletons_64577.tsv", sep="\t", index_col=0)
assert pav.shape == (64577, 39)
X = (pav > 0)
COMM = ["dura_hap1", "dura_hap2", "pisifera_hap1", "pisifera_hap2", "bk_hap1", "bk_hap2"]
OLE = ["meizhou4_hap1", "meizhou4_hap2"]; FL = ["American_hap1", "Africa_hap2"]
# accessions of the commercial parental source pools among the 27 primary assemblies (SA-EG, SEA-A; material class table)
SRC_COMM = ["houke_pa_genome", "8_pa_genome", "65_pa_genome", "176_pa_genome", "183_pa_genome"]
f = X.sum(1)
cls = np.select([f == 39, f == 38, (f >= 2) & (f <= 37)], ["Core", "Soft-core", "Shell"], "Cloud")
print("class counts", pd.Series(cls).value_counts().to_dict())
rows = []
def add(name, commcols):
    others = [c for c in X.columns if c not in commcols]
    eg_others = [c for c in others if c not in OLE + FL]
    absent = ~X[commcols].any(axis=1)
    pres = X[others].sum(1); pres_eg = X[eg_others].sum(1)
    for c in ("Shell", "Cloud"):
        m = absent & (cls == c)
        rows.append(dict(reference=name, n_ref_assemblies=len(commcols), klass=c, families_absent_from_ref=int(m.sum()),
                         of_which_in_EG_noncommercial=int((m & (pres_eg >= 1)).sum()),
                         of_which_in_ge2_EG_noncommercial=int((m & (pres_eg >= 2)).sum()),
                         of_which_only_Eoleifera_or_FL=int((m & (pres_eg == 0)).sum()),
                         total_in_class=int((cls == c).sum())))
add("TK+NS+TN haplotypes", COMM)
add("TK+NS+TN + SA-EG/SEA-A primary assemblies", COMM + SRC_COMM)
r = pd.DataFrame(rows); r.to_csv(O + "families_absent_from_commercial.tsv", sep="\t", index=False); print(r.to_string())
# family lists
absent = ~X[COMM].any(axis=1)
pd.DataFrame({"Orthogroup": pav.index[absent & np.isin(cls, ["Shell", "Cloud"])],
              "class": cls[absent & np.isin(cls, ["Shell", "Cloud"])],
              "n_assemblies": f[absent & np.isin(cls, ["Shell", "Cloud"])].values}).to_csv(O + "families_absent_TK_NS_TN.tsv.gz", sep="\t", index=False, compression="gzip")
# ---------------- SVs
sv = pd.read_csv(SVF, sep="\t", usecols=["Cluster_ID", "Chrom", "Start", "End", "SVTYPE", "SVTYPE_Group", "SVLEN_Median_bp", "Sample_Count", "Samples"])
print("SV records", len(sv))
names = collections.Counter(s for x in sv.Samples for s in x.split(";"))
print("sample names", len(names), sorted(names))
CS = {"dura_hap1", "dura_hap2", "pisifera_hap1", "pisifera_hap2", "bk_hap1", "bk_hap2"}
alias = {"EG_dura": None, "EG_pisifera": None}
S = sv.Samples.str.split(";")
inC = S.apply(lambda s: any(x in CS for x in s))
ole_fl = {"meizhou4_hap1", "meizhou4_hap2", "American_hap1"}
n_eg_other = S.apply(lambda s: sum(1 for x in s if x not in CS and x not in ole_fl))
out = []
for t, d in sv.assign(inC=inC, neg=n_eg_other).groupby("SVTYPE_Group"):
    out.append(dict(svtype=t, total=len(d), absent_from_TK_NS_TN=int((~d.inC).sum()),
                    absent_and_in_EG_noncommercial=int(((~d.inC) & (d.neg >= 1)).sum()),
                    absent_and_in_ge2_EG_noncommercial=int(((~d.inC) & (d.neg >= 2)).sum())))
out = pd.DataFrame(out); out.loc[len(out)] = ["All"] + out.iloc[:, 1:].sum().tolist()
out.to_csv(O + "svs_absent_from_commercial.tsv", sep="\t", index=False); print(out.to_string())
