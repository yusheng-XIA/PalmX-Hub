"""Build Source-Data-style tables, the SD23 add-on columns and the summary counts from s10/s11 outputs (local, small files)."""
import numpy as np, pandas as pd
from pathlib import Path
B = Path(__file__).resolve().parents[1]; O = B / "out"; R = B / "results"; R.mkdir(exist_ok=True)
D = pd.read_csv(O / "finemap_all.tsv", sep="\t")
A = pd.read_csv(O / "sv_annotation_FL_Hap2.tsv", sep="\t").set_index("sv")
F = "SNPmodel_f0.2"
def fmt(v): return np.nan if pd.isna(v) else float(f"{v:.4g}")
def build(d):
    a = A.reindex(d.sv_lead)
    t = pd.DataFrame({
        "trait": d.trait.values, "locus_id": d.locus_id.values, "chrom": d.chrom.values,
        "SV_interval_far_from_SNP_interval": d.get("far_from_SNP_interval", pd.Series([np.nan] * len(d))).values,
        "n": d.n.values, "SNPs_in_window": d.n_snp_window.values, "SVs_in_window": d.n_sv_window.values,
        "lead_SV": d.sv_lead.values, "lead_SV_pos": d.sv_lead_pos.values, "lead_SV_type": a.sv_type.values,
        "lead_SV_net_length_bp": a.sv_len.values, "lead_SV_ref_allele_bp": a.ref_len.values, "lead_SV_alt_allele_bp": a.alt_len.values,
        "lead_SV_MAC": d.sv_lead_MAC.values, "lead_SV_P_SV_model": d.sv_lead_P_SVmodel.map(fmt).values,
        "lead_SV_P_SNP_model": d.sv_lead_P_SNPmodel.map(fmt).values, "lead_SV_beta_SD": d.sv_lead_beta_SD_SNPmodel.map(fmt).values,
        "lead_SV_se_SD": d.sv_lead_se_SD_SNPmodel.map(fmt).values,
        "lead_SNP": d.snp_lead.values, "lead_SNP_MAC": d.snp_lead_MAC.values, "lead_SNP_P_SNP_model": d.snp_lead_P_SNPmodel.map(fmt).values,
        "lead_SNP_beta_SD": d.snp_lead_beta_SD_SNPmodel.map(fmt).values, "lead_SNP_se_SD": d.snp_lead_se_SD_SNPmodel.map(fmt).values,
        "r2_lead_SNP_lead_SV": d.r2_snp_sv.map(fmt).values, "distance_lead_SNP_lead_SV_bp": d.dist_snp_sv.values,
        "P_lead_SNP_given_lead_SV": d.condA_P_snp_given_sv.map(fmt).values, "P_lead_SV_given_lead_SNP": d.condB_P_sv_given_snp.map(fmt).values,
        "P_lead_SNP_given_lead_SV_SV_model": d.condA2_P_snp_given_sv_SVmodel.map(fmt).values,
        "P_lead_SV_given_lead_SNP_SNP_model": d.condB2_P_sv_given_snp_SNPmodel.map(fmt).values,
        "lead_SNP_retained_fraction_log10P": d.snp_retained_frac.map(fmt).values, "lead_SV_retained_fraction_log10P": d.sv_retained_frac.map(fmt).values,
        "reciprocal_conditioning_class": d.class_primary.values, "class_common_SNP_model": d.class_commonmodel_SNP.values,
        "class_common_SV_model": d.class_commonmodel_SV.values,
        "ABF_PP_lead_SV": d[f"PP_svlead_{F}"].map(fmt).values, "ABF_PP_lead_SNP": d[f"PP_snplead_{F}"].map(fmt).values,
        "ABF_PP_all_SVs": d[f"PP_allSV_{F}"].map(fmt).values, "CS95_size": d[f"CS95_size_{F}"].values, "CS95_n_SV": d[f"CS95_nSV_{F}"].values,
        "lead_SV_in_CS95": d[f"svlead_in_CS95_{F}"].values, "top_PP_variant": d[f"top_variant_{F}"].values, "top_PP": d[f"top_PP_{F}"].map(fmt).values,
        "ABF_PP_lead_SV_prior0.1SD": d["PP_svlead_SNPmodel_f0.1"].map(fmt).values, "ABF_PP_lead_SV_prior0.4SD": d["PP_svlead_SNPmodel_f0.4"].map(fmt).values,
        "ABF_PP_lead_SV_SV_model": d["PP_svlead_SVmodel_f0.2"].map(fmt).values, "lead_SV_in_CS95_SV_model": d["svlead_in_CS95_SVmodel_f0.2"].values,
        "CS95_SVs_with_PP": d["CS95_SVs_SNPmodel_f0.2"].values,
        "lead_SV_FL_Hap2_feature": a.feature_class.values, "lead_SV_CDS_genes": a.CDS_genes.values, "lead_SV_intron_genes": a.intron_genes.values,
        "lead_SV_promoter_2kb_genes": a.promoter_2kb_genes.values, "lead_SV_nearest_gene": a.nearest_gene.values,
        "lead_SV_nearest_gene_distance_bp": a.nearest_gene_distance_bp.values, "lead_SV_gene_function_eggNOG": a.nearest_gene_function.values})
    return t
sv = D[D.anchor == "SV"]; sn = D[D.anchor == "SNP"]
Tsv = build(sv); Tsn = build(sn).drop(columns=["SV_interval_far_from_SNP_interval"])
Tsv.to_csv(R / "SourceData_SV_intervals_reciprocal_conditioning_finemap.tsv", sep="\t", index=False)
Tsn.to_csv(R / "SourceData_SNP_peaks_reciprocal_conditioning_finemap.tsv", sep="\t", index=False)
S = pd.read_csv(O / "shell_snp_lead_conditioned_on_each_SV.tsv", sep="\t")
S = S.join(A[["feature_class", "nearest_gene", "nearest_gene_distance_bp"]], on="sv")
for c in ["snp_lead_P", "sv_P_SNPmodel", "sv_P_SVmodel", "r2_with_snp_lead", "P_snp_lead_given_sv", "snp_retained_frac"]: S[c] = S[c].map(fmt)
S.to_csv(R / "SourceData_SHELL_region_SNP_lead_conditioned_on_each_SV.tsv", sep="\t", index=False)
# SD23 add-on columns (keyed by SV_interval)
sd = pd.read_csv(B.parent / "enh_gwas/SD23_new_SV_intervals.tsv", sep="\t")
add = Tsv.set_index("locus_id").loc[sd.SV_interval, ["lead_SNP", "lead_SNP_P_SNP_model", "r2_lead_SNP_lead_SV", "P_lead_SNP_given_lead_SV",
       "P_lead_SV_given_lead_SNP", "reciprocal_conditioning_class", "ABF_PP_lead_SV", "lead_SV_in_CS95", "CS95_size",
       "lead_SV_type", "lead_SV_FL_Hap2_feature", "lead_SV_nearest_gene", "lead_SV_nearest_gene_distance_bp", "lead_SV_gene_function_eggNOG"]]
add = add.rename(columns={"lead_SNP": "top_SNP_within_500kb", "lead_SNP_P_SNP_model": "top_SNP_P", "r2_lead_SNP_lead_SV": "r2_top_SNP_lead_SV",
                          "P_lead_SNP_given_lead_SV": "P_top_SNP_given_lead_SV", "P_lead_SV_given_lead_SNP": "P_lead_SV_given_top_SNP",
                          "ABF_PP_lead_SV": "posterior_probability_lead_SV", "lead_SV_in_CS95": "lead_SV_in_95pct_credible_set",
                          "CS95_size": "credible_set_size", "lead_SV_FL_Hap2_feature": "FL_Hap2_feature",
                          "lead_SV_nearest_gene": "FL_Hap2_gene", "lead_SV_nearest_gene_distance_bp": "distance_to_gene_bp",
                          "lead_SV_gene_function_eggNOG": "gene_function_eggNOG"})
add.insert(0, "SV_interval", sd.SV_interval.values)
add.to_csv(R / "SD23_addon_columns.tsv", sep="\t", index=False)
print("SD23 rows", len(add))
# summary counts
rows = []
def add_set(name, d):
    d2 = d[~d.class_primary.isin(["no_SNP_in_window", "no_SV_in_window"])]
    r = dict(set=name, n_loci=len(d), n_with_both_classes=len(d2))
    for c in ["SV_explains_SNP", "mutual", "SNP_explains_SV", "neither"]: r[c] = int((d2.class_primary == c).sum())
    r["SV_explains_SNP_common_SNP_model"] = int((d2.class_commonmodel_SNP == "SV_explains_SNP").sum())
    r["SV_explains_SNP_common_SV_model"] = int((d2.class_commonmodel_SV == "SV_explains_SNP").sum())
    r["SV_conditioning_removes_larger_fraction"] = int((d2.snp_retained_frac < d2.sv_retained_frac).sum())
    r["median_r2_lead_SNP_SV"] = round(d2.r2_snp_sv.median(), 2)
    r["lead_SV_more_significant_same_model"] = int((d2.sv_lead_P_SNPmodel < d2.snp_lead_P_SNPmodel).sum())
    r["lead_SV_larger_abs_effect"] = int((d2.sv_lead_beta_SD_SNPmodel.abs() > d2.snp_lead_beta_SD_SNPmodel.abs()).sum())
    r["median_abs_beta_SD_SV"] = round(d2.sv_lead_beta_SD_SNPmodel.abs().median(), 2); r["median_abs_beta_SD_SNP"] = round(d2.snp_lead_beta_SD_SNPmodel.abs().median(), 2)
    for ff in ["SNPmodel_f0.2", "SNPmodel_f0.1", "SNPmodel_f0.4", "SVmodel_f0.2"]:
        r[f"leadSV_in_CS95_{ff}"] = int(d2[f"svlead_in_CS95_{ff}"].sum()); r[f"top_PP_is_SV_{ff}"] = int(d2[f"top_is_SV_{ff}"].sum())
        r[f"leadSV_PP_ge_0.5_{ff}"] = int((d2[f"PP_svlead_{ff}"] >= 0.5).sum()); r[f"median_CS95_size_{ff}"] = d2[f"CS95_size_{ff}"].median()
        r[f"median_PP_all_SVs_{ff}"] = round(d2[f"PP_allSV_{ff}"].median(), 4)
    r["median_SV_share_of_variants"] = round((d2.n_sv_window / (d2.n_sv_window + d2.n_snp_window)).median(), 4)
    rows.append(r)
add_set("SV intervals (all 115)", sv); add_set("SV intervals >250 kb from SNP interval (47)", sv[sv.far_from_SNP_interval == True])
add_set("SV intervals <=250 kb (68)", sv[sv.far_from_SNP_interval != True])
add_set("SNP peaks, 8 SV traits + nut weight (785)", sn); add_set("SNP peaks with co-located SV P<1e-5", sn[sn.sv_lead_P_SVmodel < 1e-5])
add_set("SNP peaks, excluding lauric acid", sn[sn.trait != "C12_0_Lauric_acid"])
pd.DataFrame(rows).T.to_csv(R / "summary_counts.tsv", sep="\t", header=False)
print(pd.DataFrame(rows).T.to_string())
