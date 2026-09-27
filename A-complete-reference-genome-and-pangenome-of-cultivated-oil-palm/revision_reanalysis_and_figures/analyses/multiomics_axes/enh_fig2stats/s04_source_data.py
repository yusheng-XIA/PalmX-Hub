#!/usr/bin/env python3
"""Source Data TSVs for Supplementary Fig. 2g and the Fig. 2d statistics (from out/)."""
from pathlib import Path
import pandas as pd
H = Path(__file__).resolve().parent; O = H / "out"; SD = H / "source_data"; SD.mkdir(exist_ok=True)
summ = dict(l.split("\t", 1) for l in (O / "dosage_null_summary.tsv").read_text().splitlines()[1:])
enz = pd.read_csv(O / "dosage_enzyme_matching.tsv", sep="\t")
bg = pd.read_csv(O / "dosage_background_by_copy_number.tsv", sep="\t")
eb = enz.groupby("cn_bin").identical_Cocos_and_3_Elaeis.agg(["sum", "size"])
bg["lipid_classes_n"] = bg.cn_bin.map(eb["size"]).fillna(0).astype(int)
bg["lipid_classes_identical"] = bg.cn_bin.map(eb["sum"]).fillna(0).astype(int)
bg = bg.rename(columns={"cn_bin": "Cocos_copy_number", "n_OG": "orthogroups_n",
                        "frac_identical_Cocos_3Elaeis": "orthogroups_identical_Cocos_and_3_Elaeis_fraction",
                        "frac_identical_3Elaeis": "orthogroups_identical_3_Elaeis_fraction"})
tot = dict(Cocos_copy_number="all", orthogroups_n=int(summ["locus.U_all4.n_OG"]),
           orthogroups_identical_Cocos_and_3_Elaeis_fraction=float(summ["locus.U_all4.frac_identical_Cocos_3Elaeis"]),
           orthogroups_identical_3_Elaeis_fraction=float(summ["locus.U_all4.frac_identical_3Elaeis"]),
           lipid_classes_n=24, lipid_classes_identical=17)
bg = pd.concat([bg, pd.DataFrame([tot])])
bg.to_csv(SD / "SF2g_background_by_copy_number.tsv", sep="\t", index=False, float_format="%.4f")
dist = pd.read_csv(O / "dosage_null_distribution.tsv", sep="\t")
dist.to_csv(SD / "SF2g_null_distribution.tsv", sep="\t", index=False)
keep = [k for k in summ if not k.startswith("protein_model")]
pd.DataFrame({"statistic": keep, "value": [summ[k] for k in keep]}).to_csv(SD / "SF2g_test_summary.tsv", sep="\t", index=False)
enz.drop(columns=["pool_B"]).to_csv(SD / "SF2g_enzyme_matching.tsv", sep="\t", index=False, float_format="%.4f")
st = pd.read_csv(O / "fig2d_stage_tests.tsv", sep="\t")
st = st[["axis_id", "axis", "stage", "phase", "FL_mean", "TN_mean", "FL_sd", "TN_sd", "diff_FL_minus_TN", "higher",
         "welch_t", "welch_df", "welch_P", "welch_q_axis"]]
an = pd.read_csv(O / "fig2d_twoway_anova.tsv", sep="\t")
an = an[~an.term.str.contains("HC3")]
st.to_csv(SD / "Fig2d_stage_tests.tsv", sep="\t", index=False, float_format="%.4g")
an.to_csv(SD / "Fig2d_twoway_anova.tsv", sep="\t", index=False, float_format="%.4g")
print(bg.to_string()); print(an.to_string())
