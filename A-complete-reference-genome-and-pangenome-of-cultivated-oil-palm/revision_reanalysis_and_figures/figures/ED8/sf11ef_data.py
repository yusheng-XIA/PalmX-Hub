#!/usr/bin/env python3
"""Source Data for the new Supplementary Fig. 11 panels e (confounding controls for Fig. 5e) and f (Fig. 5f).

Inputs : res/ (cluster outputs of scripts/ec_01_snv_classes.py and scripts/ec_02_confound.py, copied from
         ${CLUSTER_WORK}/enh_dsv_confound/res).
Outputs: sd/SF11e_counts.tsv, sd/SF11e_correlations.tsv, sd/SF11f_dSV_pairs.tsv, sd/SF11f_ratios.tsv
Checks : B0 background (uniform on the chromosome) reproduces Fig. 5f (All38) and Supplementary Fig. 11c
         (definition (c), 1/35) pair by pair; dSNP counts reproduce Fig. 5e / Supplementary Fig. 11b.
"""
import re
from pathlib import Path
import numpy as np, pandas as pd

H = Path(__file__).resolve().parent; R = H / "res"; SD = H / "sd"; SD.mkdir(exist_ok=True)
EB = H.parent / "enh_B/res/coloc_syri"
SETN = {"A608": "All38 (1,480 dSVs; 38 haplotypes)", "S560": "Definition (c), 1/35 (550 dSVs; 35 haplotypes)"}
CLS = {"dSNP": "dSNP (missense, stop, start, splice or conserved-proxy)",
       "NCNC": "Derived non-coding SNV (intron, UTR, intergenic or +-2 kb; outside conserved proxy)",
       "SYNall": "Derived synonymous SNV (all)",
       "SYN": "Derived synonymous SNV outside conserved proxy",
       "NCunpol": "Rare non-coding SNV, polarity not required"}
MET = {"1_raw_count": "Counts", "2_haplotype_demeaned": "Counts, haplotype mean removed",
       "3_two_way_FE": "Counts, haplotype and chromosome effects removed (two-way fixed effects)",
       "4_density_per_Mb": "Per Mb of chromosome", "5_density_per_functional_Mb": "Per Mb of CDS or conserved sequence",
       "5_density_per_functional_Mb_haplotype_demeaned": "Per Mb of CDS or conserved sequence, haplotype mean removed",
       "6_two_way_FE_plus_neutral_counts": "Two-way fixed effects plus derived non-coding and synonymous counts"}
BGN = {"B0": "Random position on the same chromosome (as Fig. 5f)",
       "B1": "Same chromosome and same decile of CDS + conserved bp within +-1 Mb",
       "B2": "As B1, with the background centre inside CDS or conserved sequence"}


def disp(h):
    m = re.fullmatch(r"(dura|pisifera|nrly|bk)_hap([12])", h)
    return f"{dict(dura='TK', pisifera='NS', nrly='Nigerian', bk='TN')[m.group(1)]}-Hap{m.group(2)}" if m else h


# ---------------- e: counts and correlations
cnt = []
for s in ("A608", "S560"):
    c = pd.read_csv(R / f"e_counts_{s}.tsv", sep="\t")
    ref = pd.read_csv(EB / (("V0_All38_formal_1480" if s == "A608" else "EG35_c_Phoenix_and_Oleifera_k1") + ".fig_e_counts.tsv"), sep="\t")
    assert (c.dSV.values == ref.dSV_Count.values).all() and (c.dSNP.values == ref.dSNP_Count.values).all(), s
    cnt.append(pd.DataFrame({"Set": SETN[s], "Haplotype": c.Sample.map(disp), "Chromosome": c.Chrom, "dSV": c.dSV,
                             "dSNP": c.dSNP, "Derived non-coding SNV": c.NCNC, "Derived synonymous SNV": c.SYNall,
                             "Derived synonymous SNV outside conserved proxy": c.SYN,
                             "Rare non-coding SNV, polarity not required": c.NCunpol,
                             "Chromosome Mb": c.Chrom_Mb.round(4), "CDS + conserved Mb": c.Functional_Mb.round(4),
                             "CDS Mb": c.CDS_Mb.round(4)}))
pd.concat(cnt).to_csv(SD / "SF11e_counts.tsv", sep="\t", index=False)
E = pd.read_csv(R / "e_stats.tsv", sep="\t")
E = E[E.Set.isin(SETN)]
rows = []
for q in E.itertuples():
    if q.Metric.startswith("7_") or q.Y.startswith("dSNP_thinned"):
        continue
    if q.Metric not in MET:
        continue
    rows.append({"Set": SETN[q.Set], "n": q.n, "SNV class": CLS[q.Y], "Adjustment": MET[q.Metric],
                 "Pearson_r": round(q.r, 4), "P": f"{q.P:.3e}", "Residual_df": q.df})
for q in E[E.Y.str.startswith("dSNP_thinned")].itertuples():
    if q.Metric in ("1_raw_count", "3_two_way_FE"):
        y = q.Y.replace("dSNP_thinned_to_", "")
        lo = E[(E.Set == q.Set) & (E.Y == q.Y) & (E.Metric == "3_two_way_FE_q025")].r.item() if q.Metric == "3_two_way_FE" else np.nan
        hi = E[(E.Set == q.Set) & (E.Y == q.Y) & (E.Metric == "3_two_way_FE_q975")].r.item() if q.Metric == "3_two_way_FE" else np.nan
        rows.append({"Set": SETN[q.Set], "n": q.n, "SNV class": f"dSNP randomly thinned to the number of '{CLS[y]}' sites ({int(q.df):,}); mean of 200 draws"
                     + (f" (2.5-97.5%: {lo:.3f}-{hi:.3f})" if q.Metric == "3_two_way_FE" else ""),
                     "Adjustment": MET[q.Metric], "Pearson_r": round(q.r, 4), "P": "", "Residual_df": ""})
for q in E[E.Metric.str.startswith("7_")].itertuples():
    y = q.Y.replace("dSNP_vs_", "")
    rows.append({"Set": SETN[q.Set], "n": q.n, "SNV class": f"r(dSV, dSNP) minus r(dSV, {CLS[y]})",
                 "Adjustment": ("Two-way fixed effects" if "FE" in q.Metric else "Counts") + "; Meng-Rosenthal-Rubin test",
                 "Pearson_r": round(q.r, 4), "P": f"{q.P:.3e}", "Residual_df": q.df})
pd.DataFrame(rows).to_csv(SD / "SF11e_correlations.tsv", sep="\t", index=False)

# ---------------- f: per-focal pairs and ratios
pp = []
for s in ("A608", "S560"):
    p = pd.read_csv(R / f"f_pairs_{s}.tsv", sep="\t")
    ref = pd.read_csv(EB / (("V0_All38_formal_1480" if s == "A608" else "EG35_c_Phoenix_and_Oleifera_k1") + ".fig_f_pairs.tsv"), sep="\t")
    assert (p.Focal_DSV_ID.values == ref.Focal_DSV_ID.values).all()
    assert np.allclose(p.Obs_dSNP, ref.Observed) and np.allclose(p.B0_dSNP, ref.Background), s
    o = pd.DataFrame({"Set": SETN[s], "Focal dSV": p.Focal_DSV_ID, "Chromosome": p.Chrom, "Focal position": p.Focal_Pos,
                      "Carrier haplotype(s)": p.Carriers.map(lambda x: ",".join(disp(h) for h in x.split(","))),
                      "CDS + conserved bp within +-1 Mb": p.Window_functional_bp, "Within-chromosome decile (0-9)": p.Decile})
    for y, lab in (("dSNP", "dSNP"), ("NCNC", "Derived non-coding SNV"), ("SYNall", "Derived synonymous SNV"),
                   ("NCunpol", "Rare non-coding SNV, polarity not required")):
        o[f"{lab}: observed"] = p[f"Obs_{y}"]
        for b in ("B0", "B1", "B2"):
            o[f"{lab}: background {b}"] = p[f"{b}_{y}"].round(4)
    pp.append(o)
pd.concat(pp).to_csv(SD / "SF11f_dSV_pairs.tsv", sep="\t", index=False)
F = pd.read_csv(R / "f_stats.tsv", sep="\t")
F = F[F.Set.isin(SETN)]
rows = []
for q in F.itertuples():
    if q.Y.startswith("dSNP_over_"):
        y = q.Y.replace("dSNP_over_", "")
        rows.append({"Set": SETN[q.Set], "Focal dSVs": q.Focal_dSVs, "SNV class": f"dSNP ratio / {CLS[y]} ratio",
                     "Background": f"{q.Background}: {BGN[q.Background]}", "Observed mean": "", "Background mean": "",
                     "Ratio": round(q.Ratio, 4), "Ratio 95% CI low": round(q.Ratio_lo95, 4), "Ratio 95% CI high": round(q.Ratio_hi95, 4),
                     "Test": "paired bootstrap of focal dSVs (1,000), two-sided", "P": f"{q.P:.3g}"})
    else:
        rows.append({"Set": SETN[q.Set], "Focal dSVs": q.Focal_dSVs, "SNV class": CLS[q.Y],
                     "Background": f"{q.Background}: {BGN[q.Background]}", "Observed mean": round(q.Obs_mean, 4),
                     "Background mean": round(q.Bg_mean, 4), "Ratio": round(q.Ratio, 4), "Ratio 95% CI low": round(q.Ratio_lo95, 4),
                     "Ratio 95% CI high": round(q.Ratio_hi95, 4), "Test": "two-sided paired Wilcoxon signed-rank",
                     "P": f"{q.P:.3e}"})
pd.DataFrame(rows).to_csv(SD / "SF11f_ratios.tsv", sep="\t", index=False)
print("written", sorted(x.name for x in SD.glob("SF11[ef]*")))
