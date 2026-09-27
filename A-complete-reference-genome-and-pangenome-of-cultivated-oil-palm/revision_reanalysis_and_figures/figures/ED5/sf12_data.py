#!/usr/bin/env python3
"""Source Data TSVs for Supplementary Fig. 12 (genome-wide vs trait-module ASE direction; DNA allele-ratio control),
built from the enh_A results (../enh_A/out, produced on the cluster by s02_part1.py, s05_dna.py, s05b_correct.py and
s06_keygenes.py). Only light tabulation is done here (binomial tests, BH adjustment, reshaping).

A/B follow Fig. 3k and Methods: FL A = FL-Hap2 (African-derived), B = FL-Hap1 (E. oleifera-derived);
TN A = TK (dura)-like, B = NS (pisifera)-like. Positive log2(A/B) = A-biased.
Outputs: sd/SF12a_genome_vs_module.tsv, sd/SF12a_gene_values.tsv, sd/SF12b_window_DNAcorrected.tsv,
         sd/SF12c_gene_DNA_RNA.tsv, sd/SF12d_key_FA_genes.tsv, work/sf12_stats.json
"""
import json
from pathlib import Path
import numpy as np, pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

H = Path(__file__).resolve().parent; O = H.parent / "enh_A/out"; SD = H / "sd"; SD.mkdir(exist_ok=True)
(H / "work").mkdir(exist_ok=True)
MODS = ["Oil biosynthesis & storage", "TAG assembly & oil body", "De-novo / saturated FA", "Unsaturated FA",
        "Lipid oxidation / antioxidant", "Shell / cell wall / lignin"]
AB = dict(zip(MODS, ["OBS", "TOF", "DSF", "UFA", "LOD", "SCL"]))
WIN = ["Days 0–65", "Days 80–140", "Days 155–185", "Hours 12–72", "All stages"]
WLAB = dict(zip(WIN, ["0–65 d", "80–140 d", "155–185 d", "12–72 h", "All stages"]))
ALLELE = {"FL": ("FL-Hap2 (African-derived)", "FL-Hap1 (E. oleifera-derived)"),
          "TN": ("TK (dura)-like", "NS (pisifera)-like")}

g = pd.read_csv(O / "gene_window_ASE_unified.tsv.gz", sep="\t")
mem = pd.read_csv(O / "module_membership.tsv", sep="\t")
p1 = pd.read_csv(O / "part1_genome_vs_module_unified.tsv", sep="\t")
p2 = pd.read_csv(O / "part2_corrected_genome_vs_module.tsv", sep="\t")
G = pd.read_csv(O / "gene_DNA_ratios.tsv.gz", sep="\t", index_col=0)
k3 = pd.read_csv(O / "part3_key_FA_genes.tsv", sep="\t")
ST = {}

# ---------------- a: genome-wide vs modules, FL and TN, five windows (part1, unified ASE table)
rows = []
for an in ("FL", "TN"):
    for w in WIN:
        x = g[(g.analysis == an) & (g.stage_group == w) & g.robust_any]
        for s in ["Genome-wide"] + MODS:
            if s == "Genome-wide":
                xs = x
            else:
                xs = x[x.gene_id.isin(set(mem[(mem.analysis == an) & (mem.trait_module == s)].gene_id))]
            r = p1[(p1.analysis == an) & (p1.window == w) & (p1.set == s)].iloc[0]
            nb, na = int((xs.med_robust < 0).sum()), int((xs.med_robust > 0).sum())
            assert len(xs) == r.n_robust and abs(100 * nb / len(xs) - r.pct_B) < 1e-3, (an, w, s)
            rows.append(dict(Hybrid=an, A_allele=ALLELE[an][0], B_allele=ALLELE[an][1], Window=WLAB[w],
                             Set=s if s == "Genome-wide" else f"{AB[s]} ({s})", Eligible_genes=int(r.n_eligible),
                             Robust_ASE_genes=len(xs), B_biased_genes=nb, A_biased_genes=na,
                             Pct_B_biased=round(100 * nb / len(xs), 3), Pct_A_biased=round(100 * na / len(xs), 3),
                             Median_log2AB=round(float(r.median_log2AB_robust), 4),
                             Binomial_P_two_sided=float(f"{stats.binomtest(nb, len(xs)).pvalue:.3g}") if s == "Genome-wide" else None,
                             Perm_P_B_one_sided=r.p_perm_lower if s != "Genome-wide" else None,
                             Perm_P_two_sided=r.p_perm_two if s != "Genome-wide" else None))
A = pd.DataFrame(rows)
for an in ("FL", "TN"):
    m = (A.Hybrid == an) & (A.Set != "Genome-wide")
    A.loc[m, "BH_P_B_one_sided"] = multipletests(A.loc[m, "Perm_P_B_one_sided"], method="fdr_bh")[1].round(4)
    A.loc[m, "BH_P_two_sided"] = multipletests(A.loc[m, "Perm_P_two_sided"], method="fdr_bh")[1].round(4)
A.to_csv(SD / "SF12a_genome_vs_module.tsv", sep="\t", index=False)
for an in ("FL", "TN"):
    gw = A[(A.Hybrid == an) & (A.Set == "Genome-wide")].set_index("Window")
    md = A[(A.Hybrid == an) & (A.Set != "Genome-wide")]
    ST[an] = dict(n_robust=int(gw.loc["All stages", "Robust_ASE_genes"]), n_B=int(gw.loc["All stages", "B_biased_genes"]),
                  n_A=int(gw.loc["All stages", "A_biased_genes"]),
                  pct_B=float(gw.loc["All stages", "Pct_B_biased"]), pct_A=float(gw.loc["All stages", "Pct_A_biased"]),
                  binom_P=float(gw.loc["All stages", "Binomial_P_two_sided"]),
                  win_pct_B=[float(gw.loc[WLAB[w], "Pct_B_biased"]) for w in WIN[:4]],
                  win_pct_A=[float(gw.loc[WLAB[w], "Pct_A_biased"]) for w in WIN[:4]],
                  win_binom_P_max=float(gw.loc[[WLAB[w] for w in WIN[:4]], "Binomial_P_two_sided"].max()),
                  n_module_tests=int(len(md)), perm_one_min=float(md.Perm_P_B_one_sided.min()),
                  perm_two_min=float(md.Perm_P_two_sided.min()), bh_one_min=float(md.BH_P_B_one_sided.min()),
                  bh_two_min=float(md.BH_P_two_sided.min()),
                  module_n_all={AB[s]: int(A[(A.Hybrid == an) & (A.Window == "All stages") & (A.Set == f"{AB[s]} ({s})")].Robust_ASE_genes.item()) for s in MODS})

# per-gene values plotted in a (all stages, genes with robust ASE)
x = g[(g.stage_group == "All stages") & g.robust_any][["analysis", "gene_id", "n_elig", "med_robust"]].copy()
for s in MODS:
    x[AB[s]] = [gid in set(mem[(mem.analysis == an) & (mem.trait_module == s)].gene_id) for an, gid in zip(x.analysis, x.gene_id)]
x["Direction"] = np.where(x.med_robust < 0, "B", np.where(x.med_robust > 0, "A", "none"))
x = x.rename(columns={"analysis": "Hybrid", "gene_id": "Gene_ID", "n_elig": "Eligible_stages",
                      "med_robust": "Median_log2AB_robust_stages"}).sort_values(["Hybrid", "Gene_ID"])
x.to_csv(SD / "SF12a_gene_values.tsv", sep="\t", index=False, float_format="%.5g")

# ---------------- b: FL windows, RNA vs RNA - DNA (short-read reciprocal), same gene set
V = {"RNA uncorrected (genes with Short-read reciprocal DNA)": "RNA",
     "RNA - DNA (Short-read reciprocal)": "RNA − DNA"}
B = p2[p2.version.isin(V)].copy()
B["Version"] = B.version.map(V); B["Window"] = B.window.map(WLAB)
B["Set"] = [s if s == "Genome-wide" else f"{AB[s]} ({s})" for s in B.set]
for v in V.values():
    m = (B.Version == v) & (B.Set != "Genome-wide")
    B.loc[m, "BH_P_B_one_sided"] = multipletests(B.loc[m, "p_perm_lower"], method="fdr_bh")[1].round(4)
B = B.rename(columns={"n_eligible": "Eligible_genes", "n_robust": "Robust_ASE_genes", "pct_B": "Pct_B_biased",
                      "median_log2AB": "Median_log2AB", "p_perm_lower": "Perm_P_B_one_sided"})
B["Window"] = pd.Categorical(B.Window, [WLAB[w] for w in WIN], ordered=True)
B = B.sort_values(["Version", "Window"], key=lambda s: s if s.name == "Window" else s.map({"RNA": 0, "RNA − DNA": 1}))
B = B[["Version", "Window", "Set", "Eligible_genes", "Robust_ASE_genes", "Pct_B_biased", "Median_log2AB",
       "Perm_P_B_one_sided", "BH_P_B_one_sided"]]
B.to_csv(SD / "SF12b_window_DNAcorrected.tsv", sep="\t", index=False, float_format="%.5g")
for v in V.values():
    gw = B[(B.Version == v) & (B.Set == "Genome-wide")].set_index("Window")
    md = B[(B.Version == v) & (B.Set != "Genome-wide")]
    ST["b_" + ("rna" if v == "RNA" else "corr")] = dict(
        n_robust_all=int(gw.loc["All stages", "Robust_ASE_genes"]), pct_B_all=float(gw.loc["All stages", "Pct_B_biased"]),
        win_pct_B=[float(gw.loc[WLAB[w], "Pct_B_biased"]) for w in WIN[:4]],
        eligible_all=int(gw.loc["All stages", "Eligible_genes"]),
        perm_one_min=float(md.Perm_P_B_one_sided.min()), bh_one_min=float(md.BH_P_B_one_sided.min()))

# ---------------- c: gene-level DNA allele ratios (short reads; FL-Hap2 and FL-Hap1 references) and RNA
C = G[["nsite_srA", "a_srA", "b_srA", "l2_srA", "nsite_srB", "a_srB", "b_srB", "l2_srB", "l2_sr_sym"]].dropna(subset=["l2_sr_sym"]).copy()
C.columns = ["Sites_on_FL-Hap2_ref", "DNA_A_reads_on_FL-Hap2_ref", "DNA_B_reads_on_FL-Hap2_ref", "DNA_log2AB_FL-Hap2_ref",
             "Sites_on_FL-Hap1_ref", "DNA_A_reads_on_FL-Hap1_ref", "DNA_B_reads_on_FL-Hap1_ref", "DNA_log2AB_FL-Hap1_ref",
             "DNA_log2AB_reciprocal_mean"]
r = g[(g.analysis == "FL") & (g.stage_group == "All stages")].set_index("gene_id")
C["RNA_median_log2AB_all_eligible_stages"] = r.med_all.reindex(C.index)
C.index.name = "Gene_ID"
C.to_csv(SD / "SF12c_gene_DNA_RNA.tsv", sep="\t", float_format="%.5g")
v = C.DNA_log2AB_reciprocal_mean; rr = C.RNA_median_log2AB_all_eligible_stages.dropna()   # RNA on the same genes
ST["c"] = dict(n_dna=int(len(v)), med_dna=float(v.median()), frac_lt05=float((v.abs() < 0.5).mean()),
               med_dnaA=float(C["DNA_log2AB_FL-Hap2_ref"].median()), med_dnaB=float(C["DNA_log2AB_FL-Hap1_ref"].median()),
               n_rna=int(rr.notna().sum()), med_rna=float(rr.median()), frac_rna_lt05=float((rr.abs() < 0.5).mean()))

# ---------------- d: key fatty-acid genes (DNA-corrected), plotted = testable and FL mean TPM >= 5
k = k3.copy()
k["Plotted"] = k.corr_median_log2AB.notna() & (k.FL_TPM >= 5)
k = k.rename(columns={"enzyme": "Enzyme", "gene_A_FLHap2": "Gene_ID_FL-Hap2", "gene_B_FLHap1": "Gene_ID_FL-Hap1",
                      "name": "eggNOG_name", "FL_TPM": "FL_mean_TPM", "n_eligible_stages": "Eligible_stages",
                      "DNA_sr_Aref": "DNA_log2AB_FL-Hap2_ref", "DNA_sr_Bref": "DNA_log2AB_FL-Hap1_ref",
                      "DNA_sr_recip": "DNA_log2AB_reciprocal_mean", "DNA_sites": "DNA_sites",
                      "RNA_median_log2AB": "RNA_median_log2AB", "corr_median_log2AB": "DNA_corrected_median_log2AB",
                      "stages_robustA_corr": "Stages_robust_A_corrected", "stages_robustB_corr": "Stages_robust_B_corrected",
                      "stages_robustA_orig": "Stages_robust_A_uncorrected", "stages_robustB_orig": "Stages_robust_B_uncorrected",
                      "verdict_corrected": "Call_after_DNA_correction"})
k["Call_after_DNA_correction"] = k.Call_after_DNA_correction.replace(
    {"E. oleifera (B) allele higher": "FL-Hap1 (B) higher", "African (A) allele higher": "FL-Hap2 (A) higher"})
k = k[["Enzyme", "Gene_ID_FL-Hap2", "Gene_ID_FL-Hap1", "eggNOG_name", "FL_mean_TPM", "Eligible_stages", "DNA_sites",
       "DNA_log2AB_FL-Hap2_ref", "DNA_log2AB_FL-Hap1_ref", "DNA_log2AB_reciprocal_mean", "RNA_median_log2AB", "DNA_corrected_median_log2AB",
       "Stages_robust_A_uncorrected", "Stages_robust_B_uncorrected", "Stages_robust_A_corrected",
       "Stages_robust_B_corrected", "Call_after_DNA_correction", "Plotted"]]
k.to_csv(SD / "SF12d_key_FA_genes.tsv", sep="\t", index=False, float_format="%.3f")
kp = k[k.Plotted]
ST["d"] = dict(n_plotted=int(len(kp)), calls=kp.Call_after_DNA_correction.value_counts().to_dict(),
               max_abs_dna=float(kp.DNA_log2AB_reciprocal_mean.abs().max()))
json.dump(ST, open(H / "work/sf12_stats.json", "w"), indent=1, ensure_ascii=False)
print(json.dumps(ST, indent=1, ensure_ascii=False))
