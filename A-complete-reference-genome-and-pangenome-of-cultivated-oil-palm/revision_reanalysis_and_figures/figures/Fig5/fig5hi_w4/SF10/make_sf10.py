#!/usr/bin/env python3
"""Supplementary Fig. 10: feasibility-constrained donor designs and the GWAS-weight sweep (African35).

Inputs: ../../enh_C/constrained_paths.tsv, ../../enh_C/weight_sweep.tsv (enh_C.py), ../by_chrom_W4.tsv,
        ../per_chrom_capture_W4.tsv (the Fig. 5h path, W = 4, P = 15).
Outputs: Supplementary_Fig_10.pdf (vector, Arial, 180 mm wide), Supplementary_Fig_10.png (600 dpi),
         ../sd/SF10a_constrained_designs.tsv, ../sd/SF10b_weight_sweep.tsv (Source Data sheets).
Style follows the existing Supplementary Figs (bold lowercase panel letters, Arial 6-7 pt, no panel titles).
"""
import re
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

H = Path(__file__).resolve().parent; W4D = H.parent; FIX = W4D.parent
C = pd.read_csv(FIX / "enh_C/constrained_paths.tsv", sep="\t")
S = pd.read_csv(FIX / "enh_C/weight_sweep.tsv", sep="\t")
S = S[S.Reward_Call_Set == "step9_calls_only"].sort_values("GWAS_Weight_W").reset_index(drop=True)
bc = pd.read_csv(W4D / "by_chrom_W4.tsv", sep="\t").set_index("Chrom").loc["TOTAL"]
cap4 = int(pd.read_csv(W4D / "per_chrom_capture_W4.tsv", sep="\t").set_index("Chrom").loc["TOTAL", "Captured_W4"])
BSL = int(C.Best_Single_Donor_Load.iloc[0])
w4 = S[S.GWAS_Weight_W == 4].iloc[0]
# the Fig. 5h path (gen_path_w4.py) and the sweep row must agree
assert (int(bc.Residual_Total_Load), int(bc.Breakpoint_Count), cap4) == (int(w4.Residual_Total_Load), int(w4.Breakpoints), int(w4.Fav_Captured))
W4_PCT = 100 * (BSL - int(bc.Residual_Total_Load)) / BSL


def disp(d):
    m = re.fullmatch(r"(dura|pisifera|nrly|bk)_hap([12])", d)
    return f"{dict(dura='TK', pisifera='NS', nrly='Nigerian', bk='TN')[m.group(1)]}-Hap{m.group(2)}" if m else d


plt.rcParams.update({"font.family": "Arial", "font.size": 6.5, "axes.linewidth": 0.6, "xtick.major.width": 0.6,
                     "ytick.major.width": 0.6, "xtick.minor.width": 0.4, "xtick.major.size": 2.5, "ytick.major.size": 2.5,
                     "xtick.minor.size": 1.5, "axes.labelsize": 7, "xtick.labelsize": 6.5, "ytick.labelsize": 6.5,
                     "pdf.fonttype": 42, "ps.fonttype": 42, "axes.spines.top": False, "axes.spines.right": False,
                     "legend.frameon": False, "legend.fontsize": 6})
BLUE, ORANGE, DARK, RED = "#3A76AF", "#E39B1B", "#222222", "#C0392B"

fig, (a, b) = plt.subplots(1, 2, figsize=(180 / 25.4, 72 / 25.4), gridspec_kw=dict(wspace=0.34, left=0.085, right=0.985,
                                                                                     bottom=0.16, top=0.92))
# ---------------- a
A = C[C.Constraint_Class == "a_per_chrom_breakpoints"].copy()
B = C[C.Constraint_Class == "b_donor_budget"].copy()
single = C[C.Design.str.startswith("Best single")].iloc[0]
Bx = pd.concat([single.to_frame().T, B])
a.plot(Bx.Donor_Count.astype(int), Bx.Reduction_vs_Best_Single_pct.astype(float), "s-", color=ORANGE, ms=3.2, lw=0.9,
       mec="white", mew=0.3, label="≤ m donors genome-wide, ≤ 2 breakpoints per chromosome")
for r in Bx.itertuples():
    if r.Donor_Count == 1:
        a.annotate("Best single donor (EG_107)", (1, 0), xytext=(4, -3), textcoords="offset points",
                   ha="left", va="top", fontsize=5.8, color=ORANGE)
    else:
        a.annotate(f"m = {r.Donor_Count}", (r.Donor_Count, r.Reduction_vs_Best_Single_pct), xytext=(3, -4),
                   textcoords="offset points", ha="left", va="top", fontsize=5.8, color=ORANGE)
a.plot(A.Donor_Count, A.Reduction_vs_Best_Single_pct, "o-", color=BLUE, ms=3.2, lw=0.9, mec="white", mew=0.3,
       label="≤ k breakpoints per chromosome, donors unrestricted")
for r in A.itertuples():
    a.annotate(f"k = {r.Max_Breakpoints_per_Chrom} ({r.Total_Breakpoints})", (r.Donor_Count, r.Reduction_vs_Best_Single_pct),
               xytext=(-5, 0), textcoords="offset points", ha="right", va="center", fontsize=5.8, color=BLUE)
a.plot([int(bc.Donor_Count)], [W4_PCT], "*", color=DARK, ms=7, mec="white", mew=0.3,
       label="Unconstrained path of Fig. 5h (W = 4, P = 15)")
a.annotate(f"{W4_PCT:.1f}% ({int(bc.Breakpoint_Count)})", (int(bc.Donor_Count), W4_PCT), xytext=(0, -6),
           textcoords="offset points", ha="center", va="top", fontsize=5.8, color=DARK)
a.set_xscale("log"); a.set_xticks([1, 2, 3, 5, 10, 20, 35]); a.set_xticklabels(["1", "2", "3", "5", "10", "20", "35"])
a.set_xlim(0.75, 48); a.set_ylim(-6, 78); a.set_yticks(range(0, 70, 10))
a.set_xlabel("Donor haplotypes in the design"); a.set_ylabel("Reduction in candidate deleterious burden\nrelative to the best single donor (%)")
a.legend(loc="upper left", handlelength=2.2, borderaxespad=0.2, fontsize=5.8)

# ---------------- b
b.plot(S.Fav_Captured, S.Residual_Total_Load, "o-", color=BLUE, ms=3.4, lw=0.9, mec="white", mew=0.3, zorder=2)
for r in S.itertuples():
    off = {0: (4, -6, "left"), 1: (0, 6, "center"), 2: (-5, 5, "right"), 4: (6, -7, "left"), 8: (-6, 0, "right")}[r.GWAS_Weight_W]
    a_ = b.annotate(f"W = {r.GWAS_Weight_W}\n{r.Breakpoints} breakpoints" if r.GWAS_Weight_W in (0, 4) else f"W = {r.GWAS_Weight_W}",
                    (r.Fav_Captured, r.Residual_Total_Load), xytext=off[:2], textcoords="offset points", ha=off[2],
                    va="center" if off[1] == 0 else ("bottom" if off[1] > 0 else "top"), fontsize=5.8,
                    color=RED if r.GWAS_Weight_W == 4 else BLUE, fontweight="bold" if r.GWAS_Weight_W == 4 else "normal")
b.plot([w4.Fav_Captured], [w4.Residual_Total_Load], "o", color=RED, ms=5.2, mec="white", mew=0.4, zorder=3)
b.set_xlim(160, 222); b.set_ylim(23505, 23590)
b.set_xlabel("SV-GWAS favourable loci captured (of 284)"); b.set_ylabel("Residual candidate deleterious burden\n(dSV + dSNP)")
for ax, l in ((a, "a"), (b, "b")):
    fig.text(0.005 if l == "a" else ax.get_position().x0 - 0.075, 0.975, l, fontsize=9, fontweight="bold", va="top")
fig.savefig(H / "Supplementary_Fig_10.pdf")
fig.savefig(H / "Supplementary_Fig_10.png", dpi=600)
plt.close(fig)

# ---------------- Source Data
SD = W4D / "sd"; SD.mkdir(exist_ok=True)
rowsA = []
for r in pd.concat([single.to_frame().T, A, B]).itertuples():
    rowsA.append(dict(Design=r.Design, Constraint=("reference" if r.Constraint_Class == "reference" else
                      ("k breakpoints per chromosome" if r.Constraint_Class.startswith("a_") else "m donors genome-wide, <= 2 breakpoints per chromosome")),
                      Donor_Haplotypes_Used=int(r.Donor_Count), Total_Breakpoints=int(r.Total_Breakpoints),
                      Max_Breakpoints_per_Chrom=int(r.Max_Breakpoints_per_Chrom), Residual_DSV=int(r.Residual_DSV),
                      Residual_DSNP=int(r.Residual_DSNP), Residual_Total_Load=int(r.Residual_Total_Load),
                      Best_Single_Donor_Load=BSL, Reduction_vs_Best_Single_pct=float(r.Reduction_vs_Best_Single_pct),
                      Donor_Set=("" if r.Constraint_Class != "b_donor_budget" else
                                 ",".join(disp(x) for x in re.search(r"Optimal donor set: ([^;]+)", r.Note).group(1).split(","))),
                      Plotted_As=("orange square" if r.Constraint_Class != "a_per_chrom_breakpoints" else "blue circle")))
rowsA.append(dict(Design="Unconstrained path, Fig. 5h (W = 4, P = 15)", Constraint="none (breakpoint penalty P = 15)",
                  Donor_Haplotypes_Used=int(bc.Donor_Count), Total_Breakpoints=int(bc.Breakpoint_Count),
                  Max_Breakpoints_per_Chrom="", Residual_DSV=int(bc.Residual_DSV), Residual_DSNP=int(bc.Residual_DSNP),
                  Residual_Total_Load=int(bc.Residual_Total_Load), Best_Single_Donor_Load=BSL,
                  Reduction_vs_Best_Single_pct=round(W4_PCT, 2), Donor_Set="", Plotted_As="black star"))
pd.DataFrame(rowsA).to_csv(SD / "SF10a_constrained_designs.tsv", sep="\t", index=False)
Sb = S[["GWAS_Weight_W", "Breakpoint_Penalty_P", "Residual_DSV", "Residual_DSNP", "Residual_Total_Load", "Breakpoints",
        "Segments", "Donor_Count", "Objective", "Reduction_vs_Best_Single_pct", "Fav_Captured", "Fav_Total"]].copy()
Sb["Plotted_As"] = ["red circle (Fig. 5h path)" if w == 4 else "blue circle" for w in Sb.GWAS_Weight_W]
Sb.to_csv(SD / "SF10b_weight_sweep.tsv", sep="\t", index=False)
print(pd.DataFrame(rowsA)[["Design", "Donor_Haplotypes_Used", "Total_Breakpoints", "Reduction_vs_Best_Single_pct", "Donor_Set"]].to_string())
print(Sb.to_string())
