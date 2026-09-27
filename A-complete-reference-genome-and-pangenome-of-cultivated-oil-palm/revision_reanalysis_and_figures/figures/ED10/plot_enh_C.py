#!/usr/bin/env python3
import os, pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
H = os.path.dirname(os.path.abspath(__file__))
plt.rcParams.update({"font.family": "Arial", "font.size": 6, "axes.linewidth": 0.5, "xtick.major.width": 0.5,
                     "ytick.major.width": 0.5, "xtick.major.size": 2, "ytick.major.size": 2, "pdf.fonttype": 42,
                     "axes.spines.top": False, "axes.spines.right": False, "legend.frameon": False})
C = pd.read_csv(f"{H}/constrained_paths.tsv", sep="\t")
S = pd.read_csv(f"{H}/weight_sweep.tsv", sep="\t")
PS = pd.read_csv(f"{H}/data/penalty_sweep.tsv", sep="\t")
BSL = int(C.Best_Single_Donor_Load.iloc[0])
BLUE, ORANGE, GREY, DARK = "#1f5fa8", "#d9751e", "#9a9a9a", "#222222"

fig, (a, b) = plt.subplots(1, 2, figsize=(180 / 25.4, 68 / 25.4), gridspec_kw=dict(wspace=0.32))
# ---- panel a
A = C[C.Constraint_Class == "a_per_chrom_breakpoints"]
B = C[C.Constraint_Class.isin(["b_donor_budget"])]
ref = C[C.Constraint_Class == "reference"]
single = ref.iloc[0]; full = ref.iloc[1]
Bx = pd.concat([single.to_frame().T, B])
ps = PS[PS.Breakpoint_Penalty_P > 15].copy()
ps["pct"] = 100 * (BSL - ps.Residual_Total_Load) / BSL
a.plot(ps.Donor_Count, ps.pct, "o--", color=GREY, ms=2.5, lw=0.6, label="Penalised DP, P = 40/100/300")
for r in ps.itertuples():
    a.annotate(f"P={r.Breakpoint_Penalty_P} ({r.Breakpoint_Count} bp)", (r.Donor_Count, r.pct), xytext=(-4, 0),
               textcoords="offset points", ha="right", va="center", fontsize=5, color=GREY)
a.plot(Bx.Donor_Count.astype(int), Bx.Reduction_vs_Best_Single_pct.astype(float), "s-", color=ORANGE, ms=3, lw=0.8,
       label="≤ m donors genome-wide, ≤ 2 bp/chromosome")
for r in Bx.itertuples():
    lab = "single\nEG_107" if r.Donor_Count == 1 else f"m={r.Donor_Count}"
    a.annotate(lab, (r.Donor_Count, r.Reduction_vs_Best_Single_pct), xytext=(0, 4), textcoords="offset points",
               ha="center", fontsize=5, color=ORANGE)
a.plot(A.Donor_Count, A.Reduction_vs_Best_Single_pct, "o-", color=BLUE, ms=3, lw=0.8,
       label="≤ k bp/chromosome, donors unrestricted")
for r in A.itertuples():
    k = r.Max_Breakpoints_per_Chrom
    a.annotate(f"k={k} ({r.Total_Breakpoints} bp)", (r.Donor_Count, r.Reduction_vs_Best_Single_pct), xytext=(5, -1), textcoords="offset points",
               ha="left", va="center", fontsize=5, color=BLUE)
a.plot([full.Donor_Count], [full.Reduction_vs_Best_Single_pct], "*", color=DARK, ms=6,
       label="Unconstrained IPH (P = 15)")
a.annotate(f"{full.Reduction_vs_Best_Single_pct:.1f}%\n{full.Total_Breakpoints} bp", (full.Donor_Count, full.Reduction_vs_Best_Single_pct),
           xytext=(-5, -2), textcoords="offset points", ha="right", va="top", fontsize=5)
a.set_xscale("log"); a.set_xticks([1, 2, 3, 5, 10, 20, 35]); a.set_xticklabels(["1", "2", "3", "5", "10", "20", "35"])
a.set_xlim(0.7, 80); a.set_ylim(-3, 68)
a.set_xlabel("Donor haplotypes used in the design"); a.set_ylabel("Load reduction vs best single donor (%)")
a.legend(loc="upper left", fontsize=5, handlelength=2.2)
a.set_title("a  Feasibility-constrained designs (African35)", loc="left", fontsize=7, fontweight="bold")

# ---- panel b
for lab, col, ls, name in (("step9_calls_only", BLUE, "-", "Reward from step9 donor calls (29/35 donors)"),
                           ("collapsed_proxy_for_6_phased_donors", ORANGE, "--", "Collapsed-proxy calls for 6 phased donors")):
    g = S[S.Reward_Call_Set == lab].sort_values("GWAS_Weight_W")
    x = g.Fav_Captured if lab == "step9_calls_only" else g.Fav_Captured_collapsed_proxy
    b.plot(x, g.Residual_Total_Load, ls, marker="o", color=col, ms=3, lw=0.8, label=name)
    for xx, r in zip(x, g.itertuples()):
        off = {("step9_calls_only", 4): (0, -9, "center"), ("collapsed_proxy_for_6_phased_donors", 2): (-6, 3, "right"), ("step9_calls_only", 8): (-6, 0, "right")}.get(
            (lab, r.GWAS_Weight_W), (0, 5, "center") if lab == "step9_calls_only" else (0, -8, "center"))
        if lab != "step9_calls_only" and r.GWAS_Weight_W == 8: off = (6, 0, "left")
        b.annotate(f"W={r.GWAS_Weight_W}", (xx, r.Residual_Total_Load), xytext=off[:2], textcoords="offset points",
                   ha=off[2], va="center" if off[1] == 0 else "baseline", fontsize=5, color=col)
b.set_xlabel("Favourable GWAS alleles captured (of 284)"); b.set_ylabel("Residual load (dSV + dSNP)")
b.set_ylim(23505, 23590); b.set_xlim(163, 230)
b.legend(loc="upper left", fontsize=5)
b.set_title("b  GWAS-weight sweep, P = 15", loc="left", fontsize=7, fontweight="bold")
for ext in ("pdf", "png"):
    fig.savefig(f"{H}/EnhC_constrained_and_weight_sweep.{ext}", dpi=600, bbox_inches="tight")
print("ok")
