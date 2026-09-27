#!/usr/bin/env python3
"""Extended Data Fig. 6 redraw (nut-weight SNP/SV GWAS at the AGL11/SHELL locus, robust n = 132 analysis).

Inputs:
  ../ED4/sd.pkl  Source Data sheets ED7a_joint_GWAS (all 370,706 SVs; SNPs every 200th + all P < 1e-5),
                 ED7d_variants, SF16ab_QQ (SNP bulk quantiles, every 250th rank)
  data/ED6c_regional_3.05-3.41Mb.tsv  regional SNP/SV statistics for chr01B:3,050,000-3,410,000 (6,701 SNPs, 72 SVs;
                                      robust n = 132 EMMAX .ps; r2 to lead SNP from GP_maf_allChr, 308 samples,
                                      same estimator as the original AGL11 build; built by ../ed6c_build_regional.py)
  src/RECOMMENDED_robust_n132_panel_h_SV_genotype_Nut_weight.tsv      genotype x phenotype (n = 132)
  src/AGL11_gene_model.gtf (gene_id evm.TU.chr01B.166, 7 exons, minus strand)
  src/SHARED_panel_h_field_shell_thickness.tsv / _field_inference.tsv  field subset (tree means)
  src/xys果实横切面_安全裁切_不含比例尺.png                              photograph (unchanged crop)
"""
import pickle
from pathlib import Path
from statistics import NormalDist

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from PIL import Image

HERE = Path(r"${WORK_DIR}/fix/redraw2/ED/work/ED6")   # restructure: original location (inputs read-only)
BASE = HERE.parents[3] / "beautify/work/ED6"          # redraw2: inputs read from the beautify build (read only)
OUT, DATA, SRC = Path(r"${WORK_DIR}/fix/restructure/ed6_no_fg"), Path(r"${WORK_DIR}/fix/restructure/ed6_no_fg") / "data", BASE / "src"   # restructure: outputs elsewhere
IN_DATA = BASE / "data"
OUT.mkdir(exist_ok=True)
DATA.mkdir(exist_ok=True)
MM = 1 / 25.4
DARK, GREY = "#242A30", "#B9C0C7"
SNPC, SVC = "#3C78A8", "#D9772B"
import sys; sys.path.insert(0, str(HERE.parents[1] / "common")); import palA
import matplotlib._mathtext as _mt; _mt.SHRINK_FACTOR = 0.8  # beautify: super/subscripts >= 4 pt (default 0.7)
TYPE_COL = {"INS": palA.SV["INS"], "DEL": palA.SV["DEL"], "MNV/COMPLEX": palA.SV["COMPLEX"], "DUP": palA.SV["DUP"],
            "INV": palA.SV["INV"]}  # beautify: scheme A SV colours (= ED5, SF9, Fig. 5)
import json as _json   # [snp repair 2026-09-27] repaired SNP call set: number of SNP tests and lambda from ed9_meta.json
_META = _json.load(open("${WORK_DIR}/fix/snp_repair_rerun/ed9/ed9_meta.json"))
SNP_BONF, SV_BONF = 0.05 / _META["n_snp_tests"], 0.05 / 370706
F1, F2 = "SV_chr01B_3231689_97849432", "SV_chr01B_3261344_43f7146f"
FOCAL_LAB = {F1: "1.76-kb complex (MNV) record", F2: "17-bp deletion"}

mpl.rcParams.update({
    "font.family": "Arial", "font.size": 7, "axes.labelsize": 7, "axes.titlesize": 7,
    "xtick.labelsize": 6, "ytick.labelsize": 6, "legend.fontsize": 6,
    "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6,
    "xtick.major.size": 2.5, "ytick.major.size": 2.5, "axes.unicode_minus": True,
    "mathtext.fontset": "custom", "mathtext.rm": "Arial", "mathtext.it": "Arial:italic",
    "mathtext.bf": "Arial:bold", "pdf.fonttype": 42, "svg.fonttype": "none",
})


def ptex(p, nd=1):
    m, e = f"{p:.{nd}e}".split("e")
    return rf"$P$ = {m} × 10$^{{{int(e)}}}$"


def clean(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def letter(fig, ax, s, dx=-0.05, dy=0.008):
    bb = ax.get_position()
    fig.text(bb.x0 + dx, bb.y1 + dy, s, fontsize=palA.LETTER_PT, fontweight="bold", va="bottom", ha="left")


d = pickle.load(open(BASE.parent / "ED4/sd.pkl", "rb"))
gw = d["ED7a_joint_GWAS"]
assert (gw.analysis_version == "robust_n132").all()
reg = pd.read_csv(IN_DATA / "ED6c_regional_3.05-3.41Mb.tsv", sep="\t")
# [snp repair 2026-09-27] SNP scan, QQ and regional SNPs/r2 from the repaired SNP call set (fix/snp_repair_rerun/ed9)
NEW = Path("${WORK_DIR}/fix/snp_repair_rerun/ed9")
off = (gw.cumulative_pos - gw.pos).groupby(gw.chrom).first()
ns = pd.read_csv(NEW / "ED9a_SNP_display.tsv", sep="	"); ns["p"] = 10 ** -ns.neglog10_p
ns["cumulative_pos"] = ns.pos + ns.chrom.map(off); ns["marker"] = ns.chrom + ":" + ns.pos.astype(str); ns["analysis_version"] = "robust_n132"
gw = pd.concat([gw[gw.variant_type == "SV"], ns], ignore_index=True)
nr = pd.read_csv(NEW / "ED9c_regional_SNPs.tsv", sep="	")
svr2 = pd.read_csv(NEW / "regional_sv_r2.tsv", sep="	", header=None, names=["marker", "pos", "r2", "n"]).set_index("marker").r2
rsv = reg[reg.variant_type == "SV"].copy(); rsv["r2_to_lead"] = rsv.marker.map(svr2).fillna(rsv.r2_to_lead)
reg = pd.concat([rsv, nr], ignore_index=True)
assert (reg.variant_type == "SNP").sum() == len(nr) and (reg.variant_type == "SV").sum() == 72
geno = pd.read_csv(SRC / "RECOMMENDED_robust_n132_panel_h_SV_genotype_Nut_weight.tsv", sep="\t")
field = pd.read_csv(SRC / "SHARED_panel_h_field_shell_thickness.tsv", sep="\t")
finf = pd.read_csv(SRC / "SHARED_panel_h_field_inference.tsv", sep="\t").iloc[0]
var = d["ED7d_variants"].copy()
gtf = pd.read_csv(SRC / "AGL11_gene_model.gtf", sep="\t", header=None, comment="#")
ex = gtf[(gtf[2] == "exon") & gtf[8].str.contains('gene_id "evm.TU.chr01B.166";', regex=False)]
assert len(ex) == 7

snp = gw[gw.variant_type == "SNP"]
sv = gw[gw.variant_type == "SV"]
assert len(sv) == 370706

fig = plt.figure(figsize=(180 * MM, 170 * MM))

# ---------------------------------------------------------------- a Manhattan (redraw2: mirrored SNP/SV layout)
# SNP scan plotted upwards, SV scan downwards on the same genomic axis; chromosomes alternate dark/light shades
# of the modality colour. Same points as before: SNP display subset (every 200th + all P < 1e-5), all SVs.
ax = fig.add_axes([0.07, 0.745, 0.915, 0.225])
chroms = sorted(gw.chrom.unique(), key=lambda c: int(c[3:5]))
bounds = gw.groupby("chrom").cumulative_pos.agg(["min", "max"]).loc[chroms]
starts = bounds["min"] - gw.groupby("chrom").pos.min().loc[chroms]
SNP_LT, SV_LT = palA.tint(SNPC, 0.5), palA.tint(SVC, 0.5)
ci = {c: i for i, c in enumerate(chroms)}
snp_col = np.where(snp.chrom.map(ci).values % 2 == 0, SNPC, SNP_LT)
sv_col = np.where(sv.chrom.map(ci).values % 2 == 0, SVC, SV_LT)
ax.scatter(snp.cumulative_pos, snp.neglog10_p, s=1.2, c=snp_col, lw=0, rasterized=True, zorder=3)
ax.scatter(sv.cumulative_pos, -sv.neglog10_p, s=1.6, marker="^", c=sv_col, lw=0, rasterized=True, zorder=3)
ax.axhline(0, color=DARK, lw=0.6, zorder=4)
ax.axhline(-np.log10(SNP_BONF), color=SNPC, ls="--", lw=0.6, zorder=2)
ax.axhline(np.log10(SV_BONF), color=SVC, ls="--", lw=0.6, zorder=2)
bb = dict(facecolor="white", edgecolor="none", pad=0.5, alpha=0.9)
ax.text(0.998, -np.log10(SNP_BONF) + 0.3, "SNP Bonferroni, " + ptex(SNP_BONF, 2),
        transform=ax.get_yaxis_transform(), ha="right", va="bottom", fontsize=5.5, color=SNPC, bbox=bb, zorder=5)
ax.text(0.998, np.log10(SV_BONF) - 0.35, "SV Bonferroni, " + ptex(SV_BONF, 2),
        transform=ax.get_yaxis_transform(), ha="right", va="top", fontsize=5.5, color=SVC, bbox=bb, zorder=5)
ax.text(0.33, 0.965, "SNP (displayed subset)", transform=ax.transAxes, ha="left", va="top", fontsize=6,
        color=SNPC)
ax.text(0.33, 0.035, "SV (all tested)", transform=ax.transAxes, ha="left", va="bottom", fontsize=6, color=SVC)
lead = snp.loc[snp.p.idxmin()]
ax.annotate(r"$\it{AGL11}$/$\it{SHELL}$ region" + "\n(chr01B, ~3.15–3.26 Mb)", xy=(lead.cumulative_pos, lead.neglog10_p),
            xytext=(lead.cumulative_pos + 1.6e8, 12.0), fontsize=5.5, va="center",
            arrowprops=dict(arrowstyle="-", lw=0.5, color="#555555"))
ax.set_xticks([(starts[c] + bounds.loc[c, "max"]) / 2 for c in chroms], [str(int(c[3:5])) for c in chroms])
ax.tick_params(axis="x", length=0, pad=2)
ax.set_xlim(0, bounds["max"].max())
ax.set_ylim(-10.2, 14.5)
yt = [-7.5, -5, -2.5, 0, 2.5, 5, 7.5, 10, 12.5]
ax.set_yticks(yt, [f"{abs(v):g}" for v in yt])
ax.set_xlabel("Chromosome (FL-Hap2 coordinates)", labelpad=2)
ax.set_ylabel(r"$-$log$_{10}$ $P$")
clean(ax)
ax.spines["bottom"].set_visible(False)
ax_a = ax

# ---------------------------------------------------------------- b QQ (exact tails)
ax = fig.add_axes([0.07, 0.405, 0.17, 0.26])
N_SNP, N_SV = _META["n_snp_tests"], 370706
# SV: all tests available -> exact
p_sv = np.sort(sv.p.values)
e_sv = -np.log10(np.arange(1, N_SV + 1) / N_SV)
o_sv = -np.log10(p_sv)
# SNP: all P < 1e-5 are present in ED7a -> exact ranks 1..m; bulk from SF16ab_QQ (every 250th rank, same i/N convention)
nq = pd.read_csv(NEW / "ED9b_SNP_QQ.tsv", sep="	")
e_snp, o_snp = nq.expected_neglog10P.values, nq.observed_neglog10P.values
m = int((nq.rank_source == "exact_rank_all_P_lt_1e-5").sum()); bulk = nq[nq.rank_source != "exact_rank_all_P_lt_1e-5"]
lam_sv = NormalDist().inv_cdf(1 - np.median(p_sv) / 2) ** 2 / 0.454936
# SNP lambda_GC from all 28,192,651 P values with the same median-P formula as the SV value
# (median chi2(1) quantile of P / 0.454936); computed on the cluster from the 16 emmax_chr*.ps files
# (work/trace/ED/out/ed6_gwas.log: 0.9539857335928091). SV value is recomputed below from all 370,706 P.
lam_snp = _META["lambda_snp"]  # [snp repair] published model, nut weight
thin = np.r_[np.arange(0, 2000), np.arange(2000, N_SV, 20)]
ax.scatter(e_snp, o_snp, s=2, color=SNPC, lw=0, rasterized=True, label=rf"SNP ($\lambda_{{GC}}$ = {lam_snp:.3f})")
ax.scatter(e_sv[thin], o_sv[thin], s=2.5, marker="^", color=SVC, lw=0, rasterized=True,
           label=rf"SV ($\lambda_{{GC}}$ = {lam_sv:.3f})")
lim = 14
ax.plot([0, 7.6], [0, 7.6], color="#777777", lw=0.6, ls="--")
ax.set_xlim(0, 7.8)
ax.set_ylim(0, lim)
ax.set_xlabel(r"Expected $-$log$_{10}$ $P$")
ax.set_ylabel(r"Observed $-$log$_{10}$ $P$")
ax.legend(frameon=False, loc="upper left", handletextpad=0.1, markerscale=2.5)
clean(ax)
ax_b = ax
pd.DataFrame({"variant_type": ["SNP"] * len(e_snp) + ["SV"] * len(thin), "expected_neglog10P": np.r_[e_snp, e_sv[thin]],
              "observed_neglog10P": np.r_[o_snp, o_sv[thin]],
              "rank_source": ["exact_rank_all_P_lt_1e-5"] * m + ["every_250th_rank"] * len(bulk) +
                             ["exact_rank_all_370706_SV_display_thinned"] * len(thin)}).to_csv(
    DATA / "ED6b_nut_weight_QQ.tsv", sep="\t", index=False)

# ---------------------------------------------------------------- c regional + gene
X0, X1 = 3.05, 3.41
ax = fig.add_axes([0.33, 0.515, 0.655, 0.15])
rs = reg[reg.variant_type == "SNP"]
cmap = mpl.colormaps["viridis_r"]
sc = ax.scatter(rs.pos / 1e6, rs.neglog10_p, c=rs.r2_to_lead, cmap=cmap, vmin=0, vmax=1, s=2.2, lw=0,
                rasterized=True, zorder=2)
rv = reg[reg.variant_type == "SV"]
for t, mk in (("INS", "^"), ("DEL", "v"), ("MNV/COMPLEX", "D")):
    s = rv[rv.svtype == t]
    ax.scatter(s.pos / 1e6, s.neglog10_p, marker=mk, s=9, facecolor="white", edgecolor=TYPE_COL[t], lw=0.7,
               zorder=3, label=f"SV: {t.lower() if t != 'MNV/COMPLEX' else 'complex (MNV)'}")
for f in (F1, F2):
    r = rv[rv.marker == f].iloc[0]
    ax.scatter(r.pos / 1e6, r.neglog10_p, marker="D" if f == F1 else "v", s=22, color=TYPE_COL[r.svtype],
               edgecolor="black", lw=0.5, zorder=4)
    ax.annotate(f"{FOCAL_LAB[f]}\n" + ptex(r.p) + f"; $r^2$ to lead = {r.r2_to_lead:.2f}",
                xy=(r.pos / 1e6, r.neglog10_p), xytext=(3.172 if f == F1 else r.pos / 1e6 + 0.004, 11.6 if f == F1 else 9.6),
                ha="left", fontsize=5.2,
                arrowprops=dict(arrowstyle="-", lw=0.45, color="#555555"))
L = rs.loc[rs.p.idxmin()]
ax.scatter(L.pos / 1e6, L.neglog10_p, marker="*", s=45, color="#C0392B", edgecolor="black", lw=0.4, zorder=5)
ax.text(L.pos / 1e6 - 0.004, L.neglog10_p - 0.9, f"lead SNP chr01B:{int(L.pos):,}\n" + ptex(L.p), fontsize=5.2,
        va="center", ha="right")  # redraw2: moved below-left of the star, away from the colour bar
ax.axhline(-np.log10(SV_BONF), color=SVC, ls="--", lw=0.6)
ax.axhline(-np.log10(SNP_BONF), color=SNPC, ls=":", lw=0.7)  # extfix: SNP genome-wide threshold

ax.set_xlim(X0, X1)
ax.set_ylim(0, 14.5)
ax.set_ylabel(r"$-$log$_{10}$ $P$")
ax.tick_params(labelbottom=False)
hl, lb = ax.get_legend_handles_labels()
hl.append(mpl.lines.Line2D([], [], color=SVC, ls="--", lw=0.6))
lb.append("SV threshold")
hl.append(mpl.lines.Line2D([], [], color=SNPC, ls=":", lw=0.7))  # extfix
lb.append("SNP threshold")
ax.legend(hl, lb, frameon=False, loc="lower right", ncol=5, handletextpad=0.3, columnspacing=0.8,
          bbox_to_anchor=(1.0, 0.985), handlelength=1.6, borderaxespad=0.0)
clean(ax)
cax = fig.add_axes([0.435, 0.672, 0.07, 0.008])  # redraw2: colour bar moved above the plotting area
cb = fig.colorbar(sc, cax=cax, orientation="horizontal", ticks=[0, 0.5, 1])
cb.ax.tick_params(labelsize=5, length=1.5, pad=1)
cax.xaxis.set_ticks_position("top")
cb.outline.set_linewidth(0.4)
cax.text(-0.06, 0.5, r"SNP $r^2$ to lead", transform=cax.transAxes, ha="right", va="center", fontsize=5.5)
cax.set_zorder(2)
ax_c = ax

axg = fig.add_axes([0.33, 0.478, 0.655, 0.03])
g0, g1 = ex[3].min() / 1e6, ex[4].max() / 1e6
axg.plot([g0, g1], [0, 0], color="#2A9D8F", lw=0.8)
for _, e in ex.iterrows():
    axg.add_patch(plt.Rectangle((e[3] / 1e6, -0.45), (e[4] - e[3]) / 1e6, 0.9, color="#2A9D8F", lw=0))
for xx in np.linspace(g0 + 0.002, g1 - 0.002, 5):
    axg.annotate("", xy=(xx - 0.0015, 0), xytext=(xx + 0.0015, 0),
                 arrowprops=dict(arrowstyle="-|>", lw=0.5, color="#2A9D8F", mutation_scale=4))
axg.text(g1 + 0.002, 0, r"$\it{AGL11}$/STK-like (evm.TU.chr01B.166, − strand)", fontsize=5.3, va="center",
         color="#1E6F66")
axg.set_xlim(X0, X1)
axg.set_ylim(-1, 1)
axg.axis("off")

# ---------------------------------------------------------------- d graph SV alleles in the locus
ax = fig.add_axes([0.33, 0.405, 0.655, 0.065])
var["mid"] = var.pos / 1e6
var = var.merge(gw[gw.variant_type == "SV"][["marker", "p"]], left_on="sv_id", right_on="marker", how="left")
for _, r in var.iterrows():
    col = TYPE_COL.get(r.svtype, "#999999")
    af = r.AF if pd.notna(r.AF) else np.nan
    focal = r.sv_id in (F1, F2)
    if pd.notna(af):
        ax.plot([r.mid, r.mid], [0, af], color=col, lw=1.3 if focal else 0.9, solid_capstyle="butt")
        ax.scatter(r.mid, af, s=10 if focal else 5, color=col, edgecolor="black" if focal else col, lw=0.5, zorder=3)
    else:
        ax.scatter(r.mid, 0.02, marker="D", s=12, color=col, edgecolor="black", lw=0.5, zorder=3)
    if abs(r.svlen) >= 1000:
        near = var[(var.mid - r.mid).abs() < 0.0035].AF.max()   # redraw2: clear neighbouring lollipop heads
        ax.text(r.mid, max(af if pd.notna(af) else 0, near if pd.notna(near) else 0) + 0.035, f"{abs(r.svlen) / 1000:.1f} kb", fontsize=5.0,
                ha="center", va="bottom", rotation=90, color="#444444")
hd = [mpl.lines.Line2D([], [], color=TYPE_COL[t], marker=m_, lw=0.9, ms=3) for t, m_ in
      (("INS", "o"), ("DEL", "o"), ("MNV/COMPLEX", "D"))]
ax.legend(hd, ["insertion", "deletion", "complex (MNV)"], frameon=False, loc="upper left",
          ncol=1, handletextpad=0.3, labelspacing=0.25, bbox_to_anchor=(0.0, 1.1), fontsize=5.0)
ax.set_xlim(X0, X1)
ax.set_ylim(0, 0.62)
ax.set_yticks([0, 0.3, 0.6])
ax.set_ylabel("Allele\nfrequency", fontsize=6)
ax.set_xlabel("chr01B position (Mb)", labelpad=1)
clean(ax)
ax_d = ax
var.drop(columns=["marker"]).to_csv(DATA / "ED6d_locus_graph_SV_alleles.tsv", sep="\t", index=False)

# ---------------------------------------------------------------- e genotype boxplots
ax_e = []
rng = np.random.default_rng(7)
for k, f in enumerate((F1, F2)):
    ax = fig.add_axes([0.07 + k * 0.505, 0.075, 0.41, 0.22])   # restructure: f,g removed; e spans the row
    s = geno[geno.marker == f]
    order = ["0/0", "0/1", "1/1"]
    vals = [s[s.genotype == g].Nut_weight_g.values for g in order]
    bp = ax.boxplot(vals, positions=[0, 1, 2], widths=0.55, whis=1.5, showfliers=False, patch_artist=True,
                    medianprops=dict(color=DARK, lw=0.9), whiskerprops=dict(lw=0.6), capprops=dict(lw=0.6),
                    boxprops=dict(lw=0.6))
    for p_, c in zip(bp["boxes"], ["#CFE3F0", "#9CC5E0", "#5C9BCB"]):
        p_.set_facecolor(c)
    for i, v in enumerate(vals):
        ax.scatter(i + rng.uniform(-0.17, 0.17, len(v)), v, s=3, color=DARK, alpha=0.55, lw=0, zorder=3)
        ax.text(i, -0.35, f"n = {len(v)}", ha="center", va="top", fontsize=5.2)
    ax.set_xticks([0, 1, 2], order)
    ax.tick_params(axis="x", pad=9)
    ax.set_ylim(0, 7.6)
    ax.set_title(FOCAL_LAB[f].replace(" record", "\nrecord") if k == 0 else FOCAL_LAB[f] + "\n", fontsize=6, pad=2)
    if k == 0:
        ax.set_ylabel("Nut weight (g)")
    else:
        ax.set_ylabel("Nut weight (g)")   # restructure: panels now far apart, both keep a y axis
    ax.set_xlabel("Genotype", labelpad=1)
    clean(ax)
    ax_e.append(ax)

# restructure: panels f (photo) and g (field shell thickness) deleted (decision 2026-09-24)

for ax, s, dx in [(ax_a, "a", -0.055), (ax_b, "b", -0.055), (ax_c, "c", -0.045), (ax_d, "d", -0.045),
                  (ax_e[0], "e", -0.055)]:
    letter(fig, ax, s, dx=dx, dy=0.012 if s not in "d" else 0.004)

for ext in ("pdf", "png"):
    fig.savefig(OUT / f"Extended_Data_Fig_06.{ext}", dpi=600, facecolor="white")

summ = geno.groupby(["marker", "genotype"]).Nut_weight_g.agg(n="size", median="median", q1=lambda x: x.quantile(.25),
                                                             q3=lambda x: x.quantile(.75)).reset_index()
summ.to_csv(DATA / "ED6e_genotype_group_summary.tsv", sep="\t", index=False)
print(summ)
print("lead", L.marker, L.p, "SV thr", SV_BONF, "SNP thr", SNP_BONF)
print("lambda SV recomputed", round(lam_sv, 4), "SNP (Source Data)", lam_snp, "SNP tail m", m)
print("SV min P", p_sv[0], "n SV below thr", (p_sv < SV_BONF).sum(), "n SNP below thr", (snp.p < SNP_BONF).sum())
