#!/usr/bin/env python3
"""Supplementary Fig. 14 (old numbering; placement decided at renumbering): homozygous derived coding variation in the
commercial source groups (c, d) and SHELL coding alleles versus fruit traits (a, b).

Inputs (all in fix/headline_pop):
  g6/per_sample_load_final.tsv  per-accession derived/homozygous/heterozygous counts (g6_checks.py; 308 accessions)
  g6/g6_checks.tsv              group contrasts, chromosome jackknife s.e., accession bootstrap CI and sensitivity checks
  g3/g3_shell_genotypes_v3.tsv  sh-MPOB / sh-AVROS dosages and GWAS phenotypes (cl/g3_shell2.py)
  g3/g3_shell_known_alleles_v3.tsv  EMMAX tests (SNP kinship) of the combined sh dosage (n: all accessions with SHELL
                                    genotypes; n_lead: those also genotyped at the lead SNP, used for the lead-SNP tests)
Outputs: out/Supplementary_Fig_14.{pdf,png} (600 dpi), docx_2400/Supplementary_Fig_14.png, source_data/SF14*.tsv
"""
import sys
from pathlib import Path
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
HP = Path("${WORK_DIR}/fix/snp_repair_rerun/sf9")
sys.path.insert(0, "${WORK_DIR}/fix/beautify/common")
import palA  # noqa: E402

MM = 1 / 25.4
mpl.rcParams.update({
    "font.family": "Arial", "font.size": 7, "axes.labelsize": 7, "axes.titlesize": 7,
    "xtick.labelsize": 6, "ytick.labelsize": 6, "legend.fontsize": 6, "axes.linewidth": 0.6,
    "xtick.major.width": 0.6, "ytick.major.width": 0.6, "xtick.major.size": 2.5, "ytick.major.size": 2.5,
    "axes.unicode_minus": True, "mathtext.fontset": "custom", "mathtext.rm": "Arial",
    "mathtext.it": "Arial:italic", "mathtext.bf": "Arial:bold", "pdf.fonttype": 42, "svg.fonttype": "none",
})
DARK = "#242A30"


def clean(ax):
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False); ax.set_axisbelow(True)


d = pd.read_csv(HP / "g6/per_sample_load_final.tsv", sep="\t")
# sfplan: colour c by the repaired-data K = 4 groups (Supplementary Data 16, K4_dominant_group); 6 accessions changed
import openpyxl  # noqa: E402
_wb = openpyxl.load_workbook(HP.parents[2] / "deliver/Supplementary_Tables/Supplementary_Tables.xlsx", read_only=True)
_rows = list(_wb["Supplementary Table 16"].iter_rows(values_only=True)); _h = _rows[1]; _i = _h.index("K4_dominant_group")
_k4 = {r[0]: str(r[_i]).replace("K4_Pop", "K4P") for r in _rows[2:] if r and isinstance(r[0], str) and r[_i]}
assert set(d["sample"]) <= set(_k4)
d["k4"] = d["sample"].map(_k4)
R = pd.read_csv(HP / "g6/g6_checks.tsv", sep="\t")
P = R[R.check == "primary_datepalm"].copy()
sh = pd.read_csv(HP / "g3/g3_shell_genotypes_v3.tsv", sep="\t")
shp = pd.read_csv(HP / "g3/g3_shell_known_alleles_v3.tsv", sep="\t").set_index("trait")
k4lab = {"K4P1": "Pop1", "K4P2": "Pop2", "K4P3": "Pop3", "K4P4": "Pop4"}
archlab = {"AFR": "AFR", "HHG": "HHG", "IDB": "IDB", "SAEG": "SA-EG", "SEAA": "SEA-A", "SEAB": "SEA-B"}

fig = plt.figure(figsize=(180 * MM, 128 * MM))
# panel order follows first citation: a, b = SHELL (Results [64]); c, d = coding zygosity (Results [66])
ax_c = fig.add_axes([0.075, 0.115, 0.40, 0.35])   # coding zygosity scatter (panel c)
ax_d = fig.add_axes([0.585, 0.115, 0.40, 0.35])   # relative differences (panel d)
ax_a = fig.add_axes([0.075, 0.60, 0.40, 0.35])    # shell thickness (panel a)
ax_b = fig.add_axes([0.585, 0.60, 0.40, 0.35])    # nut weight (panel b)

# ---- c: homozygous derived missense genotypes vs coding heterozygosity
comm = d.arch.isin(["SEAA", "SAEG"])
for i, k in enumerate(["K4P1", "K4P2", "K4P3", "K4P4"]):
    s = d[d.k4 == k]
    for m, mk in ((comm, "o"), (~comm, "^")):
        t = s[m.loc[s.index]]
        ax_c.scatter(t.het_coding, t.missense_hom / 1000, s=9, marker=mk, facecolor=palA.POP[i], edgecolor="white",
                     linewidth=0.3, alpha=0.9, zorder=3)
from scipy import stats  # noqa: E402
rho, pr = stats.spearmanr(d.het_coding, d.missense_hom)
ax_c.set_xlabel("Heterozygosity at polarized coding SNPs")
ax_c.set_ylabel("Homozygous derived missense\ngenotypes per accession (×10$^{3}$)")
h = [plt.Line2D([], [], marker="o", ls="", mfc=palA.POP[i], mec="white", ms=4, label=k4lab[k]) for i, k in enumerate(k4lab)]
h += [plt.Line2D([], [], marker="o", ls="", mfc="#9A9A9A", mec="white", ms=4, label="SEA-A or SA-EG"),
      plt.Line2D([], [], marker="^", ls="", mfc="#9A9A9A", mec="white", ms=4, label="Other groups")]
ax_c.legend(handles=h, frameon=False, loc="upper right", handletextpad=0.2, borderaxespad=0.2, labelspacing=0.3)
ax_c.text(0.03, 0.05, rf"Spearman $\rho$ = {rho:.2f}, $n$ = {len(d)}", transform=ax_c.transAxes, fontsize=6, color=DARK)
clean(ax_c)

# ---- d: relative differences (commercial source groups vs AFR + IDB + SEA-B)
cats = [("missense", "Missense"), ("synonymous", "Synonymous"), ("stop_gained", "Stop-gained")]
meas = [("derived_alleles", "Derived alleles", "#A0A0A0"), ("hom_derived", "Homozygous derived genotypes", palA.EXPAND)]
x0 = np.arange(len(cats))
for j, (mcode, mlab, col) in enumerate(meas):
    y = []; e = []
    for c, _ in cats:
        r = P[(P.cat == c) & (P.measure == mcode)].iloc[0]; y.append(100 * r.rel_diff); e.append(100 * 1.96 * r.jk_se)
    xx = x0 + (j - 0.5) * 0.32
    ax_d.bar(xx, y, width=0.3, color=col, edgecolor="none", label=mlab, zorder=2)
    ax_d.errorbar(xx, y, yerr=e, fmt="none", ecolor=DARK, elinewidth=0.6, capsize=1.5, zorder=3)
    for xi, yi, ei in zip(xx, y, e):
        ax_d.text(xi, yi + ei + 1.2 if yi >= 0 else yi - ei - 3, f"{yi:+.1f}", ha="center", va="bottom", fontsize=5.5, color=DARK)
ax_d.axhline(0, color=DARK, lw=0.6)
ax_d.set_xticks(x0); ax_d.set_xticklabels([f"{l}\n({int(P[(P.cat == c)].n_sites.iloc[0]):,} sites)" for c, l in cats])
ax_d.set_ylabel("Difference, SEA-A + SA-EG vs\nAFR + IDB + SEA-B (%)")
ax_d.set_ylim(-5, 55)
ax_d.legend(frameon=False, loc="upper left", handlelength=1.0, borderaxespad=0.2)
clean(ax_d)


# ---- a, b: SHELL coding alleles
CLS = [(0, "0", palA.A["ancestry", "dura_like"]), (1, "1", palA.A["ancestry", "mixed"]), (2, "2", palA.A["ancestry", "pisifera_like"])]
rng = np.random.default_rng(3)


def shell_panel(ax, trait, ylab, key):
    s = sh[sh[trait].notna() & sh.sh_total.notna()]
    data, pos, cols, labs = [], [], [], []
    for v, lab, col in CLS:
        y = s.loc[s.sh_total == v, trait].to_numpy()
        if len(y) == 0:
            continue
        data.append(y); pos.append(v); cols.append(col); labs.append(f"{lab}\n$n$ = {len(y)}")
    bp = ax.boxplot([x for x in data if len(x) > 1], positions=[p for p, x in zip(pos, data) if len(x) > 1], widths=0.5,
                    patch_artist=True, showfliers=False, medianprops=dict(color=DARK, lw=0.8),
                    whiskerprops=dict(color=DARK, lw=0.6), capprops=dict(color=DARK, lw=0.6), boxprops=dict(lw=0.6, color=DARK))
    for patch, col in zip(bp["boxes"], [c for c, x in zip(cols, data) if len(x) > 1]):
        patch.set_facecolor(palA.tint(col, 0.45))
    for p, y, col in zip(pos, data, cols):
        ax.scatter(p + rng.uniform(-0.17, 0.17, len(y)), y, s=5, color=col, edgecolor="none", alpha=0.85, zorder=3)
    ax.set_xticks(pos); ax.set_xticklabels(labs)
    ax.set_xlim(-0.6, 2.6)
    ax.set_xlabel("Mutant $SHELL$ alleles ($sh^{\\mathrm{MPOB}}$ + $sh^{\\mathrm{AVROS}}$)")
    ax.set_ylabel(ylab)
    p = shp.loc[trait, "P_shTotal"]
    ax.text(0.97, 0.95, rf"EMMAX $P$ = ${palA.sci_tex(p)}$", transform=ax.transAxes, ha="right", va="top", fontsize=6, color=DARK)
    clean(ax)
    return s


s_a = shell_panel(ax_a, "Shell_thickness_mm", "Shell thickness (mm)", "a")
s_b = shell_panel(ax_b, "Nut_weight_g", "Nut weight (g)", "b")

for ax, s in [(ax_a, "a"), (ax_b, "b"), (ax_c, "c"), (ax_d, "d")]:
    bb = ax.get_position()
    fig.text(bb.x0 - 0.062, bb.y1 + 0.012, s, fontsize=palA.LETTER_PT, fontweight="bold", va="bottom")

(HERE / "out").mkdir(exist_ok=True)
for ext in ("pdf", "png"):
    fig.savefig(HERE / "out" / f"Supplementary_Fig_14.{ext}", dpi=600, facecolor="white")
fig.savefig(HERE / "docx_2400" / "Supplementary_Fig_14.png", dpi=2400 / (180 * MM), facecolor="white")

# ---- Source Data
SD = HERE / "source_data"; SD.mkdir(exist_ok=True)
a = d.copy()
a["Archive_group"] = a.arch.map(archlab); a["K4_group"] = a.k4.map(k4lab)
a["Contrast_group"] = np.where(a.arch.isin(["SEAA", "SAEG"]), "Commercial source (SEA-A + SA-EG)",
                               np.where(a.arch.isin(["AFR", "IDB", "SEAB"]), "AFR + IDB + SEA-B", "HHG (not in contrast)"))
cols = ["sample", "Archive_group", "K4_group", "Contrast_group", "het_coding"] + \
       [f"{c}_{m}" for c in ("missense", "synonymous", "stop_gained") for m in ("n", "der", "hom", "het")] + ["altfrac"]
a = a[cols].rename(columns={"sample": "Sample_ID", "het_coding": "Heterozygosity_polarized_coding",
                            "altfrac": "ALT_allele_fraction_synonymous"})
a.rename(columns=lambda c: c.replace("_n", "_called_sites").replace("_der", "_derived_alleles").replace("_hom", "_hom_derived")
         .replace("_het", "_het") if c.split("_")[-1] in ("n", "der", "hom", "het") else c).to_csv(SD / "SF14c_per_accession.tsv", sep="\t", index=False)
b = R.copy()
b.columns = ["Analysis", "Class", "Polarized_sites", "Measure", "Median_SEA-A+SA-EG", "Median_AFR+IDB+SEA-B",
             "Relative_difference", "Chromosome_jackknife_se", "Accession_bootstrap_2.5pct", "Accession_bootstrap_97.5pct", "P_Mann_Whitney"]
b.to_csv(SD / "SF14d_group_contrasts.tsv", sep="\t", index=False)
for tr, nm in (("Shell_thickness_mm", "SF14a_shell_thickness"), ("Nut_weight_g", "SF14b_nut_weight")):
    s = sh[sh[tr].notna()][["id", "MPOB", "AVROS", "sh_total", "lead_A1", tr]].copy()
    s.columns = ["Sample_ID", "sh_MPOB_dosage (chr01B:3,259,207 G)", "sh_AVROS_dosage (chr01B:3,259,200 A)",
                 "Mutant_SHELL_alleles", "Lead_SNP_chr01B:3,153,030_A1_dosage", tr]
    s.to_csv(SD / f"{nm}.tsv", sep="\t", index=False)
t = shp[["n", "n_lead", "P_MPOB", "P_AVROS", "P_shTotal", "P_lead", "r2_sh_lead", "P_lead_given_sh", "P_sh_given_lead",
         "n_sh0", "n_sh1", "n_sh2", "med_sh0", "med_sh1", "med_sh2"]].reset_index()
t.to_csv(SD / "SF14ab_EMMAX_tests.tsv", sep="\t", index=False)
print("rho", rho, pr)
print(P[["cat", "measure", "rel_diff", "jk_se"]])
