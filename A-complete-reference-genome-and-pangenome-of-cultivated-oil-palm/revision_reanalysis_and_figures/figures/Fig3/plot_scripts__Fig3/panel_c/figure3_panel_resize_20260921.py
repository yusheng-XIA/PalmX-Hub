#!/usr/bin/env python3
"""Re-render Figure 3 panels at the author-requested physical dimensions.

This is a layout-only export. It reads reviewed source tables / accepted panel
inputs and does not recalculate classifications or alter the scientific values.
"""
from __future__ import annotations

import hashlib
import math
import os
import platform
import shutil
import sys
import tempfile
from pathlib import Path

import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
from matplotlib.lines import Line2D
from matplotlib.patches import Patch, PathPatch, Rectangle
from matplotlib.path import Path as MplPath
import numpy as np
import pandas as pd
from scipy import stats

OUT = Path(__file__).resolve().parent
MM = 1 / 25.4

ANALYSIS = Path("${ANALYSIS_DIR}")
REV = ANALYSIS / "22_answer_reviews/00_ms/05_MS/MS_revision_3/Final_20260823_VectorRevision/8月29日/Revised_Panels_20260831"
EF = ANALYSIS / "22_answer_reviews/00_ms/05_MS/new_revision/Figure3ef_four_classes_20260904"
REF = ANALYSIS / "22_answer_reviews/00_ms/03_V3/03_figure3/04_FL_TN_ASE_reference_panels"
BREED = ANALYSIS / "22_answer_reviews/00_ms/03_V3/03_figure3/13_comprehensive_primitive_ancestry_breeding_20260805/runs/RUN-COMP-BREEDING-DPM-20260805-001"

SRC_C = REV / "Figure3c_plotdata.tsv"
SRC_D = REV / "Figure3d_intersection_plotdata.tsv"
SRC_E = EF / "Figure3e_four_ASE_classes_plotdata.tsv"
SRC_E_TESTS = EF / "Figure3e_four_ASE_classes_tests.tsv"
SRC_F = EF / "Figure3f_four_ASE_classes_plotdata.tsv"
SRC_G = REV / "Figure3g_plotdata.tsv"
SRC_H = REF / "source_04B_TN_regulatory_classes_19stages.tsv"
SRC_I = REF / "source_04C_TN_cis_contribution.tsv"
SRC_J = REV / "Figure3j_plotdata_gene_clustered.tsv"
SRC_J_TESTS = REV / "Figure3j_gene_trajectory_permutation_tests.tsv"
SRC_K = REF / "source_06A_FL_TN_trait_complement.tsv"
SRC_L_BINS = ANALYSIS / "22_answer_reviews/00_ms/03_V3/03_figure3/11_diploid_32chrom_ancestry_breeding_20260805/results/diploid_ancestry_bins.tsv"
SRC_L_BASE = BREED / "results/P9_final_attempt003_favorable_rescue/favorable_DPM_targets.tsv"
SRC_L_RESCUE = BREED / "results/P9_final_attempt003_favorable_rescue/favorable_haplotype_rescue_candidates.tsv"

INPUTS = [SRC_C, SRC_D, SRC_E, SRC_E_TESTS, SRC_F, SRC_G, SRC_H, SRC_I,
          SRC_J, SRC_J_TESTS, SRC_K, SRC_L_BINS, SRC_L_BASE, SRC_L_RESCUE]

SIZES_MM = {
    "Fig3ab": (58, 52), "Fig3c": (123, 52), "Fig3d": (38, 42),
    "Fig3e": (72, 42), "Fig3f": (70, 42), "Fig3g": (99, 47),
    "Fig3h": (86, 47), "Fig3i": (58, 40), "Fig3j": (56, 40),
    "Fig3k": (69, 40), "Fig3l": (183, 54),
}

CLASSES = ["NoDiff", "HapDom", "Sub", "NoASE"]
CLASS_COLORS = {"NoDiff": "#F8766D", "HapDom": "#F1DC91", "Sub": "#8DD3C7", "NoASE": "#74B3D4"}
STAGES = ["0d", "15d", "35d", "50d", "65d", "80d", "95d", "110d", "125d", "140d", "155d", "170d", "185d", "12h", "24h", "36h", "48h", "60h", "72h"]
PHASES = ["Days 0–65", "Days 80–140", "Days 155–185", "Hours 12–72"]
REG_ORDER = ["I.Cis_only", "II.Trans_only", "III.Cis_trans_enhancing", "IV.Cis_trans_compensating", "V.Compensatory", "VI.Conserved", "VII.Ambiguous"]
CMAP = LinearSegmentedColormap.from_list("dpm", ["#3767A6", "#8FAED0", "#F7F5EE", "#D8B64C", "#E56565"])


def style() -> None:
    mpl.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Liberation Sans", "DejaVu Sans"],
        "font.size": 6.0, "axes.labelsize": 6.0, "axes.titlesize": 7.0,
        "xtick.labelsize": 6.0, "ytick.labelsize": 6.0,
        "legend.fontsize": 6.0, "axes.linewidth": 0.65,
        "xtick.major.width": 0.65, "ytick.major.width": 0.65,
        "xtick.major.size": 2.4, "ytick.major.size": 2.4,
        "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
        "savefig.facecolor": "white", "figure.facecolor": "white",
    })


def canvas(stem: str):
    w, h = SIZES_MM[stem]
    return plt.figure(figsize=(w * MM, h * MM), facecolor="white")


def letter(fig, text: str, x=0.012, y=0.985):
    fig.text(x, y, text, ha="left", va="top", fontsize=8, fontweight="bold")


def clean(ax):
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(direction="out")


def save(fig, stage: Path, stem: str):
    meta = {
        "Title": f"{stem} — {SIZES_MM[stem][0]} x {SIZES_MM[stem][1]} mm",
        "Subject": "Layout-only re-render from reviewed Figure 3 source data",
        "Creator": Path(__file__).name,
    }
    fig.savefig(stage / f"{stem}.pdf", metadata=meta, facecolor="white")
    fig.savefig(stage / f"{stem}_600dpi.png", dpi=600, facecolor="white")
    plt.close(fig)


def panel_ab(stage: Path):
    fig = canvas("Fig3ab")
    # a: faithful compact redraw of the accepted Hap1/Hap2 allele schematic.
    ax = fig.add_axes([0.08, 0.57, 0.88, 0.38]); ax.set_axis_off()
    y1, y2 = 0.76, 0.25
    ax.plot([0.10, 0.98], [y1, y1], color="#6C8AA4", lw=0.6)
    ax.plot([0.10, 0.98], [y2, y2], color="#B17A5B", lw=0.6)
    cols = {"Biallelic": "#79B5D4", "Allele with same CDS": "#D48472", "Haplotype-specific": "#F2C77E", "Gene absent": "white"}
    xstarts = [0.12, 0.28, 0.43, 0.59, 0.75, 0.88]
    top_types = ["Allele with same CDS", "Biallelic", "Biallelic", "Haplotype-specific", "Biallelic", "Gene absent"]
    bot_types = ["Allele with same CDS", "Biallelic", "Biallelic", "Gene absent", "Biallelic", "Haplotype-specific"]
    for y, types, edge in [(y1, top_types, "#6C8AA4"), (y2, bot_types, "#B17A5B")]:
        for x, typ in zip(xstarts, types):
            ax.add_patch(Rectangle((x, y - 0.07), 0.085, 0.14, fc=cols[typ], ec=edge if typ != "Gene absent" else "#999999", lw=0.45, ls="--" if typ == "Gene absent" else "-"))
    ax.text(0.00, y1, "Hap1", va="center", fontweight="bold")
    ax.text(0.00, y2, "Hap2", va="center", fontweight="bold")
    handles = [Patch(fc=cols[k], ec="#777777", lw=0.4, ls="--" if k == "Gene absent" else "-", label=k) for k in cols]
    fig.legend(handles=handles, frameon=False, ncol=2, loc="upper center", bbox_to_anchor=(0.60, 0.57), handlelength=1.0, columnspacing=0.8, handletextpad=0.35)
    letter(fig, "a")

    # b: exact accepted counts from plot_allele_classification.R, relabelled as in the example.
    bx = fig.add_axes([0.17, 0.12, 0.78, 0.31])
    labels = ["FL_hap1", "FL_hap2", "TN_hap1", "TN_hap2"]
    counts = np.array([[27106, 264, 3214], [27106, 264, 7430], [18439, 7747, 7720], [18439, 7747, 8486]], float)
    pct = counts / counts.sum(axis=1, keepdims=True) * 100
    bar_colors = ["#79B5D4", "#D48472", "#F1DC91"]
    bottom = np.zeros(4)
    for j in range(3):
        bx.bar(np.arange(4), pct[:, j], bottom=bottom, width=0.78, color=bar_colors[j], edgecolor="white", lw=0.35)
        bottom += pct[:, j]
    bx.set_ylim(0, 100); bx.set_yticks([0, 25, 50, 75, 100]); bx.set_ylabel("Proportion of genes (%)")
    bx.set_xticks(range(4), labels, rotation=0); bx.spines[["top", "right"]].set_visible(False); bx.tick_params(axis="x", length=0)
    fig.text(0.012, 0.435, "b", ha="left", va="top", fontsize=8, fontweight="bold")
    save(fig, stage, "Fig3ab")


def panel_c(stage: Path):
    d = pd.read_csv(SRC_C, sep="\t")
    fig = canvas("Fig3c"); ax = fig.add_axes([0.085, 0.25, 0.90, 0.56])
    x = np.arange(19); width = 0.30
    spec = [("TN", "Allele_A_biased", -width/2, 1, "#79B9D8", "TN A (Dura/TK-like)"),
            ("FL", "Allele_A_biased", width/2, 1, "#AEC9E0", "FL A (Africa hap2)"),
            ("TN", "Allele_B_biased", -width/2, -1, "#E99AA7", "TN B (Pisifera/NS-like)"),
            ("FL", "Allele_B_biased", width/2, -1, "#F1D994", "FL B (American hap1)")]
    for analysis, call, off, sign, col, lab in spec:
        vals = d[(d.analysis == analysis) & (d.ase_call == call)].set_index("stage").reindex(STAGES).genes.to_numpy()
        ax.bar(x+off, sign*vals, width=width, color=col, edgecolor="white", lw=.3, label=lab)
    lim = 5500
    ax.axhline(0, color="#333333", lw=.65); ax.axvline(12.5, color="#777777", lw=.6, ls="--")
    ax.axvspan(-.5, 12.5, color="#F2F7F5", zorder=-2); ax.axvspan(12.5, 18.5, color="#FFF7EA", zorder=-2)
    ax.set(xlim=(-.6,18.6), ylim=(-lim,lim), ylabel="Number of ASE genes")
    ax.set_xticks(x, STAGES, rotation=55, ha="right")
    ax.set_yticks([-4000,-2000,0,2000,4000], ["4,000","2,000","0","2,000","4,000"])
    ax.text(6, 5000, "Development (0–185 d)", ha="center", va="top")
    ax.text(15.5, 5000, "Postharvest (12–72 h)", ha="center", va="top")
    ax.text(18.5, 4000, "Allele A biased", ha="right", color="#555555")
    ax.text(18.5, -4250, "Allele B biased", ha="right", color="#555555")
    clean(ax)
    ax.legend(frameon=False, ncol=2, loc="lower left", bbox_to_anchor=(0.23, 1.01), columnspacing=.8, handlelength=1.1)
    letter(fig, "c")
    save(fig, stage, "Fig3c")


def flow_patch(ax, x0, x1, y0a, y0b, y1a, y1b, color):
    c = (x1-x0)*0.48
    verts = [(x0,y0a),(x0+c,y0a),(x1-c,y1a),(x1,y1a),(x1,y1b),(x1-c,y1b),(x0+c,y0b),(x0,y0b),(x0,y0a)]
    codes = [MplPath.MOVETO,MplPath.CURVE4,MplPath.CURVE4,MplPath.CURVE4,MplPath.LINETO,MplPath.CURVE4,MplPath.CURVE4,MplPath.CURVE4,MplPath.CLOSEPOLY]
    ax.add_patch(PathPatch(MplPath(verts,codes), fc=color, ec="none", alpha=.35, zorder=1))


def panel_d(stage: Path):
    d = pd.read_csv(SRC_D, sep="\t")
    fig = canvas("Fig3d"); ax = fig.add_axes([0.08, 0.12, 0.84, 0.69]); ax.set_axis_off()
    x0,x1=.10,.90; bw=.08
    for r in d.itertuples():
        flow_patch(ax, x0+bw/2, x1-bw/2, r.left_bottom, r.left_top, r.right_bottom, r.right_top, CLASS_COLORS[r.overall_class_FL])
    left = d.drop_duplicates("overall_class_FL").set_index("overall_class_FL")
    right = d.drop_duplicates("overall_class_TN").set_index("overall_class_TN")
    for cls in CLASSES:
        l = left.loc[cls]; r = right.loc[cls]
        ax.add_patch(Rectangle((x0-bw/2,l.left_bottom),bw,l.left_top-l.left_bottom,fc=CLASS_COLORS[cls],ec="#333",lw=.35,zorder=3))
        ax.add_patch(Rectangle((x1-bw/2,r.right_bottom),bw,r.right_top-r.right_bottom,fc=CLASS_COLORS[cls],ec="#333",lw=.35,zorder=3))
    ax.text(x0,-.035,"FL",ha="center",va="top"); ax.text(x1,-.035,"TN",ha="center",va="top")
    fig.text(.5,.035,"Shared ASE-eligible\northogroups, n = 7,059",ha="center",va="bottom")
    handles=[Patch(fc=CLASS_COLORS[c],ec="none",label=c) for c in CLASSES]
    fig.legend(handles=handles,frameon=False,ncol=2,loc="upper center",bbox_to_anchor=(.58,.96),columnspacing=.6,handlelength=.9,handletextpad=.3)
    letter(fig,"d")
    save(fig,stage,"Fig3d")


def p_text(p):
    return f"P = {p:.3f}" if p >= .001 else f"P = {p:.1e}"


def panel_e(stage: Path):
    d=pd.read_csv(SRC_E,sep="\t"); tests=pd.read_csv(SRC_E_TESTS,sep="\t")
    fig=canvas("Fig3e"); ax=fig.add_axes([.16,.20,.82,.60])
    regions=[("upstream2kb_per_kb","Upstream\n2 kb"),("gene_body_per_kb","Gene"),("downstream2kb_per_kb","Downstream\n2 kb")]
    pos=[]; arrays=[]
    for ri,(reg,_) in enumerate(regions):
        for ci,cls in enumerate(CLASSES):
            pos.append(ri*5+ci); arrays.append(d[(d.region==reg)&(d.overall_class==cls)].SNPs_per_1000bp.to_numpy())
    bp=ax.boxplot(arrays,positions=pos,widths=.7,showfliers=False,patch_artist=True,medianprops={"lw":.65},whiskerprops={"lw":.5},capprops={"lw":.5})
    for i,box in enumerate(bp["boxes"]):
        col=CLASS_COLORS[CLASSES[i%4]]; box.set(fc=mpl.colors.to_rgba(col,.22),ec=col,lw=.7); bp["medians"][i].set_color(col)
        for a in bp["whiskers"][2*i:2*i+2]+bp["caps"][2*i:2*i+2]: a.set_color(col)
    kw=tests[tests.test=="Kruskal-Wallis"].set_index("region")
    for ri,(reg,_) in enumerate(regions):
        letters=dict(x.split(":") for x in kw.loc[reg,"letters"].split(";"))
        ax.text(ri*5+1.5,48,p_text(float(kw.loc[reg,"pvalue"])),ha="center")
        for ci,cls in enumerate(CLASSES): ax.text(ri*5+ci,43.5,letters[cls],ha="center",fontweight="bold")
    ax.set(xlim=(-.8,13.8),ylim=(0,52),ylabel="SNPs per 1,000 bp")
    ax.set_xticks([1.5,6.5,11.5],[x[1] for x in regions]); clean(ax)
    ax.text(-.65,50,"TN",fontweight="bold",ha="left",va="bottom")
    fig.legend(handles=[Patch(fc=mpl.colors.to_rgba(CLASS_COLORS[c],.22),ec=CLASS_COLORS[c],label=c) for c in CLASSES],frameon=False,ncol=4,loc="upper center",bbox_to_anchor=(.62,.985),columnspacing=.6,handlelength=.8,handletextpad=.25)
    letter(fig,"e")
    save(fig,stage,"Fig3e")


def kde(vals,xmax):
    a=pd.to_numeric(vals,errors="coerce").to_numpy(float); a=a[np.isfinite(a)&(a>=0)&(a<=xmax)]
    x=np.linspace(0,xmax,350); return x,stats.gaussian_kde(a)(x)


def panel_f(stage: Path):
    d=pd.read_csv(SRC_F,sep="\t")
    fig=canvas("Fig3f"); ax=fig.add_axes([.15,.19,.82,.73]); ins=ax.inset_axes([.54,.54,.43,.42])
    styles={"FL":"-","TN":(0,(6,3))}
    for mat in ["FL","TN"]:
        for cls in CLASSES:
            q=d[(d.analysis==mat)&(d.ASE_type==cls)]
            x,y=kde(q.Ka_Ks,3); ax.plot(x,y,color=CLASS_COLORS[cls],lw=1.0,ls=styles[mat])
            x,y=kde(q.Ks,.1); ins.plot(x,y,color=CLASS_COLORS[cls],lw=.75,ls=styles[mat])
    ax.axvline(1,color="#777",lw=.6,ls="--"); ax.set(xlim=(0,3),ylim=(0,None),xlabel="Ka/Ks ratio",ylabel="Density"); ax.set_xticks([0,1,2,3]); clean(ax)
    ins.set(xlim=(0,.1),ylim=(0,None),xlabel="Ks",ylabel="Density"); ins.set_xticks([0,.05,.10]); ins.tick_params(labelsize=6)
    for s in ins.spines.values(): s.set_linewidth(.5); s.set_color("#888")
    lc=[Line2D([0],[0],color=CLASS_COLORS[c],lw=1.2,label=c) for c in CLASSES]
    lm=[Line2D([0],[0],color="#333",lw=1.1,ls=styles[m],label=m) for m in ["FL","TN"]]
    lg=ax.legend(handles=lc,title="ASE class",frameon=False,ncol=2,loc="lower right",bbox_to_anchor=(1.00,.01),columnspacing=.55,handlelength=1.2,handletextpad=.25)
    ax.add_artist(lg); ax.legend(handles=lm,title="Material",frameon=False,ncol=2,loc="lower right",bbox_to_anchor=(1.0,.34),columnspacing=.55,handlelength=1.6,handletextpad=.25)
    letter(fig,"f")
    save(fig,stage,"Fig3f")


def panel_g(stage: Path):
    d=pd.read_csv(SRC_G,sep="\t"); modes=[f"M{i}" for i in range(1,13)]+["Conserved"]
    val=d.pivot(index="stage",columns="mode",values="percentage").reindex(index=STAGES,columns=modes).fillna(0)
    fig=canvas("Fig3g"); ax=fig.add_axes([.10,.12,.83,.78])
    im=ax.imshow(val.to_numpy(),cmap=CMAP,aspect="auto",vmin=3,vmax=16)
    xlabels=modes[:-1]+["Cons."]
    ax.set_xticks(range(len(modes)),xlabels,fontweight="bold"); ax.xaxis.tick_top(); ax.tick_params(top=True,bottom=False,labeltop=True,labelbottom=False)
    ax.set_yticks(range(19),STAGES); ax.axhline(12.5,color="white",lw=1.2)
    cb=fig.colorbar(im,ax=ax,pad=.012,fraction=.03); cb.set_label("Within-stage genes (%)"); cb.ax.tick_params(labelsize=6)
    for s in ax.spines.values(): s.set_linewidth(.5); s.set_color("#777")
    letter(fig,"g")
    save(fig,stage,"Fig3g")


def panel_h(stage: Path):
    d=pd.read_csv(SRC_H,sep="\t"); piv=d.pivot(index="regulatory_class",columns="stage",values="percentage").reindex(index=REG_ORDER,columns=STAGES).fillna(0)
    labels=["I Cis only","II Trans only","III Cis + trans\nenhancing","IV Cis + trans\ncompensating","V Compensatory","VI Conserved","VII Ambiguous"]
    fig=canvas("Fig3h"); ax=fig.add_axes([.25,.22,.65,.62])
    xs,ys=np.meshgrid(np.arange(19),np.arange(7)); vals=piv.to_numpy()
    sc=ax.scatter(xs.ravel(),ys.ravel(),s=6+65*vals.ravel()/55,c=vals.ravel(),cmap="Reds",vmin=0,vmax=55,edgecolor="white",lw=.25)
    ax.set_xticks(range(19),STAGES,rotation=70,ha="right"); ax.set_yticks(range(7),labels); ax.invert_yaxis(); ax.set_xlim(-.6,18.6); ax.axvline(12.5,color="#777",lw=.6,ls="--")
    ax.text(6,-1.05,"Seed development",ha="center",color="#3D7A69"); ax.text(15.5,-1.05,"Postharvest",ha="center",color="#B68A25")
    clean(ax); ax.set_xlabel("Stage")
    cb=fig.colorbar(sc,ax=ax,pad=.012,fraction=.035); cb.set_label("Genes within stage (%)"); cb.ax.tick_params(labelsize=6)
    letter(fig,"h")
    save(fig,stage,"Fig3h")


def panel_i(stage: Path):
    d=pd.read_csv(SRC_I,sep="\t"); bins=["0–1","1–2","2–3","3–4","4+"]; colors=["#E56565","#D8B64C","#159D82","#3767A6"]; markers=["o","s","^","D"]
    labs=["Early (0–65 d)","Mid (80–140 d)","Late (155–185 d)","Postharvest (12–72 h)"]
    fig=canvas("Fig3i"); ax=fig.add_axes([.19,.25,.78,.53]); x=np.arange(5)
    for phase,lab,col,mk,off in zip(PHASES,labs,colors,markers,[-.075,-.025,.025,.075]):
        q=d[d.stage_group==phase].set_index("abs_A_bin").reindex(bins); y=q["median"].to_numpy(); err=np.vstack([y-q.ci_low.to_numpy(),q.ci_high.to_numpy()-y])
        ax.errorbar(x+off,y,yerr=err,color=col,marker=mk,mfc="white",ms=2.5,lw=.7,capsize=1.2,label=lab)
    ax.axhline(.5,color="#666",lw=.6,ls="--"); ax.text(4.25,.502,"Cis = trans",ha="right",va="bottom")
    ax.set(xlim=(-.3,4.3),ylim=(.28,.515),xlabel="|A| (log2 scale)",ylabel="Cis contribution |B|/(|B|+|A-B|)")
    ax.set_xticks(x,["0–1","1–2","2–3","3–4","≥4"]); clean(ax)
    ax.legend(frameon=False,ncol=2,loc="lower left",bbox_to_anchor=(-.03,1.02),columnspacing=.5,handlelength=1.1,handletextpad=.25)
    letter(fig,"i")
    save(fig,stage,"Fig3i")


def panel_j(stage: Path):
    d=pd.read_csv(SRC_J,sep="\t"); tests=pd.read_csv(SRC_J_TESTS,sep="\t"); rows=["PDO","DO","ODO"]
    short=["I","II","III","IV","V","VI","VII"]
    pct=d.pivot(index="inheritance_class",columns="regulatory_class",values="row_percentage").reindex(index=rows,columns=REG_ORDER).to_numpy()
    res=d.pivot(index="inheritance_class",columns="regulatory_class",values="trajectory_permutation_residual").reindex(index=rows,columns=REG_ORDER).to_numpy()
    vmax=75; norm=TwoSlopeNorm(vmin=-vmax,vcenter=0,vmax=vmax)
    fig=canvas("Fig3j"); ax=fig.add_axes([.14,.25,.82,.50]); im=ax.imshow(res,cmap=CMAP,norm=norm,aspect="auto")
    for r in range(3):
        for c in range(7): ax.text(c,r,f"{pct[r,c]:.1f}%",ha="center",va="center",fontsize=6,color="white" if abs(res[r,c])>42 else "#222")
    ax.set_xticks(range(7),short); ax.set_yticks(range(3),rows); ax.spines[:].set_visible(False)
    cax=fig.add_axes([.14,.82,.82,.025]); cb=fig.colorbar(im,cax=cax,orientation="horizontal"); cb.set_ticks([-75,-50,-25,0,25,50,75]); cb.ax.xaxis.set_ticks_position("top"); cb.ax.tick_params(labelsize=6,pad=1)
    fig.text(.96,.985,"Pearson residual",ha="right",va="top")
    gp=float(tests.loc[tests.scope=="global_3x7","empirical_p"].iloc[0]); fig.text(.55,.765,f"Global permutation P = {gp:.4f}",ha="center")
    letter(fig,"j")
    save(fig,stage,"Fig3j")


def panel_k(stage: Path):
    d=pd.read_csv(SRC_K,sep="\t")
    modules=["Oil biosynthesis & storage","TAG assembly & oil body","De-novo / saturated FA","Unsaturated FA","Lipid oxidation / antioxidant","Shell / cell wall / lignin"]
    abbrev=["OBS","TOF","DSF","UFA","LOD","SCL"]
    fig=canvas("Fig3k"); ax=fig.add_axes([.12,.24,.84,.59])
    vmax=max(1,float(np.nanpercentile(np.abs(d.robust_median_log2_ratio),98))); norm=TwoSlopeNorm(vmin=-vmax,vcenter=0,vmax=vmax)
    for yi,mod in enumerate(modules):
        for pi,phase in enumerate(PHASES):
            for ai,analysis in enumerate(["FL","TN"]):
                q=d[(d.trait_module==mod)&(d.stage_group==phase)&(d.analysis==analysis)]
                if q.empty: continue
                z=q.iloc[0]; x=pi*2+ai
                ax.scatter(x,yi,s=8+.22*z.robust_ASE_percentage,c=[z.robust_median_log2_ratio],cmap="coolwarm",norm=norm,marker="o" if analysis=="FL" else "s",edgecolor="#E56565" if analysis=="FL" else "#3767A6",lw=.5)
    ax.set_xlim(-.6,7.6); ax.set_ylim(5.6,-.7); ax.set_yticks(range(6),abbrev); ax.set_xticks(range(8),["FL","TN"]*4); ax.xaxis.tick_top(); ax.tick_params(top=True,bottom=False,labeltop=True,labelbottom=False,length=0)
    for xpos,name in zip([.22,.44,.66,.87],["Early","Middle","Late","Postharvest"]): fig.text(xpos,.96,name,ha="center",fontweight="bold")
    for pi in range(4): ax.axvspan(pi*2-.5,pi*2+1.5,color="#F7F8F8" if pi<3 else "#FFF7E8",zorder=-2)
    for s in ax.spines.values(): s.set_visible(False)
    sm=mpl.cm.ScalarMappable(norm=norm,cmap="coolwarm"); cb=fig.colorbar(sm,ax=ax,orientation="horizontal",pad=.16,fraction=.06); cb.set_label("Median log2(A/B)"); cb.ax.tick_params(labelsize=6)
    sizes=[50,70,90]; handles=[plt.scatter([],[],s=8+.22*x,fc="#DDE2E5",ec="#777",lw=.4,label=f"{x}%") for x in sizes]
    ax.legend(handles=handles,title="Robust ASE",frameon=False,ncol=3,loc="lower center",bbox_to_anchor=(.52,-.42),handletextpad=.2,columnspacing=.5)
    letter(fig,"k")
    save(fig,stage,"Fig3k")


def panel_l(stage: Path):
    bins=pd.read_csv(SRC_L_BINS,sep="\t"); bins=bins[bins.Individual=="FL"].copy()
    base=pd.read_csv(SRC_L_BASE,sep="\t",low_memory=False); rescue=pd.read_csv(SRC_L_RESCUE,sep="\t",low_memory=False)
    fig=canvas("Fig3l")
    axes=[fig.add_axes([.055,.35,.43,.55]),fig.add_axes([.54,.35,.43,.55])]
    anc_colors={"Dura":"#3767A6","Pisifera":"#E56565","Meizhou4":"#D8B64C"}; geno={"D":"#3767A6","P":"#E56565","M":"#D8B64C"}
    action={"INTRODUCE_OR_TUNE":("^","#159D82"),"RETAIN_FL":("^","#4F70B5"),"TIMING_SCREEN":(">","#E3A018")}
    actual={"TN_h1":"#00897B","TN_h2":"#6BC5BA","FL_HapA":"#7651A8","FL_HapB":"#B091CF"}
    for col,ax in enumerate(axes):
        chroms=[f"chr{i:02d}" for i in range(1+8*col,9+8*col)]
        for yi,chrom in enumerate(chroms):
            q=bins[bins.Chromosome==chrom].sort_values("Bin_index"); y=7-yi; lengths=[]
            for hap,dy in [(1,.13),(2,-.13)]:
                endcol=f"Hap{hap}_end0"; startcol=f"Hap{hap}_start0"; anccol=f"Hap{hap}_ancestry"; length=float(q[endcol].max())/1e6; lengths.append(length)
                for r in q.itertuples():
                    x0=getattr(r,startcol)/1e6; x1=getattr(r,endcol)/1e6; a=str(getattr(r,anccol)); known=a in anc_colors
                    ax.add_patch(Rectangle((x0,y+dy-.045),x1-x0,.09,fc=anc_colors.get(a,"white"),ec="none" if known else "#888",lw=.15,hatch=None if known else "////"))
                ax.add_patch(Rectangle((0,y+dy-.045),length,.09,fill=False,ec="#444",lw=.25))
            t=base[base.Chromosome==chrom]
            for r in t.itertuples():
                frac=float(r.chromosome_fraction); alleles=str(r.target_DPM_genotype).split("/")
                for hap,dy,L in [(0,.23,lengths[0]),(1,-.23,lengths[1])]:
                    if hap < len(alleles) and alleles[hap] in geno: ax.scatter(frac*L,y+dy,s=3.0,marker="s",c=geno[alleles[hap]],edgecolor="#222",lw=.15,zorder=5)
                mk,colr=action.get(str(r.breeding_action),("^","#999")); ax.scatter(frac*max(lengths),y+.32,s=3.5,marker=mk,c=colr,edgecolor="none",zorder=6)
            rr=rescue[rescue.Chromosome==chrom]
            for r in rr.itertuples():
                frac=float(r.chromosome_fraction); ax.scatter(frac*max(lengths),y+.32,s=3.5,marker="D",c="#8A5FB2",edgecolor="none",zorder=6)
                for hapname in [str(getattr(r,"TN_ASE_hap","")),str(getattr(r,"FL_ASE_hap",""))]:
                    if hapname in actual:
                        hi=0 if hapname in {"TN_h1","FL_HapA"} else 1; ax.scatter(frac*lengths[hi],y+(.23 if hi==0 else -.23),s=3,marker="s",c=actual[hapname],edgecolor="#222",lw=.15,zorder=6)
            ax.text(-4.5,y,chrom.replace("chr0","chr").replace("chr","chr"),ha="right",va="center")
        ax.set(xlim=(-8,202),ylim=(-.55,7.65)); ax.set_yticks([]); ax.xaxis.tick_top(); ax.set_xticks(np.arange(0,201,25)); ax.tick_params(axis="x",length=2,pad=1)
        for s in ["left","right","bottom"]: ax.spines[s].set_visible(False)
        ax.spines["top"].set_color("#777")
    fig.text(.055,.985,"Upper track: H1; lower track: H2",ha="left",va="top")
    letter(fig,"l")
    # Four compact legends; labels retain the accepted source terminology.
    leg1=[Patch(fc=anc_colors[k],ec="#444",lw=.25,label=k) for k in ["Dura","Pisifera","Meizhou4"]]+[Patch(fc="white",ec="#888",hatch="////",label="Pending")]
    leg2=[Line2D([0],[0],marker=action[k][0],color="none",markerfacecolor=action[k][1],markersize=4,label=l) for k,l in [("INTRODUCE_OR_TUNE","Introduce/tune"),("RETAIN_FL","Retain FL"),("TIMING_SCREEN","Timing screen")]]+[Line2D([0],[0],marker="D",color="none",markerfacecolor="#8A5FB2",markersize=3.5,label="Haplotype screen")]
    leg3=[Line2D([0],[0],marker="s",color="none",markerfacecolor=geno[k],markeredgecolor="#222",markeredgewidth=.2,markersize=3.5,label=f"{k} = {v}") for k,v in [("D","Dura"),("P","Pisifera"),("M","Meizhou4")]]
    leg4=[Line2D([0],[0],marker="s",color="none",markerfacecolor=v,markeredgecolor="#222",markeredgewidth=.2,markersize=3.5,label=k) for k,v in actual.items()]
    for x,handles,title,ncol in [(0.02,leg1,"FL ancestry background",2),(0.29,leg2,"Favorable action",2),(0.58,leg3,"Diploid target chips",1),(0.78,leg4,"Actual-haplotype chips",2)]:
        la=fig.add_axes([x,.01,.20,.27]); la.axis("off"); la.legend(handles=handles,title=title,frameon=False,ncol=ncol,loc="upper left",borderaxespad=0,columnspacing=.5,handlelength=1.0,handletextpad=.3)
    save(fig,stage,"Fig3l")


def sha256(path: Path) -> str:
    h=hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda:f.read(1024*1024),b""): h.update(block)
    return h.hexdigest()


def main():
    style()
    missing=[str(x) for x in INPUTS if not x.is_file()]
    if missing: raise SystemExit("Missing inputs:\n"+"\n".join(missing))
    stage=Path(tempfile.mkdtemp(prefix=".fig3_panels_build_",dir=OUT))
    try:
        for fn in [panel_ab,panel_c,panel_d,panel_e,panel_f,panel_g,panel_h,panel_i,panel_j,panel_k,panel_l]: fn(stage)
        outputs=sorted(stage.glob("Fig3*.pdf"))+sorted(stage.glob("Fig3*_600dpi.png"))
        if len(outputs)!=22: raise RuntimeError(f"Expected 22 outputs, got {len(outputs)}")
        for p in outputs:
            dest=OUT/p.name
            if dest.exists(): raise FileExistsError(f"Refusing to overwrite {dest}")
        for p in outputs: os.replace(p,OUT/p.name)
    finally:
        shutil.rmtree(stage,ignore_errors=True)
    manifest=[]
    for stem,(w,h) in SIZES_MM.items():
        manifest.append(f"{stem}\t{w}\t{h}\t{OUT/(stem+'.pdf')}\t{OUT/(stem+'_600dpi.png')}")
    (OUT/"FIGURE3_PANEL_MANIFEST.tsv").write_text("panel\twidth_mm\theight_mm\tpdf\tpng_600dpi\n"+"\n".join(manifest)+"\n")
    lines=[f"{sha256(p)}  {p.name}" for p in sorted(OUT.glob("Fig3*.pdf"))+sorted(OUT.glob("Fig3*_600dpi.png"))]
    (OUT/"FIGURE3_PANEL_CHECKSUMS.sha256").write_text("\n".join(lines)+"\n")
    inlines=[f"{sha256(p)}\t{p}" for p in INPUTS]
    (OUT/"FIGURE3_PANEL_INPUTS.sha256.tsv").write_text("sha256\tpath\n"+"\n".join(inlines)+"\n")
    (OUT/"FIGURE3_PANEL_RUN_INFO.txt").write_text(
        f"command={sys.executable} {Path(__file__).resolve()}\npython={platform.python_version()}\nmatplotlib={mpl.__version__}\npandas={pd.__version__}\nnumpy={np.__version__}\nscipy={stats.__version__ if hasattr(stats,'__version__') else 'see scipy package'}\n")

if __name__=="__main__":
    main()
