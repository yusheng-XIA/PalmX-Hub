#!/usr/bin/env python3
"""Render 01A, 01C–01H and 02A–02D in the supplied reference styles."""

from __future__ import annotations

from pathlib import Path
import math
import sys

import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.lines import Line2D
from matplotlib.patches import Arc, Circle, PathPatch, Patch, Rectangle
from matplotlib.path import Path as MplPath
import numpy as np
import pandas as pd
from scipy import stats
from scipy.ndimage import gaussian_filter1d


OUT = Path(__file__).resolve().parent
ASE_OUT = Path("${ANALYSIS_DIR}/22_answer_reviews/00_ms/03_V3/03_figure3/01_ASE/00_shared/runs/RUN-ASE-HETEROSIS-DOWNSTREAM-V2-001/output")
CURRENT = Path("${ANALYSIS_DIR}/22_answer_reviews/00_ms/03_V3/02_figure/_runs/RUN-ASE-FIGURES-CURRENT-V3-001/output")

RED = "#F8766D"
BLUE = "#4EA5DF"
TEAL = "#8DD3C7"
YELLOW = "#FFF59D"
GRAY = "#B8BDC3"
DARK = "#333333"
LIGHT = "#E7E8EA"
CLASS_ORDER = ["NoDiff", "Sub", "HapDom", "NoASE"]
CLASS_COLOR = {"NoDiff": RED, "Sub": TEAL, "HapDom": YELLOW, "NoASE": BLUE}
STAGES = ["0d", "15d", "35d", "50d", "65d", "80d", "95d", "110d", "125d",
          "140d", "155d", "170d", "185d", "12h", "24h", "36h", "48h", "60h", "72h"]


def setup() -> None:
    mpl.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 7.5,
        "axes.labelsize": 8.0, "xtick.labelsize": 6.8, "ytick.labelsize": 6.8,
        "legend.fontsize": 6.4, "axes.linewidth": 0.75,
        "xtick.major.width": 0.75, "ytick.major.width": 0.75,
        "xtick.major.size": 3.0, "ytick.major.size": 3.0,
        "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
    })


def save(fig: plt.Figure, stem: str) -> None:
    for ext in ("pdf", "svg", "png"):
        kw = {"bbox_inches": "tight", "pad_inches": 0.045, "facecolor": "white"}
        if ext == "png":
            kw["dpi"] = 600
        fig.savefig(OUT / f"{stem}.{ext}", **kw)
    plt.close(fig)


def letter(fig: plt.Figure, s: str, x: float = 0.015, y: float = 0.985) -> None:
    fig.text(x, y, s, fontsize=12.5, fontweight="bold", ha="left", va="top")


def open_axes(ax: plt.Axes) -> None:
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(direction="out")


def ribbon(ax: plt.Axes, p0, p1, color: str, width: float = 0.035, alpha: float = 0.38) -> None:
    x0, y0 = p0; x1, y1 = p1
    c1 = (0.42*x0 + 0.58*x1, 0.56)
    c2 = (0.58*x0 + 0.42*x1, 0.56)
    verts = [(x0, y0-width), c1, c2, (x1, y1-width),
             (x1, y1+width), c2, c1, (x0, y0+width), (x0, y0-width)]
    codes = [MplPath.MOVETO, MplPath.CURVE4, MplPath.CURVE4, MplPath.CURVE4,
             MplPath.LINETO, MplPath.CURVE4, MplPath.CURVE4, MplPath.CURVE4, MplPath.CLOSEPOLY]
    ax.add_patch(PathPatch(MplPath(verts, codes), facecolor=color, edgecolor="none", alpha=alpha))


def panel_01a() -> None:
    classes = pd.read_csv(OUT / "source_B_ASE_class_proportions.tsv", sep="\t")
    totals = classes.groupby("analysis").genes.sum().to_dict()
    fig, ax = plt.subplots(figsize=(5.2, 4.35))
    fig.subplots_adjust(left=0.07, right=0.98, bottom=0.08, top=0.96)
    ax.set_aspect("equal"); ax.axis("off"); ax.set_xlim(-1.25, 1.25); ax.set_ylim(-1.08, 1.15)
    letter(fig, "a")

    ax.add_patch(Arc((0,0), 2.0, 2.0, theta1=88, theta2=272, lw=5.2, color=RED))
    ax.add_patch(Arc((0,0), 2.0, 2.0, theta1=-88, theta2=88, lw=5.2, color=BLUE))
    ax.text(-1.12, 0.05, "Allele A", rotation=90, ha="center", va="center", fontsize=8)
    ax.text(1.12, 0.05, "Allele B", rotation=-90, ha="center", va="center", fontsize=8)

    angles = np.deg2rad([145, 180, 220, 55, 10, -35])
    pts = [(math.cos(a), math.sin(a)) for a in angles]
    for i, (x, y) in enumerate(pts):
        col = RED if x < 0 else BLUE
        ax.scatter(x, y, marker="*", s=400, c=col, edgecolors="black", linewidths=1.1, zorder=5)
    ribbon(ax, pts[0], pts[4], "#A6D8F4", 0.045)
    ribbon(ax, pts[1], pts[3], "#A6D8F4", 0.045)
    ribbon(ax, pts[2], pts[5], "#A6D8F4", 0.045)

    # Same-CDS pair and excluded ambiguous gene, following the reference grammar.
    ax.scatter([-0.36, 0.36], [-0.93, -0.93], marker="*", s=360,
               c=[TEAL, TEAL], edgecolors="black", linewidths=1.0, zorder=5)
    ribbon(ax, (-0.36,-0.93), (0.36,-0.93), "#CDEEE7", 0.035, 0.60)
    ax.text(0, -0.72, "paired alleles with the same CDS", ha="center", fontsize=7.1)
    ax.scatter(-0.27, 0.90, marker="*", s=340, facecolors="white", edgecolors="black",
               linewidths=1.0, linestyles="dashed", zorder=6)
    ax.plot([-0.18,0.02],[0.96,1.04], color=RED, lw=5, solid_capstyle="butt")
    ax.text(0.05, 0.96, "unpaired/ambiguous gene excluded", va="center", fontsize=6.7)

    ax.text(0, 0.37, "Graph-based 1:1 allele pairs", ha="center", fontsize=8.5)
    ax.text(0, 0.18, f"FL   Africa hap2 ↔ American hap1   n = {totals['FL']:,}", ha="center", fontsize=7.4)
    ax.text(0, 0.03, f"TN   Dura/TK-like ↔ Pisifera/NS-like   n = {totals['TN']:,}", ha="center", fontsize=7.4)
    ax.text(0, -0.15, "replicate-aware ASE on eligible allele pairs", ha="center", color="#666", fontsize=6.6)
    save(fig, "01A_allele_framework")


def panel_01c() -> None:
    d = pd.read_csv(OUT / "source_C_intact_LTR_profiles.tsv", sep="\t")
    order = [("FL","Africa hap2","FL HapA"),("FL","American hap1","FL HapB"),
             ("TN","Dura-like","TN HapA"),("TN","Pisifera-like","TN HapB")]
    # Reference panel c uses three smooth gene-type curves and one common y scale.
    classes = ["HapDom", "Sub", "NoDiff"]
    line_color = {"HapDom": "#4EA5DF", "Sub": "#8DD3C7", "NoDiff": "#F8766D"}
    data_max = float(d[d.overall_class.isin(classes)].occupancy_percent.max())
    # Give the curves visible headroom and keep an explicit zero baseline.
    ymax = max(0.8, math.ceil((data_max + 0.06) * 10) / 10)
    # A taller canvas gives each of the four profiles a reference-like vertical aspect.
    fig, axes = plt.subplots(4, 1, figsize=(3.55, 5.55), sharex=True, sharey=True)
    fig.subplots_adjust(left=0.20, right=0.98, bottom=0.09, top=0.90, hspace=0.12)
    letter(fig, "c")
    for i, (ax, (analysis, hap, short_label)) in enumerate(zip(axes, order)):
        z = d[(d.analysis==analysis)&(d.haplotype==hap)]
        for cls in classes:
            q = z[z.overall_class==cls].sort_values("plot_bin")
            smoothed = gaussian_filter1d(q.occupancy_percent.to_numpy(float), sigma=1.6,
                                         mode="nearest")
            ax.plot(q.plot_bin, smoothed, color=line_color[cls], lw=1.35)
        ax.add_patch(Rectangle((35.5, 0.008), 72, ymax-0.016, fill=False, edgecolor=DARK,
                               linewidth=0.75, linestyle=(0, (4, 2)), clip_on=True))
        ax.text(0.025, 0.78, short_label, transform=ax.transAxes, fontsize=6.7)
        ax.set_xlim(0,143); ax.set_ylim(bottom=0,top=ymax); open_axes(ax)
        ax.set_yticks(np.arange(0, ymax + 0.001, 0.2))
        ax.grid(False)
        if i < 3: ax.tick_params(labelbottom=False)
    axes[-1].set_xticks([0,35.5,107.5,143],["−5.0", "TSS", "TES", "+5.0 kb"])
    fig.supylabel("Intact-LTR occupancy (%)", x=0.055, fontsize=8)
    handles=[Line2D([0],[0],color=line_color[c],lw=1.5,label=c) for c in classes]
    fig.legend(handles=handles,title="Gene type",frameon=False,ncol=3,loc="upper center",
               bbox_to_anchor=(0.60,0.995),columnspacing=1.0,handlelength=2.0)
    save(fig, "01C_intact_LTR_profiles")


def panel_01d() -> None:
    d = pd.read_csv(OUT / "source_D_stage_mirror.tsv", sep="\t")
    fig, ax = plt.subplots(figsize=(4.25, 3.15))
    fig.subplots_adjust(left=0.16,right=0.98,bottom=0.35,top=0.80)
    letter(fig,"d")
    x=np.arange(len(STAGES)); width=0.34
    for analysis,off,col in [("FL",-width/2,RED),("TN",width/2,BLUE)]:
        q=d[d.analysis==analysis].pivot(index="stage",columns="ase_call",values="genes").reindex(STAGES)
        ax.bar(x+off,q.Allele_A_biased,width,color=col,edgecolor=col,lw=0.85)
        ax.bar(x+off,-q.Allele_B_biased,width,facecolor="white",edgecolor=col,lw=1.05)
    ax.axhline(0,color=DARK,lw=0.7); ax.set_xlim(-0.6,18.6)
    lim=ax.get_ylim()[1]; ax.set_ylim(-lim,lim)
    ticks=ax.get_yticks(); ax.set_yticks(ticks); ax.set_yticklabels([f"{abs(int(t)):,}" for t in ticks])
    ax.set_ylabel("Number"); ax.set_xticks(x,STAGES,rotation=45,ha="right")
    open_axes(ax)
    ax.plot([.02,.67],[-.20,-.20],transform=ax.transAxes,color=TEAL,lw=4,clip_on=False)
    ax.plot([.70,.98],[-.20,-.20],transform=ax.transAxes,color=YELLOW,lw=4,clip_on=False)
    ax.text(.345,-.31,"Developmental",transform=ax.transAxes,ha="center",va="top",fontsize=6.3,clip_on=False)
    ax.text(.84,-.31,"Postharvest",transform=ax.transAxes,ha="center",va="top",fontsize=6.3,clip_on=False)
    handles=[Patch(facecolor=RED,label="FL"),Patch(facecolor=BLUE,label="TN"),
             Patch(facecolor=GRAY,edgecolor=GRAY,label="A > B"),Patch(facecolor="white",edgecolor=GRAY,label="A < B")]
    fig.legend(handles=handles,frameon=False,ncol=2,loc="upper left",bbox_to_anchor=(0.14,0.98),columnspacing=1.0)
    save(fig,"01D_stage_mirror_ASE")


def panel_01e() -> None:
    d=pd.read_csv(OUT/"source_E_FL_TN_alluvial.tsv",sep="\t")
    total=d.orthogroups.sum(); order=["NoDiff","HapDom","Sub","NoASE"]
    left=d.groupby("overall_class_FL").orthogroups.sum().reindex(order).fillna(0)
    right=d.groupby("overall_class_TN").orthogroups.sum().reindex(order).fillna(0)
    gap=0.012; usable=1-gap*(len(order)-1)
    def intervals(v):
        out={}; top=1
        for c in order:
            h=usable*v[c]/v.sum(); out[c]=(top-h,top); top-=h+gap
        return out
    li,ri=intervals(left),intervals(right)
    lcur={c:li[c][0] for c in order}; rcur={c:ri[c][0] for c in order}
    fig,ax=plt.subplots(figsize=(4.0,3.35)); fig.subplots_adjust(left=.13,right=.94,bottom=.12,top=.84)
    letter(fig,"e"); ax.set_xlim(0,1); ax.set_ylim(-.04,1.03); ax.axis("off")
    for row in d.sort_values(["overall_class_FL","overall_class_TN"]).itertuples():
        a,b,n=row.overall_class_FL,row.overall_class_TN,row.orthogroups
        h=usable*n/total; y0,y1=lcur[a],rcur[b]
        verts=[(.20,y0),(.45,y0),(.55,y1),(.80,y1),(.80,y1+h),(.55,y1+h),(.45,y0+h),(.20,y0+h),(.20,y0)]
        codes=[MplPath.MOVETO,MplPath.CURVE4,MplPath.CURVE4,MplPath.CURVE4,MplPath.LINETO,MplPath.CURVE4,MplPath.CURVE4,MplPath.CURVE4,MplPath.CLOSEPOLY]
        ax.add_patch(PathPatch(MplPath(verts,codes),facecolor=CLASS_COLOR[a],edgecolor="none",alpha=.40))
        lcur[a]+=h; rcur[b]+=h
    for x,vals,ints,label in [(.14,left,li,"FL"),(.80,right,ri,"TN")]:
        for c in order:
            lo,hi=ints[c]; ax.add_patch(Rectangle((x,lo),.07,hi-lo,facecolor=CLASS_COLOR[c],edgecolor=DARK,lw=.65))
        ax.text(x+.035,-.025,label,ha="center",va="top",fontsize=8)
    ax.text(.5,-.025,f"Reviewed 1:1 orthogroups, n = {total:,}",ha="center",va="top",fontsize=6.3,color="#666")
    handles=[Patch(facecolor=CLASS_COLOR[c],edgecolor=DARK,lw=.4,label=c) for c in order]
    fig.legend(handles=handles,title="Gene type",frameon=False,ncol=4,loc="upper center",bbox_to_anchor=(.57,.99),columnspacing=.9)
    save(fig,"01E_FL_TN_alluvial")


def kde(ax, values, color, ls="-", xmax=3.0):
    v=pd.to_numeric(values,errors="coerce").dropna().to_numpy(float); v=v[(v>=0)&(v<=xmax)]
    if len(v)>5 and np.std(v)>0:
        xx=np.linspace(0,xmax,260); ax.plot(xx,stats.gaussian_kde(v)(xx),color=color,lw=1.05,ls=ls)


def panel_01f() -> None:
    d=pd.read_csv(OUT/"source_F_KaKs_by_ASE_class.tsv",sep="\t")
    fig,ax=plt.subplots(figsize=(3.55,3.15)); fig.subplots_adjust(left=.16,right=.97,bottom=.16,top=.94)
    letter(fig,"f")
    for analysis,ls in [("FL","-"),("TN","--")]:
        for cls in ["HapDom","NoDiff","Sub"]: kde(ax,d[(d.analysis==analysis)&(d.ASE_type==cls)].Ka_Ks,CLASS_COLOR[cls],ls,3)
    ax.set_xlim(0,3); ax.set_xlabel(r"$K_a/K_s$ ratio"); ax.set_ylabel("Density"); open_axes(ax)
    ax.axvline(1,color="#555",lw=.6,ls=":")
    ins=ax.inset_axes([.53,.50,.44,.44])
    for analysis,ls in [("FL","-"),("TN","--")]:
        for cls in ["HapDom","NoDiff","Sub"]: kde(ins,d[(d.analysis==analysis)&(d.ASE_type==cls)].Ks,CLASS_COLOR[cls],ls,.10)
    ins.set_xlim(0,.10); ins.set_xlabel(r"$K_s$ value",fontsize=6.2); ins.set_ylabel("Density",fontsize=6.2); ins.tick_params(labelsize=5.5); open_axes(ins)
    handles=[Line2D([0],[0],color=CLASS_COLOR[c],lw=1.2,label=c) for c in ["HapDom","NoDiff","Sub"]]
    handles += [Line2D([0],[0],color=DARK,ls="-",label="FL"),Line2D([0],[0],color=DARK,ls="--",label="TN")]
    ax.legend(handles=handles,frameon=False,ncol=2,loc="lower right",bbox_to_anchor=(1.0,.035),fontsize=5.7)
    save(fig,"01F_KaKs")


def panel_01g() -> None:
    raw=pd.concat([pd.read_csv(CURRENT/"FL_diagnostic_SNP_density_current.tsv.gz",sep="\t"),pd.read_csv(CURRENT/"TN_diagnostic_SNP_density_current.tsv.gz",sep="\t")])
    fig,axes=plt.subplots(1,2,figsize=(5.35,3.12),sharey=True); fig.subplots_adjust(left=.12,right=.98,bottom=.20,top=.80,wspace=.12)
    letter(fig,"g")
    regions=[("upstream2kb_per_kb","Upstream\n2 kb"),("gene_body_per_kb","Gene"),("downstream2kb_per_kb","Downstream\n2 kb")]
    cats=["HapDom","NoDiff","Sub"]
    tests=[]

    def p_text(p: float) -> str:
        if p < 1e-300:
            return r"$P < 10^{-300}$"
        if p < 0.001:
            exponent=int(math.floor(math.log10(p)))
            coefficient=p/(10**exponent)
            return rf"$P = {coefficient:.1f}\times10^{{{exponent}}}$"
        return rf"$P = {p:.3f}$"

    def holm_adjust(pvals: list[float]) -> np.ndarray:
        p=np.asarray(pvals,float); order=np.argsort(p); adjusted=np.empty(len(p)); running=0.0
        for rank,idx in enumerate(order):
            running=max(running,min(1.0,p[idx]*(len(p)-rank)))
            adjusted[idx]=running
        return adjusted

    def compact_letters(sig: dict[tuple[int,int],bool], medians: list[float]) -> list[str]:
        n_sig=sum(sig.values()); ranked=sorted(range(3),key=lambda i:medians[i],reverse=True)
        if n_sig==0:
            return ["a","a","a"]
        if n_sig==3:
            out=[""]*3
            for label,idx in zip(["a","b","c"],ranked): out[idx]=label
            return out
        if n_sig==1:
            i,j=next(pair for pair,is_sig in sig.items() if is_sig); third=({0,1,2}-{i,j}).pop()
            high,low=(i,j) if medians[i]>=medians[j] else (j,i)
            out=[""]*3; out[high]="a"; out[low]="b"; out[third]="ab"
            return out
        # With three groups, two significant pairs leave one non-significant pair.
        i,j=next(pair for pair,is_sig in sig.items() if not is_sig); isolated=({0,1,2}-{i,j}).pop()
        pair_letter="b" if medians[isolated]>=np.mean([medians[i],medians[j]]) else "a"
        isolated_letter="a" if pair_letter=="b" else "b"
        out=[pair_letter]*3; out[isolated]=isolated_letter
        return out

    for ax,analysis in zip(axes,["FL","TN"]):
        q=raw[raw.analysis==analysis]; data=[]; pos=[]; cols=[]
        for j,(col,_) in enumerate(regions):
            region_values=[]
            for k,c in enumerate(cats):
                values=np.log10(1+q.loc[q.overall_class==c,col].dropna().to_numpy(float))
                data.append(values); region_values.append(values); pos.append(j*4+k); cols.append(CLASS_COLOR[c])
            kw=stats.kruskal(*region_values)
            pairs=[(0,1),(0,2),(1,2)]
            raw_p=[stats.mannwhitneyu(region_values[a],region_values[b],alternative="two-sided").pvalue for a,b in pairs]
            adj_p=holm_adjust(raw_p); sig={pair:(p<0.05) for pair,p in zip(pairs,adj_p)}
            medians=[float(np.median(v)) for v in region_values]; cld=compact_letters(sig,medians)
            ax.text(j*4+1,2.49,p_text(float(kw.pvalue)),ha="center",va="bottom",fontsize=5.4)
            for k,label in enumerate(cld):
                ax.text(j*4+k,2.27,label,ha="center",va="bottom",fontsize=6.8,fontweight="bold")
            for (a,b),rp,ap in zip(pairs,raw_p,adj_p):
                tests.append({"analysis":analysis,"region":col,"test":"Mann-Whitney U",
                              "comparison":f"{cats[a]} vs {cats[b]}","pvalue":rp,
                              "holm_pvalue":ap,"significant_0.05":ap<0.05,
                              "letters":";".join(f"{c}:{l}" for c,l in zip(cats,cld))})
            tests.append({"analysis":analysis,"region":col,"test":"Kruskal-Wallis",
                          "comparison":"all three ASE classes","pvalue":float(kw.pvalue),
                          "holm_pvalue":np.nan,"significant_0.05":kw.pvalue<0.05,
                          "letters":";".join(f"{c}:{l}" for c,l in zip(cats,cld))})
        bp=ax.boxplot(data,positions=pos,widths=.72,showfliers=False,patch_artist=True,
                      medianprops=dict(color=DARK,lw=.8),whiskerprops=dict(color="#555",lw=.55),capprops=dict(color="#555",lw=.55))
        for b,c in zip(bp["boxes"],cols): b.set_facecolor("white"); b.set_edgecolor(c); b.set_linewidth(1.0)
        ax.set_xticks([1,5,9],[x[1] for x in regions]); ax.set_title(analysis,fontsize=8,fontweight="normal",pad=10); open_axes(ax)
        ax.set_ylim(-.10,2.68)
    axes[0].set_ylabel(r"log$_{10}$(1 + SNPs per 1,000 bp)")
    handles=[Patch(facecolor="white",edgecolor=CLASS_COLOR[c],lw=1,label=c) for c in cats]
    fig.legend(handles=handles,frameon=False,ncol=3,loc="upper center",bbox_to_anchor=(.56,.99))
    pd.DataFrame(tests).to_csv(OUT/"source_G_SNP_density_tests.tsv",sep="\t",index=False)
    save(fig,"01G_SNP_density")


def panel_01h() -> None:
    # Candidates are selected to represent the observed cultivar phenotypes, then
    # ranked within each relevant family by robust ASE and directional consistency.
    trait=pd.read_csv(ASE_OUT/"trait_haplotype_ASE.tsv.gz",sep="\t")
    existing=pd.read_csv(OUT/"source_H_candidate_trajectories.tsv",sep="\t")
    specs=[
        ("FL","evm.TU.chr09B.1430","MADS34","Seedlessness"),
        ("FL","evm.TU.chr08B.1796","SAD","High oleic acid"),
        ("FL","evm.TU.chr03B.1947","VTE1","Oxidative stability"),
        ("TN","evm.TU.chr10.635","DGAT","High oil content"),
        ("TN","evm.TU.chr01.3309","SAD","Low-oleic context"),
        ("TN","evm.TU.chr12.1328","LOX9","Rancidity susceptibility"),
    ]
    selected=[]; summaries=[]
    for analysis,gene_id,label,phenotype in specs:
        if label=="MADS34":
            q=existing[(existing.analysis==analysis)&(existing.gene_id==gene_id)].copy()
        else:
            q=trait[(trait.analysis==analysis)&(trait.gene_id==gene_id)].copy()
            q=q.drop_duplicates(["analysis","gene_id","stage"]).sort_values("stage_index")
            q["allele_A_fraction"]=q.ref_fraction
            q["allele_B_fraction"]=1-q.ref_fraction
        q=q.sort_values("stage_index").copy(); q["display_label"]=label; q["phenotype_module"]=phenotype
        selected.append(q)
        eligible=q[q.eligible.astype(bool)]
        median_ratio=float(eligible.log2_allele_ratio.median())
        dominant="Allele A" if median_ratio>0 else "Allele B"
        robust=eligible[eligible.robust_ase.astype(bool)]
        n_a=int((robust.ase_call=="Allele_A_biased").sum())
        n_b=int((robust.ase_call=="Allele_B_biased").sum())
        if n_a>=3 and n_b>=3:
            trajectory_label=f"Bias switching (A {n_a}; B {n_b})"
        elif n_a>=n_b:
            trajectory_label=f"A-biased ASE ({n_a}/{len(robust)} robust stages)"
        else:
            trajectory_label=f"B-biased ASE ({n_b}/{len(robust)} robust stages)"
        summaries.append({"analysis":analysis,"gene_id":gene_id,"display_label":label,
                          "phenotype_module":phenotype,"eligible_stages":len(eligible),
                          "robust_ASE_stages":len(robust),"robust_A_stages":n_a,
                          "robust_B_stages":n_b,"median_log2_A_over_B":median_ratio,
                          "dominant_haplotype":dominant,"trajectory_label":trajectory_label})
    d=pd.concat(selected,ignore_index=True)
    keep=[c for c in ["analysis","gene_id","display_label","phenotype_module","stage","stage_index",
                       "stage_group","allele_A_fraction","allele_B_fraction","log2_allele_ratio","padj",
                       "eligible","robust_ase","ase_call","allele_A_label","allele_B_label"] if c in d.columns]
    d[keep].to_csv(OUT/"source_H_phenotype_candidate_trajectories.tsv",sep="\t",index=False)
    summary=pd.DataFrame(summaries)
    summary.to_csv(OUT/"source_H_phenotype_candidate_summary.tsv",sep="\t",index=False)

    fig,axes=plt.subplots(2,3,figsize=(7.15,4.15),sharex=True,sharey=True)
    fig.subplots_adjust(left=.09,right=.99,bottom=.16,top=.83,wspace=.18,hspace=.50)
    letter(fig,"h")
    x=np.arange(19)
    for ax,(analysis,gene_id,gene,phenotype) in zip(axes.flat,specs):
        q=d[(d.analysis==analysis)&(d.gene_id==gene_id)].sort_values("stage_index")
        trajectory_label=summary.loc[(summary.analysis==analysis)&(summary.gene_id==gene_id),"trajectory_label"].iat[0]
        ax.plot(x,q.allele_A_fraction,"o-",color=RED,ms=2.5,lw=1.0)
        ax.plot(x,q.allele_B_fraction,"o-",color=BLUE,ms=2.5,lw=1.0)
        ax.axhline(.5,color="#666",lw=.55,ls="--"); ax.set_ylim(-.03,1.03); ax.set_xlim(-.5,18.5)
        ax.set_title(f"$\\it{{{gene}}}$ ({analysis})  |  {phenotype}\n{trajectory_label}",
                     fontsize=6.8,fontweight="normal",pad=4)
        ax.set_xticks([0,4,8,12,15,18],[STAGES[i] for i in [0,4,8,12,15,18]],rotation=45,ha="right")
        open_axes(ax)
        ax.plot([0,12],[-.12,-.12],color=TEAL,lw=3,clip_on=False); ax.plot([13,18],[-.12,-.12],color=YELLOW,lw=3,clip_on=False)
    axes[0,0].set_ylabel("Allelic fraction"); axes[1,0].set_ylabel("Allelic fraction")
    fig.legend([Line2D([0],[0],color=RED,marker="o",lw=1),Line2D([0],[0],color=BLUE,marker="o",lw=1)],
               ["Allele A","Allele B"],title="Haplotype contribution",frameon=False,ncol=2,
               loc="upper center",bbox_to_anchor=(.55,.995))
    save(fig,"01H_candidate_trajectories")


def panel_02a() -> None:
    d=pd.read_csv(OUT/"source2_A_ASE_genes_per_1Mb.tsv",sep="\t")
    maxmb=int(d.window_mb.max())+1; vmax=float(d.ASE_genes.quantile(.98))
    cm_fl=LinearSegmentedColormap.from_list("fl",["#F2F7EB","#73B943","#2E6B45"])
    cm_tn=LinearSegmentedColormap.from_list("tn",["#F8F3EF","#ED9B7B","#D9543D"])
    fig,ax=plt.subplots(figsize=(5.9,6.7)); fig.subplots_adjust(left=.11,right=.92,bottom=.09,top=.91)
    letter(fig,"A")
    for chrom in range(1,17):
        y=16-chrom; ax.text(-4.5,y,f"Chr{chrom}",ha="right",va="center",fontsize=6.8)
        for analysis,dy,cmap in [("FL",.13,cm_fl),("TN",-.13,cm_tn)]:
            q=d[(d.analysis==analysis)&(d.chrom_num==chrom)]; lookup=dict(zip(q.window_mb,q.ASE_genes)); length=int(q.window_mb.max())+1 if len(q) else 0
            for mb in range(length): ax.add_patch(Rectangle((mb,y+dy-.10),1,.19,facecolor=cmap(min(1,lookup.get(mb,0)/vmax)),edgecolor="none"))
    ax.set_xlim(-6,maxmb); ax.set_ylim(-.65,15.65); ax.set_yticks([]); ax.xaxis.tick_top(); ax.xaxis.set_label_position("top")
    ax.set_xlabel("The number of ASE genes within 1-Mb windows",labelpad=8,fontsize=9)
    ax.set_xticks(np.linspace(0,maxmb,8),[f"{int(v)}Mb" for v in np.linspace(0,maxmb,8)])
    ax.spines[:].set_visible(False); ax.tick_params(top=True,bottom=False,direction="out")
    ax.legend([Patch(color="#73B943"),Patch(color="#ED9B7B")],["FL","TN"],frameon=False,ncol=2,loc="lower right")
    save(fig,"02A_chromosome_1Mb")


def panel_02b() -> None:
    p=pd.read_csv(OUT/"source2_B_upset_intersections.tsv",sep="\t").head(18)
    sizes=pd.read_csv(OUT/"source2_B_upset_set_sizes.tsv",sep="\t")
    labels=["FL Early","FL Mid","FL Late","FL Post","TN Early","TN Mid","TN Late","TN Post"]
    fig=plt.figure(figsize=(4.55,3.75)); letter(fig,"B")
    gs=fig.add_gridspec(2,2,width_ratios=[.31,.69],height_ratios=[.48,.52],left=.18,right=.98,bottom=.12,top=.91,wspace=.05,hspace=.04)
    blank=fig.add_subplot(gs[0,0]); blank.axis("off"); top=fig.add_subplot(gs[0,1]); mat=fig.add_subplot(gs[1,1],sharex=top); left=fig.add_subplot(gs[1,0],sharey=mat)
    x=np.arange(len(p)); top.bar(x,p.orthogroups,color=DARK,width=.68); top.set_ylabel("Intersection Size"); top.set_xticks([]); open_axes(top)
    for j,pat in enumerate(p.pattern.astype(str)):
        bits=[int(v) for v in pat.zfill(8)]; on=[i for i,b in enumerate(bits) if b]
        for i,b in enumerate(bits): mat.scatter(j,i,s=10,color=DARK if b else "#E0E1E2")
        if len(on)>1: mat.plot([j,j],[min(on),max(on)],color=DARK,lw=.7)
    mat.set_yticks(range(8)); mat.set_yticklabels([]); mat.invert_yaxis(); mat.set_xticks([]); mat.spines[:].set_visible(False)
    ss=sizes.set_index("set").reindex(["FL_Early","FL_Mid","FL_Late","FL_Postharvest","TN_Early","TN_Mid","TN_Late","TN_Postharvest"])
    if ss.orthogroups.isna().any(): ss=sizes.set_index("set").iloc[:8]
    left.barh(range(8),ss.orthogroups.to_numpy(),color=["#9E9AC8"]*4+["#C6E9EF"]*4,height=.55)
    left.set_yticks(range(8),labels); left.invert_xaxis(); left.set_xlabel("Set Size"); left.spines[["top","right","left"]].set_visible(False)
    save(fig,"02B_phase_upset")


def panel_02c() -> None:
    d=pd.read_csv(OUT/"source2_C_volcano_95d.tsv",sep="\t")
    fig,axes=plt.subplots(1,2,figsize=(4.55,3.05),sharey=True); fig.subplots_adjust(left=.13,right=.98,bottom=.16,top=.83,wspace=.08)
    letter(fig,"C")
    for ax,analysis in zip(axes,["FL","TN"]):
        q=d[d.analysis==analysis].copy(); q["y"]=-np.log10(q.padj.clip(lower=1e-300)).clip(upper=82)
        cols=np.where(q.ase_call.eq("Allele_A_biased"),"#DCEFF1",np.where(q.ase_call.eq("Allele_B_biased"),"#D8D7E8","#ECEDEE"))
        ax.scatter(q.log2_allele_ratio,q.y,s=3,c=cols,alpha=.75,rasterized=True,edgecolors="none")
        robust=q[q.robust_ase & (q.y>5)].copy()
        picked=[]
        for target in (-8,-5,-2,2,5,8):
            if len(robust): picked.append((robust.log2_allele_ratio-target).abs().idxmin())
        hi=robust.loc[picked].drop_duplicates("gene_id") if picked else robust.head(0)
        ax.scatter(hi.log2_allele_ratio,hi.y,s=28,c="#D7191C",edgecolors="white",
                   linewidths=.45,zorder=10)
        for r in hi.itertuples(): ax.annotate(r.gene_id.split(".")[-1],(r.log2_allele_ratio,r.y),xytext=(2,2),textcoords="offset points",fontsize=4.8)
        ax.axvline(-.5,color=DARK,lw=.55,ls=":"); ax.axvline(.5,color=DARK,lw=.55,ls=":"); ax.axhline(-math.log10(.05),color=DARK,lw=.55,ls=":")
        ax.set_xlim(-11,11); ax.set_ylim(0,90); ax.set_xlabel(r"log$_2$(FoldChange)"); open_axes(ax)
        ax.add_patch(Rectangle((0,1.0),1,.11,transform=ax.transAxes,facecolor="#D9D9D9",edgecolor=DARK,lw=.65,clip_on=False))
        ax.text(.5,1.055,analysis,transform=ax.transAxes,ha="center",va="center",fontweight="bold",fontsize=7)
    axes[0].set_ylabel(r"−log$_{10}$(FDR)")
    save(fig,"02C_volcano_95d")


def panel_02d() -> None:
    # The comprehensive renderer includes FAD2 and every eligible gene in the
    # named oil-production, unsaturation and rancidity families.
    from render_trait_function_panels import main as render_trait_function_panels
    render_trait_function_panels()


def main() -> None:
    setup()
    renderers={
        "01a":panel_01a,"01c":panel_01c,"01d":panel_01d,"01e":panel_01e,
        "01f":panel_01f,"01g":panel_01g,"01h":panel_01h,"02a":panel_02a,
        "02b":panel_02b,"02c":panel_02c,"02d":panel_02d,
    }
    targets=[x.lower() for x in sys.argv[1:]] or list(renderers)
    unknown=[x for x in targets if x not in renderers]
    if unknown:
        raise SystemExit(f"Unknown panel(s): {', '.join(unknown)}")
    for target in targets:
        renderers[target]()


if __name__ == "__main__":
    main()
