#!/usr/bin/env python3
"""Render individual FL/TN ASE, cis/trans and heterosis panels.

Every panel is exported as PDF, editable SVG, 600-dpi PNG and a source TSV.
The script intentionally keeps all outputs in this directory (no subfolders).
"""

from __future__ import annotations

from pathlib import Path
import math
import warnings

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize, TwoSlopeNorm
from matplotlib.lines import Line2D
from matplotlib.patches import FancyBboxPatch, Rectangle
import numpy as np
import pandas as pd
from scipy.stats import chi2_contingency, fisher_exact


OUT = Path(__file__).resolve().parent
SHARED = OUT.parent / "01_ASE" / "00_shared"
ASE_RUN = SHARED / "runs" / "RUN-ASE-HETEROSIS-DOWNSTREAM-V2-001"
HET_RUN = SHARED / "runs" / "RUN-TN-EXPR-HETEROSIS-V2-001"

STAGES = ["0d", "15d", "35d", "50d", "65d", "80d", "95d", "110d", "125d",
          "140d", "155d", "170d", "185d", "12h", "24h", "36h", "48h", "60h", "72h"]
PHASES = ["Days 0–65", "Days 80–140", "Days 155–185", "Hours 12–72"]
REG_ORDER = ["I.Cis_only", "II.Trans_only", "III.Cis_trans_enhancing",
             "IV.Cis_trans_compensating", "V.Compensatory", "VI.Conserved", "VII.Ambiguous"]
REG_LABELS = ["I  Cis only", "II  Trans only", "III  Cis + trans\nenhancing",
              "IV  Cis + trans\ncompensating", "V  Compensatory", "VI  Conserved", "VII  Ambiguous"]
REG_COLORS = ["#E76F51", "#F4A261", "#2A9D8F", "#5AB4AC", "#457B9D", "#4A5568", "#B7BEC7"]

# Figure 3 palette harmonized to the favorable D/P/M breeding blueprint.
DPM_BLUE = "#3767A6"
DPM_RED = "#E56565"
DPM_GOLD = "#D8B64C"
DPM_TEAL = "#159D82"
DPM_PURPLE = "#8A5FB2"
DPM_HEATMAP_CMAP = LinearSegmentedColormap.from_list(
    "dpm_blue_gold_red",
    [DPM_BLUE, "#8FAED0", "#F7F5EE", DPM_GOLD, DPM_RED],
)


def style() -> None:
    mpl.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 8.0,
        "axes.labelsize": 8.5, "axes.titlesize": 9.5,
        "xtick.labelsize": 7.2, "ytick.labelsize": 7.2,
        "legend.fontsize": 7.2, "axes.linewidth": 0.9,
        "xtick.major.width": 0.9, "ytick.major.width": 0.9,
        "xtick.major.size": 3.0, "ytick.major.size": 3.0,
        "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
        "savefig.facecolor": "white", "figure.facecolor": "white",
    })


def clean(ax, left=True, bottom=True):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(left)
    ax.spines["bottom"].set_visible(bottom)


def panel_letter(fig, letter: str, x=0.012, y=0.985):
    fig.text(x, y, letter, ha="left", va="top", fontsize=16, fontweight="bold")


def save(fig, stem: str):
    fig.savefig(OUT / f"{stem}.pdf", bbox_inches="tight")
    fig.savefig(OUT / f"{stem}.svg", bbox_inches="tight")
    fig.savefig(OUT / f"{stem}.png", dpi=600, bbox_inches="tight")
    plt.close(fig)


def write_source(df: pd.DataFrame, stem: str):
    df.to_csv(OUT / f"source_{stem}.tsv", sep="\t", index=False)


def chromosome_panel(analysis: str, letter: str, stem: str):
    src = pd.read_csv(OUT / "source2_A_ASE_genes_per_1Mb.tsv", sep="\t")
    d = src.loc[src.analysis.eq(analysis)].copy()
    d["DBA_genes"] = d["eligible_genes"]
    d["DEBA_genes"] = d["ASE_genes"]
    d["coordinate_anchor"] = "Africa hap2" if analysis == "FL" else "Dura/TK-like"
    d["DBA_definition"] = "diagnostic bi-allelic gene-pair candidates"
    d["DEBA_definition"] = "robust ASE subset of DBA candidates"
    write_source(d, stem)

    green = LinearSegmentedColormap.from_list("dba", ["#F2F2F2", "#B7D77B", "#2C6E5B"])
    salmon = LinearSegmentedColormap.from_list("deba", ["#F2F2F2", "#F2B39A", "#E76F51"])
    vmax_dba = max(1, np.nanpercentile(d.DBA_genes, 99))
    vmax_deba = max(1, np.nanpercentile(d.DEBA_genes, 99))
    norm_dba, norm_deba = Normalize(0, vmax_dba), Normalize(0, vmax_deba)
    chr_ids = sorted(d.chrom_num.unique())
    max_mb = int(d.window_mb.max() + 1)

    fig, ax = plt.subplots(figsize=(7.2, 6.1))
    y_positions = np.arange(len(chr_ids))[::-1]
    for y, chrom in zip(y_positions, chr_ids):
        sub = d[d.chrom_num.eq(chrom)].sort_values("window_mb")
        chr_end = int(sub.window_mb.max() + 1)
        ax.add_patch(Rectangle((0, y + 0.09), chr_end, 0.26, facecolor="#F3F3F3", edgecolor="#D0D0D0", lw=0.35))
        ax.add_patch(Rectangle((0, y - 0.35), chr_end, 0.26, facecolor="#F3F3F3", edgecolor="#D0D0D0", lw=0.35))
        for row in sub.itertuples():
            ax.add_patch(Rectangle((row.window_mb, y + 0.09), 1.02, 0.26,
                                   facecolor=green(norm_dba(row.DBA_genes)), edgecolor="none"))
            ax.add_patch(Rectangle((row.window_mb, y - 0.35), 1.02, 0.26,
                                   facecolor=salmon(norm_deba(row.DEBA_genes)), edgecolor="none"))
    ax.set_xlim(0, max_mb)
    ax.set_ylim(-0.75, len(chr_ids) - 0.25)
    ax.set_yticks(y_positions)
    ax.set_yticklabels([f"Chr{int(c)}" for c in chr_ids])
    ax.xaxis.set_ticks_position("top")
    ax.xaxis.set_label_position("top")
    ax.set_xlabel("Chromosomal position (Mb)", labelpad=7)
    ax.tick_params(axis="x", top=True, labeltop=True, bottom=False, labelbottom=False)
    clean(ax, left=False, bottom=False)
    ax.tick_params(axis="y", length=0, pad=6)
    ax.set_title(f"{analysis}: DBA and DEBA genes within 1-Mb windows", pad=28, fontweight="bold")
    ax.text(1.0, -0.055, f"Coordinate anchor: {d.coordinate_anchor.iloc[0]}", transform=ax.transAxes,
            ha="right", va="top", color="#555555", fontsize=6.7)

    cax1 = fig.add_axes([0.765, 0.105, 0.018, 0.13])
    cb1 = mpl.colorbar.ColorbarBase(cax1, cmap=green, norm=norm_dba, orientation="vertical")
    cb1.set_label("DBA genes", fontsize=6.7); cb1.ax.tick_params(labelsize=6)
    cax2 = fig.add_axes([0.855, 0.105, 0.018, 0.13])
    cb2 = mpl.colorbar.ColorbarBase(cax2, cmap=salmon, norm=norm_deba, orientation="vertical")
    cb2.set_label("DEBA genes", fontsize=6.7); cb2.ax.tick_params(labelsize=6)
    panel_letter(fig, letter)
    fig.subplots_adjust(left=0.105, right=0.97, top=0.87, bottom=0.07)
    save(fig, stem)


def decision_tree():
    rows = [
        ("A", "log₂(TK / NS) in parental expression", "between-parent total-expression contrast"),
        ("B", "log₂(TK-like / NS-like allele) in TN F₁", "allelic contrast in replicated F₁ ASE"),
        ("A − B", "trans component", "residual parental contrast after cis contribution"),
        ("A cutoff", "|A| ≥ 0.75", "exploratory effect threshold"),
        ("B cutoff", "|B| ≥ 0.50", "F₁ ASE effect threshold"),
        ("A − B cutoff", "|A − B| ≥ 0.75", "exploratory trans-effect threshold"),
    ]
    write_source(pd.DataFrame(rows, columns=["parameter", "display", "definition"]), "04A_TN_cis_trans_decision_tree")
    fig, ax = plt.subplots(figsize=(7.2, 4.7)); ax.set_axis_off()
    ax.text(0.02, 0.95, r"$A=\log_2(TK/NS)$", fontsize=9.5)
    ax.text(0.02, 0.88, r"$B=\log_2(TK\!\!-like/NS\!\!-like)$ in TN F$_1$", fontsize=9.5)
    ax.text(0.02, 0.81, r"Trans component: $A-B$", fontsize=9.5)

    def box(x, y, w, h, text, fc="#FFFFFF", ec="#AAB2B8", fs=8.3, weight="normal"):
        p = FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.012,rounding_size=0.016",
                           fc=fc, ec=ec, lw=1.0)
        ax.add_patch(p); ax.text(x+w/2, y+h/2, text, ha="center", va="center", fontsize=fs, fontweight=weight)
        return (x, y, w, h)

    def connect(a, b):
        ax.plot([a[0]+a[2], b[0]], [a[1]+a[3]/2, b[1]+b[3]/2], color="#6B7075", lw=0.8, zorder=0)

    root1 = box(.03, .53, .12, .09, r"$|A|\geq0.75$", "#405E91", "#405E91", 8.5, "bold")
    root0 = box(.03, .22, .12, .09, r"$|A|<0.75$", "#FFFFFF")
    b1 = box(.23, .60, .13, .08, r"$|B|\geq0.50$", "#8799BD")
    b0 = box(.23, .46, .13, .08, r"$|B|<0.50$", "#FFFFFF")
    b1b = box(.23, .27, .13, .08, r"$|B|\geq0.50$", "#8799BD")
    b0b = box(.23, .13, .13, .08, r"$|B|<0.50$", "#FFFFFF")
    for z in (b1,b0): connect(root1,z)
    for z in (b1b,b0b): connect(root0,z)

    c1 = box(.44, .69, .16, .075, r"$|A-B|<0.75$", "#E0E0E0")
    c2 = box(.44, .59, .16, .075, "$|A-B|\\geq0.75$\n$B(A-B)>0$")
    c3 = box(.44, .49, .16, .075, "$|A-B|\\geq0.75$\n$B(A-B)<0$")
    c4 = box(.44, .39, .16, .075, r"$|A-B|\geq0.75$")
    c5 = box(.44, .27, .16, .075, r"$|A-B|\geq0.75$")
    c6 = box(.44, .15, .16, .075, r"$|A-B|<0.75$", "#E0E0E0")
    c7 = box(.44, .05, .16, .065, "Other")
    for z in (c1,c2,c3): connect(b1,z)
    connect(b0,c4); connect(b1b,c5); connect(b0b,c6)
    labels = ["I  Cis only", "III  Cis + trans enhancing", "IV  Cis + trans compensating",
              "II  Trans only", "V  Compensatory", "VI  Conserved", "VII  Ambiguous"]
    colors = [REG_COLORS[0],REG_COLORS[2],REG_COLORS[3],REG_COLORS[1],REG_COLORS[4],REG_COLORS[5],REG_COLORS[6]]
    for z, lab, col in zip((c1,c2,c3,c4,c5,c6,c7), labels, colors):
        ax.text(.65, z[1]+z[3]/2, lab, va="center", ha="left", fontsize=8.1, color=col, fontweight="bold")
    ax.text(.02, .005, "Exploratory parental-expression layer (TK and NS: n = 1 per stage); TN F₁ ASE is replicate-supported.",
            ha="left", va="bottom", fontsize=6.7, color="#555555")
    ax.set_xlim(0,1); ax.set_ylim(0,1)
    panel_letter(fig, "A")
    fig.suptitle("TN cis/trans regulatory classification", x=.53, y=.985, fontsize=10.5, fontweight="bold")
    fig.subplots_adjust(left=.03,right=.99,top=.91,bottom=.04)
    save(fig, "04A_TN_cis_trans_decision_tree")


def regulatory_bubble():
    df = pd.read_csv(ASE_RUN / "output" / "TN_cis_trans_classified.tsv", sep="\t")
    tab = (df.groupby(["stage_index","stage","regulatory_class"], observed=True).size().rename("genes").reset_index())
    totals = df.groupby(["stage_index","stage"], observed=True).size().rename("stage_total").reset_index()
    tab = tab.merge(totals, on=["stage_index","stage"])
    tab["percentage"] = 100 * tab.genes / tab.stage_total
    tab["inference_level"] = "exploratory_parent_n1"
    write_source(tab, "04B_TN_regulatory_classes_19stages")
    piv = tab.pivot(index="regulatory_class", columns="stage", values="percentage").reindex(index=REG_ORDER, columns=STAGES).fillna(0)
    maxp = float(piv.to_numpy().max())
    display_max = max(55.0, 5.0 * math.ceil(maxp / 5.0))

    # Matplotlib scatter ``s`` is marker area (pt^2).  Keep this mapping in
    # one function so the data bubbles and size legend are quantitatively
    # identical.
    def bubble_size(pct):
        arr = np.clip(np.asarray(pct, dtype=float), 0.0, display_max)
        sizes = 16.0 + 410.0 * arr / display_max
        return float(sizes) if np.ndim(sizes) == 0 else sizes

    fig, ax = plt.subplots(figsize=(8.65, 3.85))
    fig.subplots_adjust(left=.205, right=.795, bottom=.255, top=.835)
    xs, ys, vals = [], [], []
    for yi, reg in enumerate(REG_ORDER[::-1]):
        for xi, st in enumerate(STAGES):
            xs.append(xi); ys.append(yi); vals.append(piv.loc[reg,st])
    vals = np.asarray(vals)
    sc = ax.scatter(xs, ys, s=bubble_size(vals), c=vals, cmap="Reds",
                    vmin=0, vmax=display_max, edgecolor="white", linewidth=.45)
    ax.set_xticks(range(19)); ax.set_xticklabels(STAGES, rotation=45, ha="right")
    ax.set_yticks(range(7)); ax.set_yticklabels(REG_LABELS[::-1])
    ax.set_xlim(-.7,18.7); ax.set_ylim(-.65,6.65)
    ax.axvline(12.5, color="#777777", lw=.8, ls="--")
    ax.text(6,6.47,"Seed development", ha="center", va="bottom", fontsize=7, color="#3D7A69")
    ax.text(15.5,6.47,"Postharvest", ha="center", va="bottom", fontsize=7, color="#B68A25")
    clean(ax)
    ax.set_xlabel("Stage"); ax.set_ylabel("Regulatory class")
    cax = fig.add_axes([.830, .320, .017, .500])
    cbar = fig.colorbar(sc, cax=cax)
    cbar.set_label("Genes within stage (%)", fontsize=8.0, labelpad=7)
    cbar.set_ticks([v for v in [0, 10, 20, 30, 40, 50] if v <= display_max])
    cbar.ax.tick_params(labelsize=7.1, width=.75, length=3.0, pad=2.5)
    cbar.outline.set_linewidth(.8)

    # Dedicated axis makes the area legend stable across PDF/SVG/PNG exports.
    lax = fig.add_axes([.805, .075, .175, .175])
    lax.set(xlim=(0, 1), ylim=(0, 1)); lax.axis("off")
    lax.text(.5, .94, "Bubble size (%)", ha="center", va="top",
             fontsize=7.7, fontweight="medium", color="#222222")
    size_levels = [10, 30, 50]
    legend_x = [.18, .50, .82]
    lax.scatter(legend_x, [.55] * 3, s=[bubble_size(v) for v in size_levels],
                facecolor="#D9D9D9", edgecolor="#666666", linewidth=.55, zorder=3)
    for x, value in zip(legend_x, size_levels):
        lax.text(x, .16, f"{value}%", ha="center", va="center",
                 fontsize=7.0, color="#333333")

    fig.text(.795, .060, "Exploratory: TK/NS parental expression n = 1 per stage",
             ha="right", va="bottom", fontsize=6.35, color="#555555")
    panel_letter(fig,"B")
    ax.set_title("TN regulatory-class landscape across 19 stages", fontweight="bold", pad=10)
    save(fig,"04B_TN_regulatory_classes_19stages")


def clustered_bootstrap(values_by_gene, rng, nboot=500):
    keys = list(values_by_gene)
    if len(keys) < 2: return (np.nan, np.nan)
    out = np.empty(nboot)
    for i in range(nboot):
        chosen = rng.choice(keys, size=len(keys), replace=True)
        out[i] = np.nanmedian(np.concatenate([values_by_gene[k] for k in chosen]))
    return tuple(np.nanpercentile(out, [2.5,97.5]))


def cis_contribution():
    df = pd.read_csv(ASE_RUN / "output" / "TN_cis_trans_classified.tsv", sep="\t")
    denom = df.B_log2fc.abs() + df.trans_log2fc.abs()
    df["cis_contribution"] = np.where(denom > 0, df.B_log2fc.abs()/denom, np.nan)
    df["A_bin"] = pd.cut(df.A_log2fc.abs(), [-np.inf,1,2,3,4,np.inf], labels=["0–1","1–2","2–3","3–4","4+"])
    rng = np.random.default_rng(20260728)
    out=[]
    for phase in PHASES:
        for abin in ["0–1","1–2","2–3","3–4","4+"]:
            s=df[(df.stage_group==phase)&(df.A_bin.astype(str)==abin)].dropna(subset=["cis_contribution"])
            by={g:x.cis_contribution.to_numpy() for g,x in s.groupby("orthogroup")}
            lo,hi=clustered_bootstrap(by,rng)
            out.append([phase,abin,float(s.cis_contribution.median()),lo,hi,len(by),len(s),"exploratory_parent_n1"])
    out=pd.DataFrame(out,columns=["stage_group","abs_A_bin","median","ci_low","ci_high","genes","gene_stage_rows","inference_level"])
    write_source(out,"04C_TN_cis_contribution")
    plot_cis_contribution(out)


def plot_cis_contribution(out: pd.DataFrame):
    """Compact reference-style rendering of the precomputed cis summary."""
    fig,ax=plt.subplots(figsize=(3.045,3.25)); x=np.arange(5)
    phase_colors=[DPM_RED,DPM_GOLD,DPM_TEAL,DPM_BLUE]
    phase_labels=["Early (0–65 d)","Mid (80–140 d)","Late (155–185 d)","Postharvest (12–72 h)"]
    markers=["o","s","^","D"]
    offsets=np.array([-.075,-.025,.025,.075])
    for phase,label,col,marker,offset in zip(PHASES,phase_labels,phase_colors,markers,offsets):
        s=out[out.stage_group.eq(phase)].set_index("abs_A_bin").loc[["0–1","1–2","2–3","3–4","4+"]]
        y=s["median"].to_numpy(); yerr=np.vstack([y-s.ci_low.to_numpy(),s.ci_high.to_numpy()-y])
        ax.errorbar(x+offset,y,yerr=yerr,marker=marker,ms=4.0,lw=1.15,capsize=2.0,
                    capthick=.85,color=col,mfc="white",mew=1.0,label=label,zorder=3)
    ax.axhline(.5,color="#6A6A6A",lw=.75,ls=(0,(4,3)),zorder=1)
    ax.set_ylim(.28,.515); ax.set_yticks([.30,.35,.40,.45,.50])
    ax.set_xlim(-.35,4.35); ax.set_xticks(x); ax.set_xticklabels(["0–1","1–2","2–3","3–4","≥4"])
    ax.set_xlabel(r"Parental divergence, $|A|$ (log$_2$ scale)")
    ax.set_ylabel(r"Cis contribution  $|B|/(|B|+|A-B|)$")
    clean(ax)
    ax.legend(frameon=False,ncol=2,loc="lower left",bbox_to_anchor=(-.01,1.015),
              fontsize=6.25,handlelength=1.25,columnspacing=.65,handletextpad=.30,borderaxespad=0)
    ax.set_title("TN cis contribution",fontweight="bold",pad=42,loc="left")
    ax.text(4.27,.502,"Cis = trans",ha="right",va="bottom",fontsize=6.2,color="#666")
    fig.text(.5,.016,"Median ± gene-clustered bootstrap 95% CI\nExploratory parental layer (n = 1 per stage)",
             ha="center",fontsize=5.85,color="#555",linespacing=1.15)
    panel_letter(fig,"c",x=.018,y=.985)
    fig.subplots_adjust(left=.18,right=.98,bottom=.27,top=.74)
    save(fig,"04C_TN_cis_contribution")


def bh(p):
    p=np.asarray(p,float); order=np.argsort(p); ranked=p[order]
    q=np.minimum.accumulate((ranked*len(p)/np.arange(1,len(p)+1))[::-1])[::-1]
    out=np.empty_like(q); out[order]=np.minimum(q,1); return out


def inheritance_regulatory():
    reg=pd.read_csv(ASE_RUN/"output"/"TN_cis_trans_classified.tsv",sep="\t",
                    usecols=["gene_africa","stage","regulatory_class"])
    het=pd.read_csv(HET_RUN/"output"/"modes_classified.tsv",sep="\t",
                    usecols=["gene_id_africa_hap2","stage","model3"])
    m=reg.merge(het,left_on=["gene_africa","stage"],right_on=["gene_id_africa_hap2","stage"],how="inner")
    rows=["PDO","DO","ODO"]
    obs=pd.crosstab(m.model3,m.regulatory_class).reindex(index=rows,columns=REG_ORDER,fill_value=0)
    _,_,_,expected=chi2_contingency(obs.to_numpy())
    residual=(obs.to_numpy()-expected)/np.sqrt(expected)
    ps=[]; records=[]
    grand=obs.to_numpy().sum()
    for i,r in enumerate(rows):
        for j,c in enumerate(REG_ORDER):
            a=int(obs.iloc[i,j]); b=int(obs.iloc[i,:].sum()-a); cc=int(obs.iloc[:,j].sum()-a); d=int(grand-a-b-cc)
            p=fisher_exact([[a,b],[cc,d]],alternative="two-sided").pvalue; ps.append(p)
            records.append([r,c,a,expected[i,j],100*a/obs.iloc[i,:].sum(),residual[i,j],p])
    q=bh(ps)
    out=pd.DataFrame(records,columns=["inheritance_class","regulatory_class","genes","expected","row_percentage","pearson_residual","fisher_p"])
    out["fisher_q_BH"]=q
    out["significance"]=np.select([q<.001,q<.01,q<.05],["***","**","*"],default="")
    out["inference_level"]="exploratory_parent_n1"
    write_source(out,"04D_TN_inheritance_regulatory_association")
    arr=out.pivot(index="inheritance_class",columns="regulatory_class",values="pearson_residual").reindex(index=rows,columns=REG_ORDER).to_numpy()
    vmax=max(2,float(np.nanpercentile(np.abs(arr),95)))
    fig,ax=plt.subplots(figsize=(4.70,3.05))
    im=ax.imshow(arr,cmap=DPM_HEATMAP_CMAP,norm=TwoSlopeNorm(vmin=-vmax,vcenter=0,vmax=vmax),aspect="auto")
    for i,r in enumerate(rows):
        for j,c in enumerate(REG_ORDER):
            z=out[(out.inheritance_class==r)&(out.regulatory_class==c)].iloc[0]
            ax.text(j,i,f"{z.row_percentage:.1f}%\n{z.significance}",ha="center",va="center",
                    fontsize=6.7,color="white" if abs(arr[i,j])>.55*vmax else "#222",fontweight="bold" if z.significance else "normal")
    short_labels=["I  Cis","II  Trans","III  Enh.","IV  Cis+trans comp.","V  Compensatory","VI  Conserved","VII  Ambiguous"]
    ax.set_xticks(range(7)); ax.set_xticklabels(short_labels,rotation=35,ha="right",rotation_mode="anchor",fontsize=6.15)
    ax.set_yticks(range(3)); ax.set_yticklabels(rows)
    for s in ax.spines.values(): s.set_visible(False)
    ax.set_title("TN inheritance × regulatory-class association",fontweight="bold",pad=10)
    cb=fig.colorbar(im,ax=ax,pad=.015,fraction=.025);cb.set_label("Pearson residual")
    fig.text(.98,.018,"Cell: within-inheritance percentage; Fisher exact test, BH-adjusted\n* q<0.05, ** q<0.01, *** q<0.001\nExploratory parent n = 1 per stage",
             ha="right",fontsize=5.7,color="#555",linespacing=1.12)
    panel_letter(fig,"D")
    fig.subplots_adjust(left=.13,right=.93,bottom=.35,top=.82)
    save(fig,"04D_TN_inheritance_regulatory_association")


def mode_heatmap(model: str, modes: list[str], stem: str, letter: str, title: str, figsize=(7.2,4.6), divider=None):
    df=pd.read_csv(HET_RUN/"output"/"stage_mode_summary.tsv",sep="\t")
    d=df[(df.model==model)&(df["mode"].isin(modes))].copy()
    write_source(d,stem)
    val=d.pivot(index="stage",columns="mode",values="percentage").reindex(index=STAGES,columns=modes).fillna(0)
    cnt=d.pivot(index="stage",columns="mode",values="genes").reindex(index=STAGES,columns=modes).fillna(0)
    arr=val.to_numpy(); vmin,vmax=float(arr.min()),float(arr.max()); center=float(np.nanmedian(arr))
    fig,ax=plt.subplots(figsize=figsize)
    im=ax.imshow(arr,cmap=DPM_HEATMAP_CMAP,norm=TwoSlopeNorm(vmin=vmin,vcenter=center,vmax=vmax),aspect="auto")
    for i in range(len(STAGES)):
        for j in range(len(modes)):
            ax.text(j,i,f"{int(cnt.iloc[i,j]):,}",ha="center",va="center",fontsize=5.1,
                    color="white" if abs(arr[i,j]-center)>.36*(vmax-vmin) else "#222")
    ax.set_xticks(range(len(modes)));ax.set_xticklabels(modes,fontweight="bold")
    ax.xaxis.tick_top();ax.tick_params(top=True,bottom=False,labeltop=True,labelbottom=False)
    ax.set_yticks(range(len(STAGES)));ax.set_yticklabels(STAGES)
    if divider is not None: ax.axvline(divider-.5,color="white",lw=2.6)
    ax.axhline(12.5,color="#333",lw=.8,ls="--")
    for s in ax.spines.values():s.set_linewidth(.7);s.set_color("#777")
    cb=fig.colorbar(im,ax=ax,pad=.014,fraction=.028);cb.set_label("Within-stage genes (%)")
    ax.set_title(title,fontweight="bold",pad=28)
    ax.text(1,-.075,"Numbers are genes; dashed line separates development and postharvest. Exploratory parent n = 1.",
            transform=ax.transAxes,ha="right",fontsize=6.3,color="#555")
    panel_letter(fig,letter)
    fig.subplots_adjust(left=.10,right=.94,bottom=.08,top=.84)
    save(fig,stem)


def trait_complement():
    df=pd.read_csv(ASE_RUN/"output"/"trait_module_ASE_summary.tsv",sep="\t")
    write_source(df,"06A_FL_TN_trait_complement")
    modules=["Oil biosynthesis & storage","De-novo / saturated FA","Unsaturated FA",
             "TAG assembly & oil body","Lipid oxidation / antioxidant","Shell / cell wall / lignin"]
    labels=["Oil biosynthesis\n& storage","De-novo /\nsaturated FA","Unsaturated FA",
            "TAG assembly\n& oil body","Lipid oxidation /\nantioxidant","Shell / cell wall\n/ lignin"]
    vmax=max(1,float(np.nanpercentile(np.abs(df.robust_median_log2_ratio),98)))
    norm=TwoSlopeNorm(vmin=-vmax,vcenter=0,vmax=vmax)
    fig,axs=plt.subplots(1,2,figsize=(7.4,4.1),sharey=True)
    for ax,analysis in zip(axs,["FL","TN"]):
        d=df[df.analysis.eq(analysis)]
        for yi,mod in enumerate(modules[::-1]):
            for xi,phase in enumerate(PHASES):
                z=d[(d.trait_module==mod)&(d.stage_group==phase)]
                if z.empty: continue
                z=z.iloc[0]; size=20+2.7*z.robust_ASE_percentage
                ax.scatter(xi,yi,s=size,c=[z.robust_median_log2_ratio],cmap="coolwarm",norm=norm,
                           edgecolor="white",lw=.55)
                ax.text(xi,yi,f"{z.robust_ASE_percentage:.0f}",ha="center",va="center",fontsize=5.9,
                        color="white" if abs(z.robust_median_log2_ratio)>.55*vmax else "#222")
        ax.set_xticks(range(4));ax.set_xticklabels(["0–65 d","80–140 d","155–185 d","12–72 h"],rotation=35,ha="right")
        ax.set_yticks(range(6));ax.set_yticklabels(labels[::-1])
        ax.set_xlim(-.6,3.6);ax.set_ylim(-.65,5.65);clean(ax)
        ax.set_title(analysis,fontweight="bold",pad=8);ax.grid(color="#E5E7E9",lw=.55,zorder=0)
        ax.set_xlabel("Stage group")
    sm=mpl.cm.ScalarMappable(norm=norm,cmap="coolwarm")
    cb=fig.colorbar(sm,ax=axs,pad=.02,fraction=.028);cb.set_label("Robust ASE median log₂(A/B)")
    fig.suptitle("Trait-linked haplotype-expression complementarity",fontweight="bold",y=.98)
    fig.text(.5,.025,"Bubble area and label: robust ASE (%). A/B: FL Africa hap2/American hap1; TN Dura(TK-like)/Pisifera(NS-like).\nExpression support is not a causal trait assignment.",
             ha="center",fontsize=6.5,color="#555")
    panel_letter(fig,"A")
    fig.subplots_adjust(left=.22,right=.91,bottom=.25,top=.84,wspace=.16)
    save(fig,"06A_FL_TN_trait_complement")


def write_index():
    rows=[
      ["03A","03A_FL_DBA_DEBA_chromosomes","FL chromosome-level DBA/DEBA density","DBA candidates and robust ASE-derived DEBA subset; Africa hap2 coordinates"],
      ["03B","03B_TN_DBA_DEBA_chromosomes","TN chromosome-level DBA/DEBA density","DBA candidates and robust ASE-derived DEBA subset; Dura/TK-like coordinates"],
      ["04A","04A_TN_cis_trans_decision_tree","TN cis/trans classification logic","Exploratory parent n=1; replicated F1 ASE"],
      ["04B","04B_TN_regulatory_classes_19stages","TN seven regulatory classes across 19 stages","Bubble area encodes within-stage percentage"],
      ["04C","04C_TN_cis_contribution","TN cis contribution","Median and gene-clustered bootstrap 95% CI"],
      ["04D","04D_TN_inheritance_regulatory_association","Inheritance by regulatory class","Pearson residual; Fisher exact BH q-values"],
      ["05A","05A_TN_12mode_matrix","TN Model 1, 12 modes","Exploratory parent n=1"],
      ["05B","05B_TN_5mode_matrix","TN Model 2, five modes","Five parental-expression modes"],
      ["05C","05C_TN_3mode_matrix","TN Model 3, three modes","Three inheritance modes"],
      ["06A","06A_FL_trait_candidate_ASE","FL trait-linked candidate ASE detail across 19 stages","Detailed candidate view; A=Africa hap2; B=American hap1; expression support only"],
      ["06B","06B_TN_trait_candidate_ASE","TN trait-linked candidate ASE detail across 19 stages","Detailed candidate view; A=Dura/TK-like; B=Pisifera/NS-like; expression support only"],
      ["06C","06C_FL_parent_independent_complementarity","FL single-panel functional haplotype complementarity","Family-level stage-group medians; ranges are descriptive, not confidence intervals"],
      ["06D","06A_FL_TN_trait_complement","Combined FL/TN trait-module line summary","Preferred main complementarity overview; expression support only; no causal claim"],
      ["07A","07A_FL_GO_GSEA","FL stage-group allele-direction GO-GSEA","Exploratory trend only; no GO term passed BH FDR < 0.05"],
      ["07B","07B_TN_GO_GSEA","TN stage-group allele-direction GO-GSEA","Exploratory trend only; no GO term passed BH FDR < 0.05"],
      ["07C","07C_FL_ASE_class_GO","FL ASE-class GO ORA","Exploratory trend only; no GO term passed BH FDR < 0.05"],
      ["07D","07D_TN_ASE_class_GO","TN ASE-class GO ORA","Exploratory trend only; no GO term passed BH FDR < 0.05"],
    ]
    pd.DataFrame(rows,columns=["panel","stem","content","interpretation_guardrail"]).to_csv(OUT/"NEW_PANEL_INDEX.tsv",sep="\t",index=False)


def main():
    style()
    chromosome_panel("FL","A","03A_FL_DBA_DEBA_chromosomes")
    chromosome_panel("TN","B","03B_TN_DBA_DEBA_chromosomes")
    decision_tree()
    regulatory_bubble()
    cis_contribution()
    inheritance_regulatory()
    mode_heatmap("model1",[f"M{i}" for i in range(1,13)],"05A_TN_12mode_matrix","A","TN 12 expression-inheritance modes",(7.4,4.7))
    mode_heatmap("model2",["H2P","L2P","B2P","CHP","CLP"],"05B_TN_5mode_matrix","B","TN five parental-expression modes",(4.8,4.7))
    mode_heatmap("model3",["ODO","PDO","DO"],"05C_TN_3mode_matrix","C","TN three inheritance modes",(3.9,4.7))
    trait_complement()
    write_index()


if __name__ == "__main__":
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        main()
