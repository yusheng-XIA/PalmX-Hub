#!/usr/bin/env python3
"""Extended Data Fig. 5 redraw from Source Data ED6a_dSV_counts / ED6b_dSV_positions."""
from pathlib import Path
import matplotlib as mpl, matplotlib.pyplot as plt, pandas as pd, numpy as np
HERE = Path(__file__).resolve().parent; MM = 1/25.4
mpl.rcParams.update({"font.family":"Arial","font.size":7,"axes.labelsize":7,"xtick.labelsize":6,"ytick.labelsize":6,
  "legend.fontsize":6,"axes.linewidth":0.6,"xtick.major.width":0.6,"ytick.major.width":0.6,"xtick.major.size":2.5,
  "ytick.major.size":2.5,"axes.unicode_minus":True,"pdf.fonttype":42})
a = pd.read_csv(HERE/"_raw_ED6a.tsv", sep="\t"); b = pd.read_csv(HERE/"_raw_ED6b.tsv", sep="\t")
assert a.DSV_Count.sum()==1480 and len(b)==1480
COL = {"INS":"#3C78A8","DEL":"#E69F00","DUP":"#CC79A7","INV":"#009E73"}
fig = plt.figure(figsize=(180*MM,150*MM))
ax = fig.add_axes([0.075,0.70,0.9,0.25])
x = np.arange(len(a)); ax.bar(x, a.DSV_Count, color="#2A9D8F", width=0.72)
for xi,v in zip(x,a.DSV_Count): ax.text(xi, v+3, str(v), ha="center", va="bottom", fontsize=5.5)
ax.set_xticks(x, [c.replace("chr","") for c in a.Chrom]); ax.set_xlim(-0.6,len(a)-0.4)
ax.set_ylabel("Candidate dSVs (n)"); ax.set_xlabel("Chromosome", labelpad=2)
for s in ("top","right"): ax.spines[s].set_visible(False)
ax2 = fig.add_axes([0.075,0.07,0.9,0.52])
chroms = list(a.Chrom); L = dict(zip(a.Chrom, a.Chrom_Length_bp/1e6))
for i,c in enumerate(chroms):
    y = len(chroms)-1-i
    ax2.add_patch(mpl.patches.FancyBboxPatch((0,y-0.18),L[c],0.36,boxstyle="round,pad=0,rounding_size=0.18",fc="#E3E5E8",ec="none",mutation_aspect=0.05))
    s = b[b.Chrom==c]
    for t in ("INS","DEL","DUP","INV"):
        p = s[s.SVTYPE==t].Position_Mb
        ax2.vlines(p, y-0.28, y+0.28, color=COL[t], lw=0.45, zorder=3 if t in("DUP","INV") else 2)
ax2.set_yticks(range(len(chroms)), chroms[::-1]); ax2.tick_params(axis="y", length=0)
ax2.set_xlim(0,180); ax2.set_ylim(-0.7,len(chroms)-0.3); ax2.set_xlabel("Position (Mb)")
for s in ("top","right","left"): ax2.spines[s].set_visible(False)
h=[mpl.lines.Line2D([],[],color=COL[t],lw=1.5,label=f"{t} ({(b.SVTYPE==t).sum():,})") for t in COL]
ax2.legend(handles=h, frameon=False, ncol=4, loc="lower right", bbox_to_anchor=(1.0,1.0), handlelength=1.0, columnspacing=1.2)
for axx,l in ((ax,"a"),(ax2,"b")):
    p=axx.get_position(); fig.text(0.008,p.y1+0.012,l,fontsize=9,fontweight="bold",va="bottom")
for ext in ("pdf","png"): fig.savefig(HERE/f"Extended_Data_Fig_05.{ext}", dpi=600, facecolor="white")
a[["Chrom","Chrom_Length_bp","DSV_Count"]].to_csv(HERE/"ED5a_dSV_counts.tsv",sep="\t",index=False)
b[["Chrom","Position_bp","SVTYPE","dSV_v5_Evidence_Class"]].to_csv(HERE/"ED5b_dSV_positions.tsv",sep="\t",index=False)
print(b.SVTYPE.value_counts().to_dict(), a.DSV_Count.max())
