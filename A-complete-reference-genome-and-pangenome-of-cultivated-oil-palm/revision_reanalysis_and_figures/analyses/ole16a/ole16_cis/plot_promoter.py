import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Polygon
plt.rcParams.update({'font.family':'Arial','font.size':6,'pdf.fonttype':42})
# classes (own coordinates relative to ATG; kb)
classes=[
 ('E. oleifera-type','MZ4-Hap1, MZ4-Hap2, FL-Hap1',3,-23.0,[(-18.47,-7.39,'eo')],'FL-Hap1: 99.9% of FL OLE16a reads'),
 ('E. guineensis A','20 haplotypes incl. FL-Hap2, TK-Hap1/2, TN-Hap1',20,-11.95,[],'FL-Hap2, TN-Hap1: silent'),
 ('E. guineensis B','12 haplotypes incl. TN-Hap2, NS-Hap1',12,-21.4,[(-11.35,-1.91,'b')],'TN-Hap2: silent'),
 ('E. guineensis C','NS-Hap2, EG_071, EG_072',3,-34.03,[(-26.68,-15.60,'rel'),(-12.91,-1.91,'rel')],''),
 ('E. guineensis D','EG_025',1,-14.42,[(-5.82,-3.29,'d')],''),
]
col={'eo':'#1b7f5f','b':'#9aa3ad','rel':'#6fb59b','d':'#c7ccd1'}
fig,ax=plt.subplots(figsize=(150/25.4,62/25.4))
y=0
for name,members,n,L,ins,expr in classes:
    ax.plot([L,0],[y,y],color='#444',lw=0.8,zorder=1)
    ax.add_patch(Rectangle((L-2.2,y-0.18),2.2,0.36,fc='#d9d9d9',ec='#666',lw=0.4))  # neighbour gene
    for a,b,c in ins:
        ax.add_patch(Rectangle((a,y-0.22),b-a,0.44,fc=col[c],ec='none',zorder=2))
        if c in ('eo','rel'):
            ax.add_patch(Rectangle((a,y-0.22),1.9,0.44,fc='none',ec='k',lw=0.3,hatch='////',zorder=3))
            ax.add_patch(Rectangle((b-1.9,y-0.22),1.9,0.44,fc='none',ec='k',lw=0.3,hatch='////',zorder=3))
    ax.add_patch(Polygon([[0,y-0.25],[0,y+0.25],[1.2,y]],fc='#c0392b',ec='none'))  # OLE16a
    for p in (-0.97,-0.16): ax.plot([p,p],[y+0.2,y+0.34],color='#b8860b',lw=0.6)
    ax.text(-60,y+0.05,f'{name} (n = {n})',ha='left',va='bottom',fontsize=6,fontweight='bold')
    ax.text(-60,y-0.05,members,ha='left',va='top',fontsize=5,color='#444')
    if expr: ax.text(1.8,y,expr,ha='left',va='center',fontsize=5)
    y-=1
ax.set_xlim(-60,16); ax.set_ylim(y+0.3,0.7)
ax.set_yticks([]); ax.set_xticks([-30,-20,-10,0]); ax.spines['bottom'].set_bounds(-36,0); ax.set_xticklabels(['-30','-20','-10','ATG'])
ax.set_xlabel('Distance upstream of OLE16a start codon (kb)')
for s in ['left','right','top']: ax.spines[s].set_visible(False)
from matplotlib.lines import Line2D
h=[Rectangle((0,0),1,1,fc=col['eo']),Rectangle((0,0),1,1,fc=col['rel']),Rectangle((0,0),1,1,fc=col['b']),Rectangle((0,0),1,1,fc='white',ec='k',hatch='////'),Line2D([0],[0],color='#b8860b',lw=0.8),Rectangle((0,0),1,1,fc='#d9d9d9',ec='#666')]
ax.legend(h,['E. oleifera-type LTR-RT insertion (11.1 kb)','Related LTR-RT insertions (91% identity)','Other insertions','LTR (~2 kb)','RY elements (conserved)','Upstream neighbour gene'],loc='upper center',fontsize=5,frameon=False,ncol=2,bbox_to_anchor=(0.5,-0.55))
plt.tight_layout()
plt.savefig('fig_promoter_draft.pdf',bbox_inches='tight'); plt.savefig('fig_promoter_draft.png',dpi=300,bbox_inches='tight')
