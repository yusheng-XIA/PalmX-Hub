#!/usr/bin/env python3
"""Final source-locked redraw of Figure 3b, 3k, and 3l."""
from __future__ import annotations
import hashlib, os, tempfile
from pathlib import Path
import fitz
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from matplotlib.lines import Line2D
from matplotlib.patches import Patch, Rectangle
import numpy as np
import pandas as pd

OUT=Path(__file__).resolve().parent
MM=1/25.4
REF=Path('${ANALYSIS_DIR}/22_answer_reviews/00_ms/03_V3/03_figure3/04_FL_TN_ASE_reference_panels')
SRC_B=Path('${ANALYSIS_DIR}/21_MS/02_result/01_figure/plot_allele_classification.R')
SRC_K=REF/'source_06A_FL_TN_trait_complement.tsv'
SRC_BINS=Path('${ANALYSIS_DIR}/22_answer_reviews/00_ms/03_V3/03_figure3/11_diploid_32chrom_ancestry_breeding_20260805/results/diploid_ancestry_bins.tsv')
RUN=Path('${ANALYSIS_DIR}/22_answer_reviews/00_ms/03_V3/03_figure3/13_comprehensive_primitive_ancestry_breeding_20260805/runs/RUN-COMP-BREEDING-DPM-20260805-001')
SRC_BASE=RUN/'results/P9_final_attempt003_favorable_rescue/favorable_DPM_targets.tsv'
SRC_RESCUE=RUN/'results/P9_final_attempt003_favorable_rescue/favorable_haplotype_rescue_candidates.tsv'
SOURCES=[SRC_B,SRC_K,SRC_BINS,SRC_BASE,SRC_RESCUE]

CL={"NoDiff":"#F8766D","Sub":"#8DD3C7","HapDom":"#FFF19A","NoASE":"#58A6CF"}
PHASES=["Days 0–65","Days 80–140","Days 155–185","Hours 12–72"]
MODULES=["Oil biosynthesis & storage","TAG assembly & oil body","De-novo / saturated FA","Unsaturated FA","Lipid oxidation / antioxidant","Shell / cell wall / lignin"]
ABBR=["OBS","TOF","DSF","UFA","LOD","SCL"]


def style():
    mpl.rcParams.update({
        'font.family':'sans-serif','font.sans-serif':['Arial','Liberation Sans','DejaVu Sans'],
        'font.size':6,'axes.labelsize':6,'xtick.labelsize':6,'ytick.labelsize':6,
        'legend.fontsize':6,'legend.title_fontsize':6,'axes.linewidth':.65,
        'xtick.major.width':.65,'ytick.major.width':.65,'pdf.fonttype':42,'ps.fonttype':42,
        'savefig.facecolor':'white','figure.facecolor':'white'
    })


def letter(fig,s): fig.text(.012,.985,s,ha='left',va='top',fontsize=8,fontweight='bold')
def save(fig,p): fig.savefig(p,facecolor='white',metadata={'Title':p.stem,'Creator':Path(__file__).name}); plt.close(fig)
def raster(pdf,png):
    d=fitz.open(pdf); pix=d[0].get_pixmap(matrix=fitz.Matrix(600/72,600/72),alpha=False); pix.save(png); d.close()
def digest(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1<<20),b''): h.update(b)
    return h.hexdigest()


def allele_statistics():
    # Exact source-script counts in category order:
    # Biallelic, allele with same CDS, haplotype-specific.
    expected='27106,264,3214, 27106,264,7430, 18439,7747,7720, 18439,7747,8486'
    compact=''.join(SRC_B.read_text().split())
    assert expected.replace(' ','') in compact, 'Reviewed allele counts not found in source R script'
    labels=['FL_hap1','FL_hap2','TN_hap1','TN_hap2']
    counts=np.array([[27106,264,3214],[27106,264,7430],[18439,7747,7720],[18439,7747,8486]],dtype=int)
    pct=counts/counts.sum(axis=1,keepdims=True)*100.0
    return labels,counts,pct


def render_b(path):
    labels,counts,pct=allele_statistics()
    colors=['#76B2D4','#D58774','#F1D992']
    fig=plt.figure(figsize=(58*MM,30*MM)); letter(fig,'b')
    ax=fig.add_axes([.17,.24,.80,.67])
    x=np.arange(4); bottom=np.zeros(4)
    for j,color in enumerate(colors):
        vals=pct[:,j]
        ax.bar(x,vals,bottom=bottom,width=.72,color=color,edgecolor='white',linewidth=.45)
        # Label only segments large enough to contain 6-pt text; the two
        # 0.8% same-CDS slivers remain visible but are not overprinted.
        for i,v in enumerate(vals):
            if v>=8:
                ax.text(i,bottom[i]+v/2,f'{v:.1f}',ha='center',va='center',fontsize=6)
        bottom+=vals
    ax.axvline(1.5,color='#777777',lw=.65,ls=(0,(1.5,2.4)))
    ax.set_xlim(-.55,3.55); ax.set_ylim(0,102)
    ax.set_yticks([0,25,50,75,100]); ax.set_ylabel('Proportion of genes (%)',labelpad=2)
    ax.set_xticks(x,labels); ax.tick_params(axis='x',length=0,pad=2); ax.tick_params(axis='y',pad=1.5)
    ax.spines[['top','right']].set_visible(False)
    save(fig,path)


def combine_ab(a,b,out):
    da,db=fitz.open(a),fitz.open(b); w=58*72/25.4; ha=22*72/25.4; hb=30*72/25.4
    doc=fitz.open(); p=doc.new_page(width=w,height=ha+hb)
    p.show_pdf_page(fitz.Rect(0,0,w,ha),da,0); p.show_pdf_page(fitz.Rect(0,ha,w,ha+hb),db,0)
    doc.set_metadata({'title':'Fig3ab','subject':'Source-locked final panels'}); doc.save(out,garbage=4,deflate=True)
    doc.close(); da.close(); db.close()


def render_k(path):
    d=pd.read_csv(SRC_K,sep='\t')
    fig=plt.figure(figsize=(69*MM,40*MM)); letter(fig,'k')
    ax=fig.add_axes([.12,.34,.84,.48])
    norm=TwoSlopeNorm(vmin=-2.5,vcenter=0,vmax=2.5)
    for yi,mod in enumerate(MODULES):
        for pi,phase in enumerate(PHASES):
            for ai,analysis in enumerate(['FL','TN']):
                q=d[(d.trait_module==mod)&(d.stage_group==phase)&(d.analysis==analysis)]
                if q.empty: continue
                r=q.iloc[0]; x=pi*2+ai
                ax.scatter(x,yi,s=7+.21*float(r.robust_ASE_percentage),
                           c=[float(r.robust_median_log2_ratio)],cmap='coolwarm',norm=norm,
                           marker='o' if analysis=='FL' else 's',
                           edgecolor='#E75F5F' if analysis=='FL' else '#2F68A2',linewidth=.55)
    ax.set_xlim(-.6,7.6); ax.set_ylim(5.6,-.7); ax.set_yticks(range(6),ABBR)
    ax.set_xticks(range(8),['FL','TN']*4); ax.xaxis.tick_top(); ax.tick_params(top=True,bottom=False,labeltop=True,labelbottom=False,length=0,pad=1.5)
    for pi in range(4): ax.axvspan(pi*2-.5,pi*2+1.5,color='#F5F7F7' if pi<3 else '#FFF7E5',zorder=-2)
    for s in ax.spines.values(): s.set_visible(False)
    for x,name in zip([.225,.435,.645,.855],['Early','Middle','Late','Postharvest']):
        fig.text(x,.965,name,ha='center',va='top',fontweight='bold')

    # Dedicated, non-overlapping lower-left size legend.
    lax=fig.add_axes([.12,.035,.35,.19]); lax.axis('off')
    handles=[Line2D([0],[0],marker='o',ls='',mfc='#DDE2E5',mec='#777777',mew=.45,
                    markersize=np.sqrt(7+.21*v),label=f'{v}') for v in [50,70,90]]
    lax.legend(handles=handles,title='Robust ASE (%)',frameon=False,ncol=3,loc='upper left',
               bbox_to_anchor=(0,1),borderaxespad=0,columnspacing=.40,handletextpad=.18)
    # Dedicated lower-right colorbar, separated from size legend and labels.
    cax=fig.add_axes([.61,.15,.33,.045])
    cb=fig.colorbar(mpl.cm.ScalarMappable(norm=norm,cmap='coolwarm'),cax=cax,orientation='horizontal')
    cb.set_ticks([-2.5,0,2.5]); cb.ax.tick_params(labelsize=6,pad=1,length=2)
    fig.text(.775,.035,'Median log2(A/B)',ha='center',va='bottom')
    save(fig,path)


def render_l(path):
    bins=pd.read_csv(SRC_BINS,sep='\t'); bins=bins[bins.Individual.eq('FL')].copy()
    base=pd.read_csv(SRC_BASE,sep='\t',low_memory=False); rescue=pd.read_csv(SRC_RESCUE,sep='\t',low_memory=False)
    assert len(base)==293 and len(rescue)==43
    fig=plt.figure(figsize=(183*MM,72*MM)); letter(fig,'l')
    fig.text(.5,.975,'Physical position (Mb)',ha='center',va='top',fontweight='bold')
    axes=[fig.add_axes([.055,.315,.425,.57]),fig.add_axes([.545,.315,.425,.57])]
    anc={'Dura':'#3767A6','Pisifera':'#E56565','Meizhou4':'#D8B64C'}
    geno={'D':'#3767A6','P':'#E56565','M':'#D8B64C'}
    action={'INTRODUCE_OR_TUNE':('^','#159D82'),'RETAIN_FL':('^','#4F70B5'),'TIMING_SCREEN':('>','#E3A018')}
    actual={'TN_h1':'#00897B','TN_h2':'#6BC5BA','FL_HapA':'#7651A8','FL_HapB':'#B091CF'}

    for col,ax in enumerate(axes):
        chroms=[f'chr{i:02d}' for i in range(1+8*col,9+8*col)]
        for yi,chrom in enumerate(chroms):
            y=7-yi
            if yi%2==0: ax.axhspan(y-.42,y+.42,color='#F6F7F7',zorder=-5)
            q=bins[bins.Chromosome.eq(chrom)].sort_values('Bin_index'); lengths=[]
            for hap,dy in [(1,.135),(2,-.135)]:
                ec=f'Hap{hap}_end0'; sc=f'Hap{hap}_start0'; ac=f'Hap{hap}_ancestry'
                L=float(q[ec].max())/1e6; lengths.append(L)
                for r in q.itertuples():
                    x0=float(getattr(r,sc))/1e6; x1=float(getattr(r,ec))/1e6; a=str(getattr(r,ac)); known=a in anc
                    ax.add_patch(Rectangle((x0,y+dy-.055),x1-x0,.11,fc=anc.get(a,'white'),
                                           ec='none' if known else '#8A8A8A',lw=.2,
                                           hatch=None if known else '////',zorder=1))
                ax.add_patch(Rectangle((0,y+dy-.055),L,.11,fill=False,ec='#333333',lw=.35,zorder=2))
            # Full, source-locked target set.
            for r in base[base.Chromosome.eq(chrom)].itertuples():
                frac=float(r.chromosome_fraction); alleles=str(r.target_DPM_genotype).split('/')
                for hi,(dy,L) in enumerate([(.135,lengths[0]),(-.135,lengths[1])]):
                    if hi<len(alleles) and alleles[hi] in geno:
                        ax.scatter(frac*L,y+dy,s=4.0,marker='s',c=geno[alleles[hi]],edgecolor='#222',lw=.18,zorder=5)
                mk,co=action.get(str(r.breeding_action),('^','#888'))
                ax.scatter(frac*max(lengths),y+.29,s=6.0,marker=mk,c=co,edgecolor='none',zorder=6)
            for r in rescue[rescue.Chromosome.eq(chrom)].itertuples():
                frac=float(r.chromosome_fraction)
                ax.scatter(frac*max(lengths),y+.29,s=5.5,marker='D',c='#8A5FB2',edgecolor='none',zorder=6)
                for hn in [str(getattr(r,'TN_ASE_hap','')),str(getattr(r,'FL_ASE_hap',''))]:
                    if hn in actual:
                        hi=0 if hn in {'TN_h1','FL_HapA'} else 1
                        ax.scatter(frac*lengths[hi],y+(.135 if hi==0 else -.135),s=4.2,marker='s',c=actual[hn],edgecolor='#222',lw=.18,zorder=7)
            nice=str(int(chrom[-2:]))
            ax.text(-16.8,y,f'chr{nice}',ha='right',va='center',fontweight='bold')
            ax.text(-3.1,y+.135,'H1',ha='right',va='center')
            ax.text(-3.1,y-.135,'H2',ha='right',va='center')
        ax.set_xlim(-20,202); ax.set_ylim(-.55,7.55); ax.set_yticks([])
        ax.xaxis.tick_top(); ax.set_xticks(np.arange(0,201,25)); ax.tick_params(axis='x',pad=1,length=2.5,direction='out')
        ax.grid(axis='x',color='#DDDDDD',lw=.35,zorder=-6)
        for s in ['left','right','bottom']: ax.spines[s].set_visible(False)
        ax.spines['top'].set_color('#777777'); ax.spines['top'].set_linewidth(.6)

    # Four separated legend blocks in the reserved lower band.
    groups=[]
    groups.append((.035,.035,.23,.22,[Patch(fc=anc[k],ec='#444',lw=.3,label=k) for k in ['Dura','Pisifera','Meizhou4']]+[Patch(fc='white',ec='#888',hatch='////',label='Pending')],'FL ancestry background',2))
    groups.append((.29,.035,.25,.22,[Line2D([0],[0],marker=action[k][0],ls='',mfc=action[k][1],mec='none',markersize=4.5,label=l) for k,l in [('INTRODUCE_OR_TUNE','Introduce/tune'),('RETAIN_FL','Retain FL'),('TIMING_SCREEN','Timing screen')]]+[Line2D([0],[0],marker='D',ls='',mfc='#8A5FB2',mec='none',markersize=4,label='Haplotype screen')],'Favorable action (293 targets)',2))
    groups.append((.57,.035,.18,.22,[Line2D([0],[0],marker='s',ls='',mfc=geno[k],mec='#222',mew=.2,markersize=4,label=f'{k} = {v}') for k,v in [('D','Dura'),('P','Pisifera'),('M','Meizhou4')]],'Diploid target chips',1))
    groups.append((.78,.035,.19,.22,[Line2D([0],[0],marker='s',ls='',mfc=v,mec='#222',mew=.2,markersize=4,label=k) for k,v in actual.items()],'Actual-haplotype chips (n = 43)',2))
    for x,y,w,h,handles,title,ncol in groups:
        a=fig.add_axes([x,y,w,h]); a.axis('off'); a.legend(handles=handles,title=title,frameon=False,ncol=ncol,loc='upper left',borderaxespad=0,columnspacing=.65,handlelength=1.0,handletextpad=.35,labelspacing=.35)
    save(fig,path)


def main():
    style()
    for p in SOURCES:
        if not p.is_file(): raise FileNotFoundError(p)
    if not (OUT/'Fig3a.pdf').is_file(): raise FileNotFoundError(OUT/'Fig3a.pdf')
    with tempfile.TemporaryDirectory(prefix='.final_bkl_',dir=OUT) as td:
        td=Path(td)
        render_b(td/'Fig3b.pdf'); combine_ab(OUT/'Fig3a.pdf',td/'Fig3b.pdf',td/'Fig3ab.pdf')
        render_k(td/'Fig3k.pdf'); render_l(td/'Fig3l.pdf')
        for s in ['Fig3b','Fig3ab','Fig3k','Fig3l']: raster(td/f'{s}.pdf',td/f'{s}_600dpi.png')
        for p in td.glob('Fig3*'): os.replace(p,OUT/p.name)
    outs=[OUT/f'{s}{e}' for s in ['Fig3b','Fig3ab','Fig3k','Fig3l'] for e in ['.pdf','_600dpi.png']]
    (OUT/'FIGURE3_BKL_FINAL_CHECKSUMS.sha256').write_text('\n'.join(f'{digest(p)}  {p.name}' for p in outs)+'\n')
    (OUT/'FIGURE3_BKL_FINAL_SOURCES.tsv').write_text('sha256\tpath\n'+'\n'.join(f'{digest(p)}\t{p}' for p in SOURCES)+'\n')
    labels,counts,pct=allele_statistics()
    rows=['haplotype\tcategory\tcount\tpercentage']
    cats=['Biallelic','Allele_with_same_CDS','Haplotype_specific']
    for i,hap in enumerate(labels):
        for j,cat in enumerate(cats): rows.append(f'{hap}\t{cat}\t{counts[i,j]}\t{pct[i,j]:.6f}')
    (OUT/'Fig3b_allele_statistics_source_data.tsv').write_text('\n'.join(rows)+'\n')

if __name__=='__main__': main()
