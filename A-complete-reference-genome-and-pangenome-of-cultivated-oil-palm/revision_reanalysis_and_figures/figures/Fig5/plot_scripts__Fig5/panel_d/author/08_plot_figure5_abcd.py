#!/usr/bin/env python3
from __future__ import annotations

import csv
import concurrent.futures as cf
import math
import re
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Wedge
import numpy as np
import pandas as pd
from PIL import Image, ImageDraw, ImageFont
from scipy.interpolate import PchipInterpolator
from scipy.optimize import curve_fit

WORK = Path('${ANALYSIS_DIR}/22_answer_reviews/00_ms/03_V3/05_figure/08_SV_hap39_Figure5abcd_20260808')
RES = WORK / 'results'
TABLE = WORK / 'figure5/tables'
PANELS = WORK / 'figure5/panels'
FINAL = WORK / 'figure5/final'
for d in (TABLE, PANELS, FINAL): d.mkdir(parents=True, exist_ok=True)

VARIANT_COUNTS = Path('${ANALYSIS_DIR}/21_MS/06_result/variant_counts.tsv')
SWAVE_PLOTDATA = Path('${ANALYSIS_DIR}/22_answer_reviews/00_ms/03_V3/05_figure/03_SWave_complexSV_maintext/tables/SWave_complexSV_bar_inset_nested_donut_20260626.plotdata.tsv')

SV_ORDER = ['DEL', 'INS', 'INV', 'DUP', 'TRA']
PLOT_ORDER = ['INS', 'DEL', 'INV', 'TRA', 'DUP']
COLORS = {'DEL':'#D94B4B','INS':'#43BBD2','INV':'#009B77','DUP':'#8063B0','TRA':'#E39D21'}
FREQ_LEVELS = ['Private','Low','Intermediate','High','Core']
FREQ_COLORS = {'Private':'#F4C9D4','Low':'#D9B3D9','Intermediate':'#B6A8D6','High':'#8FA8CE','Core':'#5B7DB1'}
SAMPLES = [
    'nrly_hap1','nrly_hap2','EG_037','EG_015','EG_041','American_hap1','meizhou4_hap1','meizhou4_hap2',
    'EG_008','EG_035','EG_033','EG_102','EG_075','EG_095','EG_062','EG_houke','EG_083','EG_107','EG_113',
    'bk_hap2','EG_086','bk_hap1','dura_hap1','dura_hap2','pisifera_hap1','pisifera_hap2','EG_072',
    'EG_071','EG_067','EG_057','EG_058','EG_090','EG_025','EG_017','EG_183','EG_065','EG_176','EG_146'
]
DISPLAY = {
    'nrly_hap1':'Nrly-hap1','nrly_hap2':'Nrly-hap2','American_hap1':'FL-hap1',
    'meizhou4_hap1':'MZ4-hap1','meizhou4_hap2':'MZ4-hap2','EG_houke':'Houke',
    'bk_hap1':'BK-hap1','bk_hap2':'BK-hap2','dura_hap1':'Dura-hap1','dura_hap2':'Dura-hap2',
    'pisifera_hap1':'Pisi-hap1','pisifera_hap2':'Pisi-hap2'
}

mpl.rcParams.update({'pdf.fonttype':42,'ps.fonttype':42,'svg.fonttype':'none','font.family':'DejaVu Sans','font.size':7,'axes.linewidth':0.65})

def save(fig, stem, dpi=600):
    for ext in ('pdf','svg'): fig.savefig(PANELS / f'{stem}.{ext}', bbox_inches='tight', facecolor='white')
    fig.savefig(PANELS / f'{stem}.png', bbox_inches='tight', facecolor='white', dpi=dpi)
    plt.close(fig)

def load_catalog(path):
    d = pd.read_csv(path, sep='\t')
    d['SVLEN_abs_bp'] = pd.to_numeric(d['SVLEN_Median_bp'], errors='coerce').abs()
    return d

def load_inputs():
    tier = load_catalog(RES / 'population_repeat/sv_population_catalog.tsv')
    syri = load_catalog(RES / 'syri_large_rearrangements_population/sv_population_catalog.tsv')
    tier_mem = pd.read_csv(RES / 'population_repeat/sv_cluster_membership.tsv', sep='\t')
    syri_mem = pd.read_csv(RES / 'syri_large_rearrangements_population/sv_cluster_membership.tsv', sep='\t')
    tier_mem['SVLEN_abs_bp'] = pd.to_numeric(tier_mem['SVLEN_bp'], errors='coerce').abs()
    syri_mem['SVLEN_abs_bp'] = pd.to_numeric(syri_mem['SVLEN_bp'], errors='coerce').abs()
    mixed = pd.concat([tier[tier.SVTYPE.isin(['DEL','INS'])], syri[syri.SVTYPE.isin(['INV','DUP','TRA'])]], ignore_index=True)
    mixed_mem = pd.concat([tier_mem[tier_mem.SVTYPE.isin(['DEL','INS'])], syri_mem[syri_mem.SVTYPE.isin(['INV','DUP','TRA'])]], ignore_index=True)
    return tier, mixed, mixed_mem

def new_syri_path(sample):
    if sample.startswith('meizhou4_'):
        return RES / f'new_callers/{sample}/syri/syri.vcf'
    row = pd.read_csv(WORK / 'manifests/new8.tsv', sep='\t').set_index('Sample').loc[sample]
    return Path(row['Existing_SyRI_VCF'])

def count_snp_indel(path):
    snp = indel = 0
    with open(path, errors='replace') as h:
        for line in h:
            if line.startswith('#'): continue
            f = line.rstrip('\n').split('\t')
            if len(f) < 8: continue
            rid = f[2].upper()
            if rid.startswith('SNP'):
                snp += 1
            elif rid.startswith('INS') or rid.startswith('DEL'):
                delta = abs(len(f[4].split(',')[0]) - len(f[3]))
                if 0 < delta < 50: indel += 1
    return snp, indel

def variant_table():
    cached_path = TABLE / 'Figure5a_variant_counts_hap38.tsv'
    if cached_path.exists():
        cached = pd.read_csv(cached_path, sep='\t')
        if set(cached.Sample) == set(SAMPLES) and len(cached) == len(SAMPLES):
            return cached.set_index('Sample')
    old = pd.read_csv(VARIANT_COUNTS, sep='\t').set_index('Sample')
    rows = []
    old_alias = {'EG_houke':'EG_houke'}
    new = {'dura_hap1','dura_hap2','pisifera_hap1','pisifera_hap2','nrly_hap1','nrly_hap2','meizhou4_hap1','meizhou4_hap2'}
    new_order = [s for s in SAMPLES if s in new]
    new_paths = [new_syri_path(s) for s in new_order]
    with cf.ProcessPoolExecutor(max_workers=4) as pool:
        new_values = list(pool.map(count_snp_indel, new_paths))
    new_counts = dict(zip(new_order, new_values))
    for s in SAMPLES:
        if s in new:
            snp, indel = new_counts[s]
        else:
            key = old_alias.get(s, s)
            snp, indel = int(old.loc[key,'SNP']), int(old.loc[key,'InDel'])
        rows.append({'Sample':s,'SNP':snp,'InDel':indel,'Source':'SyRI_direct' if s in new else 'published_variant_counts'})
    out = pd.DataFrame(rows)
    out.to_csv(TABLE / 'Figure5a_variant_counts_hap38.tsv', sep='\t', index=False)
    return out.set_index('Sample')

def draw_pie(ax, x, y, vals, radius, balanced=True):
    vals = np.asarray(vals, float)
    if balanced:
        vals = np.sqrt(np.clip(vals, 0, None))
        pos = vals > 0
        if pos.any(): vals[pos] = np.maximum(vals[pos], vals.sum()*0.025)
    if vals.sum() <= 0: return
    start = 90
    for value, typ in zip(vals, SV_ORDER):
        if value <= 0: continue
        theta = 360*value/vals.sum()
        ax.add_patch(Wedge((x,y), radius, start, start+theta, facecolor=COLORS[typ], edgecolor='white', linewidth=.25))
        start += theta
    ax.add_patch(Circle((x,y), radius, fill=False, edgecolor='#222222', linewidth=.4))

def plot_a(mixed, mem):
    variants = variant_table()
    agg = mem.groupby(['Sample','SVTYPE']).agg(SV_Count=('Cluster_ID','count'),SV_Length_bp=('SVLEN_abs_bp','sum')).reset_index()
    full = pd.MultiIndex.from_product([SAMPLES,SV_ORDER],names=['Sample','SVTYPE'])
    agg = agg.set_index(['Sample','SVTYPE']).reindex(full,fill_value=0).reset_index()
    summary = agg.groupby('Sample',as_index=False).agg(Total_Count=('SV_Count','sum'),Total_Length=('SV_Length_bp','sum')).set_index('Sample').loc[SAMPLES]
    agg.to_csv(TABLE/'Figure5a_plotdata_by_sample_svtype.tsv',sep='\t',index=False)
    summary.reset_index().to_csv(TABLE/'Figure5a_sample_summary.tsv',sep='\t',index=False)
    xs=np.arange(len(SAMPLES))*.75; base=.30; catx=xs[-1]+1.15; legx=catx+.90
    fig,ax=plt.subplots(figsize=(21.4,4.6)); ax.set_xlim(-1.10,legx+1.15); ax.set_ylim(-1.20,3.16); ax.set_aspect('equal'); ax.axis('off')
    ylen,ycount,yname,ysv,ysnp,yindel=2.35,1.15,.42,-.42,-.72,-1.00
    maxlen=summary.Total_Length.max(); maxcount=summary.Total_Count.max()
    for i,s in enumerate(SAMPLES):
        sub=agg[agg.Sample==s].set_index('SVTYPE').reindex(SV_ORDER,fill_value=0)
        lv=sub.SV_Length_bp.to_numpy(float); cv=sub.SV_Count.to_numpy(float)
        draw_pie(ax,xs[i],ylen,lv,base*math.sqrt(lv.sum()/maxlen)); draw_pie(ax,xs[i],ycount,cv,base*math.sqrt(cv.sum()/maxcount))
        ax.text(xs[i],yname,DISPLAY.get(s,s),ha='center',va='top',rotation=90,fontsize=4.85)
        ax.text(xs[i],ysv,f'{cv.sum()/1000:.1f}k',ha='center',va='center',fontsize=4.1,color='#555')
        ax.text(xs[i],ysnp,f'{variants.loc[s,"SNP"]/1e6:.1f}M',ha='center',va='center',fontsize=4.0,color='#555')
        ax.text(xs[i],yindel,f'{variants.loc[s,"InDel"]/1000:.0f}K',ha='center',va='center',fontsize=4.0,color='#666')
    ax.text(-.92,ylen,'Length',rotation=90,ha='center',va='center',fontsize=7.5); ax.text(-.92,ycount,'Number',rotation=90,ha='center',va='center',fontsize=7.5)
    ax.text(-.38,ysv,'SV records',ha='right',va='center',fontsize=4.8,color='#555'); ax.text(-.38,ysnp,'SNP',ha='right',va='center',fontsize=4.8,color='#555'); ax.text(-.38,yindel,'InDel',ha='right',va='center',fontsize=4.8,color='#555')
    for j,t in enumerate(SV_ORDER):
        x=xs[15]+(j-2)*1.25; ax.add_patch(Circle((x,3.00),.05,facecolor=COLORS[t],edgecolor='none')); ax.text(x+.08,3.00,t,va='center',fontsize=6.2)
    cc=mixed.groupby('SVTYPE').size().reindex(SV_ORDER,fill_value=0).to_numpy(float); cl=mixed.groupby('SVTYPE').SVLEN_abs_bp.sum().reindex(SV_ORDER,fill_value=0).to_numpy(float)
    draw_pie(ax,catx,ylen,cl,.37); draw_pie(ax,catx,ycount,cc,.37)
    ax.text(catx,yname,'Final\ncatalog',ha='center',va='top',fontsize=5.4); ax.text(catx,ysv,f'{cc.sum()/1000:.1f}k loci',ha='center',va='center',fontsize=4.8,fontweight='bold'); ax.text(catx,ysnp,'SV loci',ha='center',va='center',fontsize=4.2,color='#555')
    ax.plot([catx-.62,catx-.62],[-1.05,2.92],color='#ddd',lw=.5,ls=(0,(2,2)))
    for yy,maximum,refs,title,fmt in [(ylen,maxlen,[100e6,200e6,300e6],'Sample length',lambda v:f'{v/1e6:.0f} Mb'),(ycount,maxcount,[10000,20000,35000],'Sample count',lambda v:f'{v/1000:.0f}k')]:
        bottom=yy-base
        for ref in sorted(refs,reverse=True):
            if ref>maximum*1.1: continue
            r=base*math.sqrt(ref/maximum); ax.add_patch(Circle((legx,bottom+r),r,facecolor='white',edgecolor='#aaa',linewidth=.5)); ax.text(legx+base+.32,bottom+2*r,fmt(ref),va='center',fontsize=5.2)
        ax.text(legx,yy-base-.28,title,ha='center',va='top',fontsize=5.8)
    save(fig,'Figure5a_hap38_finalSV')

def freq_level(x): return str(x).split(' ')[0]

def plot_b(tier):
    d=tier[tier.SVTYPE.isin(SV_ORDER)].copy(); d['Level']=d.Frequency_Bin.map(freq_level)
    tab=d.groupby(['SVTYPE','Level']).size().unstack(fill_value=0).reindex(index=PLOT_ORDER,columns=FREQ_LEVELS,fill_value=0)
    tab.to_csv(TABLE/'Figure5b_frequency_classes_hap38.tsv',sep='\t')
    totals=tab.sum(axis=1)
    fig=plt.figure(figsize=(6.45,4.45)); ax=fig.add_axes([.11,.12,.72,.82]); bottom=np.zeros(len(PLOT_ORDER)); x=np.arange(len(PLOT_ORDER))
    for level in FREQ_LEVELS:
        vals=tab[level].to_numpy(); ax.bar(x,vals,bottom=bottom,width=.68,color=FREQ_COLORS[level],edgecolor='white',linewidth=.3); bottom+=vals
    for i,n in enumerate(totals): ax.text(i,n+max(totals)*.018,f'{int(n):,}',ha='center',va='bottom',fontsize=7.2)
    ax.set_xticks(x,labels=PLOT_ORDER); ax.set_ylabel('Number of SVs'); ax.set_xlim(-.62,4.62); ax.set_ylim(0,max(totals)*1.12); ax.spines[['top','right']].set_visible(False); ax.yaxis.set_major_formatter(lambda v,p:f'{int(v):,}')
    all_counts=d.Level.value_counts().reindex(FREQ_LEVELS,fill_value=0); all_total=all_counts.sum()
    ax.text(1.98,max(totals)*.99,'Frequency class',fontsize=8.0)
    for j,level in enumerate(FREQ_LEVELS):
        y=max(totals)*(.925-j*.055); ax.add_patch(plt.Rectangle((1.98,y-.018*max(totals)),.13,.036*max(totals),color=FREQ_COLORS[level],clip_on=False)); ax.text(2.15,y,f'{level} ({int(all_counts[level]):,}, {all_counts[level]/all_total*100:.1f}%)',va='center',fontsize=5.7)
    ins=fig.add_axes([.46,.31,.32,.43]); minor=['INV','TRA','DUP']; bot=np.zeros(3)
    for level in FREQ_LEVELS:
        vals=tab.loc[minor,level].to_numpy(); ins.bar(range(3),vals,bottom=bot,width=.65,color=FREQ_COLORS[level],edgecolor='white',linewidth=.25); bot+=vals
    for i,n in enumerate(bot): ins.text(i,n+max(bot)*.035,f'{int(n):,}',ha='center',fontsize=6)
    ins.set_xticks(range(3),minor); ins.tick_params(labelsize=5.5,length=2); ins.spines[['top','right']].set_visible(False); ins.set_ylim(0,max(bot)*1.22)
    save(fig,'Figure5b_hap38_frequency_spectrum')

def mm(k,a,b): return a*k/(b+k)
def extrap(x,y,maxx=100):
    x=np.asarray(x,float); y=np.maximum.accumulate(np.asarray(y,float)); anchor=x.max(); ay=y[-1]
    smooth=y.copy()
    for _ in range(2):
        p=np.pad(smooth,(1,1),mode='edge'); smooth=.22*p[:-2]+.56*p[1:-1]+.22*p[2:]; smooth[0]=y[0]; smooth[-1]=y[-1]
    xo=np.linspace(x.min(),anchor,240); yo=np.maximum.accumulate(PchipInterpolator(x,smooth)(xo)); yo[-1]=ay; xe=np.linspace(anchor,maxx,260)
    try:
        par,_=curve_fit(mm,x,y,p0=[ay*1.5,8],maxfev=10000); ye=ay+np.maximum(mm(xe,*par)-mm(np.array([anchor]),*par)[0],0)
    except Exception: ye=ay+(ay*1.2-ay)*(1-np.exp(-.04*(xe-anchor)))
    ye[0]=ay; return xo,yo,xe,np.maximum.accumulate(ye)

def exact_rarefaction(tier, n_total=38):
    """Expected loci observed after sampling n genomes, retaining exact SV types.

    The legacy saturation table groups TRA with several rare SyRI classes as
    ``Other``.  Figure 5c needs TRA alone, so compute the finite-population
    expectation directly from each locus' Sample_Count.
    """
    rows=[]
    denominator={n: math.comb(n_total,n) for n in range(1,n_total+1)}
    for typ in ['INS','DEL','INV','DUP','TRA']:
        counts=(tier.loc[tier.SVTYPE.eq(typ),'Sample_Count']
                    .astype(int).clip(1,n_total).value_counts())
        for n in range(1,n_total+1):
            expected=0.0
            for k,n_loci in counts.items():
                missed=0.0 if n_total-k<n else math.comb(n_total-k,n)/denominator[n]
                expected += float(n_loci)*(1.0-missed)
            rows.append({'SVTYPE_Group':typ,'N_Samples':n,
                         'Mean_Clusters':expected,'Method':'exact_hypergeometric_expectation'})
    return pd.DataFrame(rows)

def plot_c(tier):
    sat=exact_rarefaction(tier); sat.to_csv(TABLE/'Figure5c_rarefaction_source_hap38.tsv',sep='\t',index=False)
    fig,ax=plt.subplots(figsize=(6.4,4.45)); curves={}
    for t in ['INS','DEL','INV','DUP','TRA']:
        sub=sat[sat.SVTYPE_Group==t].sort_values('N_Samples'); curves[t]=extrap(sub.N_Samples,sub.Mean_Clusters)
        if t in {'INS','DEL'}:
            xo,yo,xe,ye=curves[t]; ax.plot(xo,yo/1000,color=COLORS[t],lw=2.0); ax.plot(xe,ye/1000,color=COLORS[t],lw=1.55,ls=(0,(4,3)),alpha=.78)
    ax.axvline(38,color='#aaa',lw=.7,ls=':'); ax.text(40,ax.get_ylim()[1]*.92,'n = 38',fontsize=7,color='#555'); ax.set_xlim(0,105); ax.set_xlabel('Number of genomes'); ax.set_ylabel('Number of SVs (x 10³)'); ax.spines[['top','right']].set_visible(False)
    for t in ['INS','DEL']:
        ax.text(101.5,curves[t][3][-1]/1000,t,color=COLORS[t],va='center',fontsize=7.2)
    ins=fig.add_axes([.54,.30,.34,.42])
    for t in ['INV','DUP','TRA']:
        xo,yo,xe,ye=curves[t]; ins.plot(xo,yo,color=COLORS[t],lw=1.35); ins.plot(xe,ye,color=COLORS[t],lw=1.1,ls=(0,(4,3)),alpha=.78); ins.text(101,ye[-1],t,color=COLORS[t],va='center',fontsize=5.7)
    ins.axvline(38,color='#bbb',lw=.5,ls=':'); ins.set_xlim(0,110); ins.tick_params(labelsize=5.4,length=2); ins.spines[['top','right']].set_visible(False); ins.set_xlabel('Genomes',fontsize=5.5); ins.set_ylabel('SVs (count)',fontsize=5.5)
    save(fig,'Figure5c_hap38_rarefaction')

def plot_d():
    d=pd.read_csv(SWAVE_PLOTDATA,sep='\t'); d.to_csv(TABLE/'Figure5d_SWave_plotdata.tsv',sep='\t',index=False)
    y=np.arange(len(d))[::-1]; break_l,break_r,gap=18,60,3
    mapx=lambda v: v if v<=break_l else break_l+gap+(v-break_r)
    fig=plt.figure(figsize=(5.15,3.05)); ax=fig.add_axes([.13,.16,.84,.78])
    for yy,row in zip(y,d.itertuples()):
        v=float(row.percent_all); c=row.color
        if v<=break_l: ax.barh(yy,v,color=c,height=.56)
        else:
            ax.barh(yy,break_l,color=c,height=.56); ax.barh(yy,v-break_r,left=break_l+gap,color=c,height=.56)
        ax.text(mapx(v)+.8,yy,f'{v:.1f}% ({int(row.count):,})',va='center',fontsize=7.0)
    ax.set_yticks(y,d.label); ax.set_xlim(0,mapx(72)+5); ax.set_xticks([0,5,10,15,mapx(60),mapx(65),mapx(70)],[0,5,10,15,60,65,70]); ax.tick_params(axis='y',length=0); ax.spines[['top','right','left']].set_visible(False)
    x0=break_l+gap/2-.35
    for off in (0,.42): ax.plot([x0+off,x0+.28+off],[-.02,.045],transform=ax.get_xaxis_transform(),color='#111',lw=.65,clip_on=False)
    inset=fig.add_axes([.54,.39,.35,.56]); inset.set_aspect('equal'); colors=d.color.tolist()
    inset.pie(d.percent_all,radius=1,colors=colors,startangle=90,counterclock=False,wedgeprops={'width':.22,'edgecolor':'white','linewidth':.8})
    inset.pie(d.high_percent_of_high,radius=.73,colors=colors,startangle=90,counterclock=False,wedgeprops={'width':.22,'edgecolor':'white','linewidth':.8})
    cds=d[d.region=='CDS'].iloc[0]; inset.text(0,.02,f'CDS\n{cds.percent_all:.1f} to {cds.high_percent_of_high:.1f}%',ha='center',va='center',fontsize=7.0,fontweight='bold'); inset.axis('off')
    save(fig,'Figure5d_SWave_complexSV')

def composite():
    names=['Figure5a_hap38_finalSV.png','Figure5b_hap38_frequency_spectrum.png','Figure5c_hap38_rarefaction.png','Figure5d_SWave_complexSV.png']
    imgs=[Image.open(PANELS/n).convert('RGB') for n in names]
    top_w=imgs[0].width; target=top_w//3
    bots=[]
    for im in imgs[1:]:
        ratio=target/im.width; bots.append(im.resize((target,int(im.height*ratio)),Image.Resampling.LANCZOS))
    bottom_h=max(i.height for i in bots)
    canvas=Image.new('RGB',(top_w,imgs[0].height+bottom_h),(255,255,255)); canvas.paste(imgs[0],(0,0))
    y0=imgs[0].height
    for i,im in enumerate(bots): canvas.paste(im,(i*target,y0))
    draw=ImageDraw.Draw(canvas)
    font_path=Path(mpl.get_data_path())/'fonts/ttf/DejaVuSans-Bold.ttf'
    label_font=ImageFont.truetype(str(font_path), max(42,top_w//150))
    pad=max(10,top_w//850)
    draw.text((pad,pad),'a',fill='black',font=label_font)
    draw.text((pad,y0+pad),'b',fill='black',font=label_font)
    draw.text((target+pad,y0+pad),'c',fill='black',font=label_font)
    draw.text((2*target+pad,y0+pad),'d',fill='black',font=label_font)
    canvas.save(FINAL/'Figure5_abcd_hap38_20260808.png',dpi=(600,600)); canvas.save(FINAL/'Figure5_abcd_hap38_20260808.pdf',resolution=600)

def main():
    tier,mixed,mem=load_inputs(); plot_a(mixed,mem); plot_b(tier); plot_c(tier); plot_d(); composite()
    print(f'[DONE] panels={PANELS} final={FINAL}')

if __name__=='__main__': main()
