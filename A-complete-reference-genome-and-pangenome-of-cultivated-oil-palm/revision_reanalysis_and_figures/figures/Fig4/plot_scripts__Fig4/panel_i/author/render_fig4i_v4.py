#!/usr/bin/env python3
"""Render-only refinement for Figure4i: exact canvas and unclipped endpoint labels."""
from pathlib import Path
import json,hashlib,platform,sys
from datetime import datetime,timezone
import os
os.environ.setdefault('MPLBACKEND','Agg');os.environ.setdefault('MPLCONFIGDIR','/tmp/mpl-fig4i-render-v2')
import cairosvg,matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np,pandas as pd
from scipy.signal import savgol_filter
RUN=Path('${ANALYSIS_DIR}/22_answer_reviews/00_ms/05_MS/0918_revision/final_ms/02_figure4/runs/RUN-FIG4I-39ASSEMBLY-20260923-001')
SRC=RUN/'outputs/TE_density_class_position_curve_genome39_material33.tsv'
OUT=RUN/'render_attempt4';OUT.mkdir(exist_ok=True)
CLASSES=('Core','Soft-core','Shell','Cloud');COLORS={'Core':'#82C7B8','Soft-core':'#E65D6D','Shell':'#29AFD4','Cloud':'#7DCDF8'}
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
 return h.hexdigest()
def smooth(v):
 r=np.asarray(v,float).copy()
 for a,b,w in ((0,20,7),(20,120,11),(120,140,7)):r[a:b]=savgol_filter(r[a:b],w,2,mode='interp')
 return r
f=pd.read_csv(SRC,sep='\t');assert len(f)==560
plt.rcParams.update({'font.family':'sans-serif','font.sans-serif':['Liberation Sans','Arial','DejaVu Sans'],'font.size':5.5,'axes.labelsize':6.5,'axes.linewidth':.65,'svg.fonttype':'none'})
fig,ax=plt.subplots(figsize=(60/25.4,38.5/25.4),facecolor='white')
for c in CLASSES:
 d=f[f.Pangenome_class.eq(c)].sort_values('Profile_bin');assert len(d)==140
 x=d.Profile_bin.to_numpy(float)+.5;y=smooth(d.Mean_material_TE_percent);lo=np.clip(smooth(d.CI95_low),0,100);hi=np.clip(smooth(d.CI95_high),0,100)
 ax.plot(x,y,color=COLORS[c],lw=.9,label=c,zorder=3,solid_capstyle='round',solid_joinstyle='round');ax.fill_between(x,np.minimum(lo,y),np.maximum(hi,y),color=COLORS[c],alpha=.16,linewidth=0)
ax.axvline(20,color='#333333',ls=(0,(4,3)),lw=.55);ax.axvline(120,color='#333333',ls=(0,(4,3)),lw=.55)
ax.set_xlim(0,140);ax.set_ylim(0,42.5);ax.set_yticks(np.arange(0,41,5));ax.set_xticks([0,20,120,140]);ax.set_xticklabels(['−2 kb','TSS','TES','+2 kb'])
xt=ax.get_xticklabels();xt[-1].set_ha('right')
ax.set_ylabel('TE density (%)');ax.legend(frameon=False,ncol=1,loc='upper center',bbox_to_anchor=(.56,1.0),handlelength=1.8,fontsize=5.2)
ax.spines['top'].set_visible(False);ax.spines['right'].set_visible(False);ax.tick_params(width=.55,length=2,labelsize=5.3)
fig.subplots_adjust(left=.155,right=.975,bottom=.19,top=.95);fig.text(.012,.985,'i',ha='left',va='top',fontsize=9.5,fontweight='bold')
svg=OUT/'Figure4i_39assembly_classes_60x38.5mm_vector.svg';pdf=OUT/'Figure4i_39assembly_classes_60x38.5mm_vector.pdf';png=OUT/'Figure4i_39assembly_classes_60x38.5mm_600dpi.png'
fig.savefig(svg,facecolor='white');fig.savefig(png,dpi=600,facecolor='white');plt.close(fig)
cairosvg.svg2pdf(url=str(svg),write_to=str(pdf))
rec={'status':'RENDER_COMPLETE_PENDING_AUDIT','parent_run':'RUN-FIG4I-39ASSEMBLY-20260923-001','operation':'render-only exact-size correction v4; scientific table unchanged','completed_utc':datetime.now(timezone.utc).isoformat(),'host':platform.node(),'python':sys.version,'source_table':str(SRC),'source_sha256':sha(SRC),'outputs':{str(x):sha(x) for x in (svg,pdf,png)}}
(OUT/'render_record.json').write_text(json.dumps(rec,indent=2)+'\n');print(json.dumps(rec,indent=2))
