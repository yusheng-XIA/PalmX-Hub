#!/usr/bin/env python3
"""Recompute Figure 4i with 39-assembly occupancy classes and 33-material averaging.

Immutable scientific inputs are reused from the audited Figure 4 workspace.
Family classes are defined exactly as Figure 4e/f: Core=39/39,
Soft-core=38/39, Shell=2-37/39, Cloud=1/39. TE profiles are calculated per
assembly, phased haplotypes are folded to 33 biological materials, and material
means receive equal weight in the displayed profile.
"""
from __future__ import annotations
import csv, gzip, hashlib, importlib.util, json, os, platform, shutil, subprocess, sys, warnings
from collections import Counter
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime, timezone
from pathlib import Path
os.environ.setdefault('OPENBLAS_NUM_THREADS','1'); os.environ.setdefault('OMP_NUM_THREADS','1'); os.environ.setdefault('MKL_NUM_THREADS','1')
os.environ.setdefault('MPLBACKEND','Agg'); os.environ.setdefault('MPLCONFIGDIR','/tmp/matplotlib-fig4i-genome39')
import cairosvg
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy
from scipy.signal import savgol_filter

RUN=Path(os.environ['FIG4I_RUN'])
OUT=RUN/'outputs'; STAGE=RUN/'input_stage'; SHARDS=RUN/'shards'; PROV=RUN/'provenance'; LOGS=RUN/'logs'
BASE=Path('${ANALYSIS_DIR}/22_answer_reviews/00_ms/03_V3/04_figure4/Fig4_d_i_pan39_material33_singletons_20260811')
SOURCE_RUN=Path('${ANALYSIS_DIR}/22_answer_reviews/00_ms/03_V3/04_figure4/Fig4_d_i_pan39_material33_redraw_20260808')
SOURCE_SCRIPT=SOURCE_RUN/'scripts/05_run_i_te_material33.py'
MEMBERS=BASE/'tables/Orthogroups.members.GO_singletons_64577.tsv'
PAV39=BASE/'pan39_genome_redraw_20260812/tables/Orthogroups.PAV.genome39.GO_singletons_64577.tsv'
MANIFEST=Path('${ANALYSIS_DIR}/14_pan_genome/11_new_pan/config/sample_manifest.tsv')
CLASSES=('Core','Soft-core','Shell','Cloud')
COLORS={'Core':'#82C7B8','Soft-core':'#E65D6D','Shell':'#29AFD4','Cloud':'#7DCDF8'}
WORKERS=int(os.environ.get('FIG4I_WORKERS','4'))

def sha256(path:Path):
 h=hashlib.sha256()
 with path.open('rb') as f:
  for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
 return h.hexdigest()

def import_source():
 name=f'fig4i_source_{os.getpid()}'
 spec=importlib.util.spec_from_file_location(name,SOURCE_SCRIPT)
 if spec is None or spec.loader is None:raise ImportError(SOURCE_SCRIPT)
 m=importlib.util.module_from_spec(spec);sys.modules[name]=m;spec.loader.exec_module(m);return m

def load_class_map():
 f=pd.read_csv(PAV39,sep='\t')
 if len(f)!=64577 or f['Orthogroup'].duplicated().any():raise AssertionError('PAV39 identity failed')
 freq=f.iloc[:,1:].sum(axis=1).astype(int)
 cls=np.select([freq.eq(39),freq.eq(38),freq.between(2,37)],['Core','Soft-core','Shell'],default='Cloud')
 counts=pd.Series(cls).value_counts().reindex(CLASSES).astype(int).to_dict()
 expected={'Core':19305,'Soft-core':1879,'Shell':27450,'Cloud':15943}
 if counts!=expected:raise AssertionError(f'class counts {counts} != {expected}')
 return dict(zip(f['Orthogroup'],cls)),counts

def process_one(payload):
 assembly,source_info,new_samples=payload
 warnings.filterwarnings('ignore',message='Mean of empty slice')
 source=import_source();source.RUN=RUN;source.OUT=OUT;source.STAGE=STAGE;source.MEMBER_FILE=MEMBERS
 original=source.import_original()
 original.load_edta_gff=lambda path,allowed:source.load_edta_gff_order_independent(original,path,allowed)
 original.PREP=STAGE;original.MEMBER_FILE=MEMBERS
 original_loader=original.load_members_for_assembly
 def member_loader(path,selected_assembly,class_mapping):
  recs=original_loader(path,selected_assembly,class_mapping);result=[];prefix=selected_assembly+'__'
  for gene_id,og,cat in recs:
   if not gene_id.startswith(prefix):raise ValueError(f'Unexpected member prefix {gene_id}')
   gene_id=gene_id[len(prefix):]
   if selected_assembly in new_samples and gene_id.startswith('evm.model'):gene_id=gene_id.replace('evm.model','evm.TU',1)
   result.append((gene_id,og,cat))
  return result
 original.load_members_for_assembly=member_loader;original.source_for=lambda selected_assembly:source_info
 class_map,_=load_class_map()
 shard=SHARDS/assembly;shard.mkdir(parents=True,exist_ok=True)
 gene_path=shard/f'{assembly}.gene_level.tsv.gz'
 with gzip.open(gene_path,'wt',newline='',compresslevel=4) as handle:
  w=csv.writer(handle,delimiter='\t',lineterminator='\n')
  w.writerow(['Assembly','Orthogroup','Gene_ID','Pangenome_class','Sequence','Gene_start_0based','Gene_end_0based_exclusive','Strand','Gene_length_bp','Upstream_effective_bp','Downstream_effective_bp','Upstream_2kb_TE_percent','Gene_body_TE_percent','Downstream_2kb_TE_percent','Whole_gene_plus_flanks_TE_percent','TE_source_method'])
  qa,summaries,curves=original.process_assembly(assembly,class_map,w)
 if len(summaries)!=4 or len(curves)!=560:raise AssertionError(f'incomplete {assembly}: {len(summaries)}, {len(curves)}')
 (shard/'SUCCESS.json').write_text(json.dumps({'assembly':assembly,'qa':qa,'gene_level_sha256':sha256(gene_path)},indent=2)+'\n')
 return assembly,qa,summaries,curves,str(gene_path),sha256(gene_path)

def smooth(v):
 r=np.asarray(v,float).copy()
 for a,b,w in ((0,20,7),(20,120,11),(120,140,7)):r[a:b]=savgol_filter(r[a:b],window_length=w,polyorder=2,mode='interp')
 return r

def plot_final(frame):
 plt.rcParams.update({'font.family':'sans-serif','font.sans-serif':['Liberation Sans','Arial','DejaVu Sans'],'font.size':5.5,'axes.labelsize':6.5,'axes.linewidth':0.65,'pdf.fonttype':42,'ps.fonttype':42,'svg.fonttype':'none'})
 fig,ax=plt.subplots(figsize=(60/25.4,38.5/25.4),facecolor='white')
 for c in CLASSES:
  d=frame.loc[frame.Pangenome_class.eq(c)].sort_values('Profile_bin')
  if len(d)!=140:raise AssertionError(f'{c}: {len(d)} bins')
  x=d.Profile_bin.to_numpy(float)+.5;y=smooth(d.Mean_material_TE_percent);lo=np.clip(smooth(d.CI95_low),0,100);hi=np.clip(smooth(d.CI95_high),0,100)
  lo=np.minimum(lo,y);hi=np.maximum(hi,y)
  ax.plot(x,y,color=COLORS[c],lw=.9,label=c,zorder=3,solid_capstyle='round',solid_joinstyle='round')
  ax.fill_between(x,lo,hi,color=COLORS[c],alpha=.16,linewidth=0)
 ax.axvline(20,color='#333333',ls=(0,(4,3)),lw=.55);ax.axvline(120,color='#333333',ls=(0,(4,3)),lw=.55)
 ax.set_xlim(0,140);ax.set_ylim(bottom=0);ax.set_xticks([0,20,120,140]);ax.set_xticklabels(['−2 kb','TSS','TES','+2 kb'])
 ax.set_ylabel('TE density (%)');ax.legend(frameon=False,ncol=1,loc='upper center',bbox_to_anchor=(.56,1.0),handlelength=1.8,fontsize=5.2)
 ax.spines['top'].set_visible(False);ax.spines['right'].set_visible(False);ax.tick_params(width=.55,length=2,labelsize=5.3)
 fig.subplots_adjust(left=.155,right=.975,bottom=.19,top=.95);fig.text(.012,.985,'i',ha='left',va='top',fontsize=9.5,fontweight='bold')
 svg=OUT/'Figure4i_39assembly_classes_60x38.5mm_vector.svg';pdf=OUT/'Figure4i_39assembly_classes_60x38.5mm_vector.pdf';png=OUT/'Figure4i_39assembly_classes_60x38.5mm_600dpi.png'
 fig.savefig(svg,facecolor='white');fig.savefig(png,dpi=600,facecolor='white');plt.close(fig)
 cairosvg.svg2pdf(url=str(svg),write_to=str(pdf),output_width=60/25.4*72,output_height=38.5/25.4*72)
 return svg,pdf,png

def main():
 started=datetime.now(timezone.utc)
 for x in (OUT,STAGE,SHARDS,PROV,LOGS):x.mkdir(parents=True,exist_ok=True)
 class_map,class_counts=load_class_map()
 source=import_source();source.RUN=RUN;source.OUT=OUT;source.STAGE=STAGE;source.MEMBER_FILE=MEMBERS
 original=source.import_original();design,material_map=source.load_design();source_map,new_samples=source.prepare_stage(original,design)
 assemblies=[r['sample_id'] for r in design]
 completed={};payloads=[(a,source_map[a],new_samples) for a in assemblies]
 with ProcessPoolExecutor(max_workers=WORKERS) as ex:
  fs={ex.submit(process_one,p):p[0] for p in payloads}
  for n,f in enumerate(as_completed(fs),1):
   a,qa,summaries,curves,gp,gh=f.result();completed[a]=(qa,summaries,curves,gp,gh);print(f'[{n:02d}/39] {a}',flush=True)
 qa_rows=[];summary_rows=[];curve_rows=[];gene_manifest=[]
 for a in assemblies:
  qa,s,c,gp,gh=completed[a];qa_rows.append(qa);summary_rows.extend(s);curve_rows.extend(c);gene_manifest.append((gh,gp))
 qa=pd.DataFrame(qa_rows);hap_summary=pd.DataFrame(summary_rows);hap_curve=pd.DataFrame(curve_rows)
 hap_summary['Material']=hap_summary.Assembly.map(material_map)
 numeric=[c for c in hap_summary.columns if c not in {'Assembly','Material','Pangenome_class'}]
 material_summary=hap_summary.groupby(['Material','Pangenome_class'],sort=False)[numeric].mean().reset_index();material_summary['Haplotype_count']=material_summary.Material.map(Counter(material_map.values()))
 material_curve,final_curve=source.summarize_material_curve(hap_curve,material_map)
 qa.to_csv(OUT/'TE_density_source_and_coordinate_QA_genome39.tsv',sep='\t',index=False)
 hap_summary.to_csv(OUT/'TE_density_haplotype_class_summary_genome39.tsv',sep='\t',index=False)
 material_summary.to_csv(OUT/'TE_density_material_class_summary_genome39.tsv',sep='\t',index=False)
 hap_curve.to_csv(OUT/'TE_density_haplotype_class_position_curve_genome39.tsv.gz',sep='\t',index=False,compression='gzip')
 material_curve.to_csv(OUT/'TE_density_material_class_position_curve_genome39.tsv.gz',sep='\t',index=False,compression='gzip')
 final_curve.to_csv(OUT/'TE_density_class_position_curve_genome39_material33.tsv',sep='\t',index=False)
 source.write_statistics(material_summary)
 # Rename source statistics to make the 39-class identity explicit.
 for old,new in [('TE_density_primary_paired_inference_material33.tsv','TE_density_primary_paired_inference_genome39_material33.tsv'),('TE_density_sensitivity_Kruskal_Wallis_material33.tsv','TE_density_sensitivity_Kruskal_Wallis_genome39_material33.tsv'),('TE_density_sensitivity_Mann_Whitney_Holm_material33.tsv','TE_density_sensitivity_Mann_Whitney_Holm_genome39_material33.tsv')]:
  (OUT/old).rename(OUT/new)
 svg,pdf,png=plot_final(final_curve)
 completed_at=datetime.now(timezone.utc)
 checks={'status':'EXECUTION_COMPLETE_PENDING_INDEPENDENT_AUDIT','started_utc':started.isoformat(),'completed_utc':completed_at.isoformat(),'elapsed_seconds':(completed_at-started).total_seconds(),'host':platform.node(),'python':sys.version,'pandas':pd.__version__,'numpy':np.__version__,'scipy':scipy.__version__,'matplotlib':matplotlib.__version__,'cairosvg':cairosvg.__version__,'workers':WORKERS,'family_class_definition':'39 assemblies: Core=39, Soft-core=38, Shell=2-37, Cloud=1','family_class_counts':class_counts,'haplotypes':len(qa),'materials':material_summary.Material.nunique(),'haplotype_class_rows':len(hap_summary),'material_class_rows':len(material_summary),'haplotype_curve_rows':len(hap_curve),'material_curve_rows':len(material_curve),'final_curve_rows':len(final_curve),'min_mapping_rate':float(qa.Mapping_rate.min()),'positive_te_haplotypes':int(qa.repeat_bp_merged.gt(0).sum()),'display_smoothing':'Savitzky-Golay order 2, windows 7/11/7, regions independent','outputs':{str(x):sha256(x) for x in [svg,pdf,png]},'gene_level_shards':dict(gene_manifest)}
 expected={'haplotypes':39,'materials':33,'haplotype_class_rows':156,'material_class_rows':132,'haplotype_curve_rows':21840,'material_curve_rows':18480,'final_curve_rows':560,'positive_te_haplotypes':39}
 failures=[f'{k}={checks[k]} expected {v}' for k,v in expected.items() if checks[k]!=v]
 if checks['min_mapping_rate']<.99:failures.append('mapping rate <0.99')
 if failures:checks['status']='FAILED';checks['failures']=failures
 (PROV/'execution_record.json').write_text(json.dumps(checks,indent=2)+'\n')
 with (PROV/'gene_level_output_manifest.sha256.tsv').open('w') as h:
  h.write('sha256\tpath\n');[h.write(f'{d}\t{p}\n') for d,p in gene_manifest]
 if failures:raise RuntimeError('; '.join(failures))
 print(json.dumps({k:v for k,v in checks.items() if k!='gene_level_shards'},indent=2))
if __name__=='__main__':main()
