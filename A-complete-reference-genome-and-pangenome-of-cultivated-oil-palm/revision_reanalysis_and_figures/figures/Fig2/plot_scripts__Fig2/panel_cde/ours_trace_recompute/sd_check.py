import pandas as pd, numpy as np
S='${WORK_DIR}/deliver/Source_Data_split/Source_Data_Fig2.xlsx'
rna=pd.read_csv('up/Figure2f_RNA_latest114_joint_source.tsv',sep='\t')
prot=pd.read_csv('up/Figure2f_Astral114_stage_detection.tsv',sep='\t')
sdr=pd.read_excel(S,'Fig2e_RNA_114'); sdp=pd.read_excel(S,'Fig2e_protein_114')
print('RNA markers',rna.marker.unique()); print('prot fam',prot.family.unique())
st=['125d','170d','185d','12h','36h']
m=rna.merge(sdr,on=['marker','genotype','stage'],suffixes=('','_sd'))
print('RNA SD vs up max absdiff',(m.raw_aggregate-m.raw_aggregate_sd).abs().max(),len(m),len(sdr))
m=prot.merge(sdp,on=['family','genotype','stage'],suffixes=('','_sd'))
print('prot SD vs up max reldiff',((m.mean_detected_replicate_abundance-m.mean_detected_replicate_abundance_sd)/m.mean_detected_replicate_abundance).abs().max(),len(m),len(sdp))
p=rna.pivot_table(index=['marker','stage'],columns='genotype',values='raw_aggregate')
r=(p.FL/p.TN).unstack()[st]; print('RNA FL/TN\n',r.round(2))
q=prot.pivot_table(index=['family','stage'],columns='genotype',values='mean_detected_replicate_abundance')
print(prot[prot.stage.isin(st)].pivot_table(index=['family','stage'],columns='genotype',values='detected_replicate_n').unstack())
r2=(q.FL/q.TN).unstack()[st]; print('Prot FL/TN\n',r2.round(2))
