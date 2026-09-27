import pandas as pd, numpy as np
pd.set_option('display.width',250)
ST=["0d","15d","35d","50d","65d","80d","95d","110d","125d","140d","155d","170d","185d","12h","24h","36h","48h","60h","72h"]
L=pd.read_csv('res/ld_precursor_rows.tsv.gz',sep='\t')
Q=pd.read_csv('res/run_precursor_quantiles.tsv',sep='\t',index_col=0)
L[['mat','stage','rep']]=L.key.str.split('|',expand=True)
GRP={'VPGSEQLEQAR':'OLE16a_chr11_shared','RVPGSEQLEQAR':'OLE16a_chr11_shared',
     'RPPGFEQLEQAR':'OLE16b_chr04_shared','RPPGSEQLEQAR':'OLE16b_chr04_FLhap1only'}
O=L[L['Stripped.Sequence'].isin(GRP)&L.accepted].copy(); O['grp']=O['Stripped.Sequence'].map(GRP)
print(O.groupby(['Stripped.Sequence','Precursor.Id','mat']).size().unstack())
# detection table per precursor x stage (n reps detected)
det=O.groupby(['Precursor.Id','mat','stage']).key.nunique().unstack('stage').reindex(columns=ST).fillna(0).astype(int)
print(det.to_string())
det.to_csv('res/ole16_precursor_detection.tsv',sep='\t')
# per-sample sum by group
S=O.groupby(['grp','key','mat','stage'])['Precursor.Quantity'].sum().reset_index()
S.to_csv('res/ole16_group_sample_sums.tsv',sep='\t',index=False)
M=S.groupby(['grp','mat','stage'])['Precursor.Quantity'].agg(['count','mean']).unstack('stage')
print(M['count'].reindex(columns=ST).to_string()); 
mm=M['mean'].reindex(columns=ST)
print((mm/1e6).round(2).to_string())
# FL/TN ratios for shared chr11 precursors, per precursor, detected in both
rr=[]
for p,g in O[O.grp=='OLE16a_chr11_shared'].groupby('Precursor.Id'):
    for s in ST:
        f=g[(g.mat=='FL')&(g.stage==s)]['Precursor.Quantity']; t=g[(g.mat=='TN')&(g.stage==s)]['Precursor.Quantity']
        rr.append(dict(prec=p,stage=s,nFL=len(f),nTN=len(t),FL=f.mean() if len(f) else np.nan,TN=t.mean() if len(t) else np.nan))
R=pd.DataFrame(rr); R['ratio']=R.FL/R.TN; print(R[R.stage.isin(ST[11:])].to_string())
R.to_csv('res/ole16_chr11_precursor_ratio.tsv',sep='\t',index=False)
# TN intensity vs run percentiles
T=O[(O.mat=='TN')].merge(Q,left_on='key',right_index=True)
T['vs_p05']=T['Precursor.Quantity']/T.p05; T['vs_p50']=T['Precursor.Quantity']/T.p50
print(T[T.stage.isin(ST[11:])].groupby(['Precursor.Id','stage'])[['Precursor.Quantity','p05','p50','vs_p05','vs_p50']].median().to_string())
# percentile of each OLE16 precursor within run distribution: approximate using quantiles table -> need full; report ratios
# unaccepted rows of OLE16 in TN (sub-threshold evidence)
U=L[L['Stripped.Sequence'].isin(GRP)&~L.accepted]; print('unaccepted OLE16 rows',len(U)); print(U.groupby(['mat','stage']).size())
