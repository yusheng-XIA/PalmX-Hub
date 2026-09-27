"""Build corrected SF7 Source Data tables: swap 50 d primer 3 wells F7/G7 back to instrument CSV values."""
import pandas as pd, numpy as np
X='../../deliver/Source_Data_split/Source_Data_Supplementary_Figs.xlsx'
d=pd.read_excel(X,sheet_name=None)
raw=d['SF7_raw_dCt'].copy()
L='log₂(target/ACTIN) = −ΔCt'
fix={('50d','FL'):[30.75,27.65,28.2,27.93],   # CSV G5-G8
     ('50d','TN'):[27.36,28.21,27.07,27.42]}  # CSV F5-F8
for (st,mt),v in fix.items():
    i=raw.index[(raw.Stage==st)&(raw.Material==mt)&(raw['Primer assay']==3)]
    assert len(i)==1; i=i[0]
    assert raw.at[i,'Finite target-gene wells']=='1,2,3,4'
    raw.at[i,'Target-gene Ct values']=','.join(f'{x:g}' for x in v)
    mt_=float(np.mean(v)); raw.at[i,'Mean target-gene Ct']=mt_
    dct=mt_-raw.at[i,'Mean ACTIN Ct']
    raw.at[i,'ΔCt (target − ACTIN)']=dct; raw.at[i,L]=-dct; raw.at[i,'target/ACTIN = 2^(−ΔCt)']=2**(-dct)
# stage summary: update primer-3 columns + means from raw
ss=d['SF7_stage_summary'].copy()
for i,r in ss.iterrows():
    for a in (1,2,3):
        rr=raw[(raw.Stage==r.Stage)&(raw.Material==r.Material)&(raw['Primer assay']==a)].iloc[0]
        for col,src in [(f'Primer {a} log₂(target/ACTIN)',L),(f'Primer {a} target/ACTIN','target/ACTIN = 2^(−ΔCt)')]:
            if abs(ss.at[i,col]-rr[src])>1e-12: ss.at[i,col]=rr[src]
    m=np.mean([ss.at[i,f'Primer {a} log₂(target/ACTIN)'] for a in (1,2,3)])
    if abs(ss.at[i,'Mean log₂(target/ACTIN)']-m)>1e-12:
        ss.at[i,'Mean log₂(target/ACTIN)']=m; ss.at[i,'Geometric mean target/ACTIN']=2**m
for n,t in [('SF7_raw_dCt',raw),('SF7_stage_summary',ss)]:
    t.to_csv(f'SD_{n}.tsv',sep='\t',index=False)
    o=d[n]; diff=[(i,c,o.at[i,c],t.at[i,c]) for i in t.index for c in t.columns
                  if not (o.at[i,c]==t.at[i,c] or (isinstance(o.at[i,c],float) and abs(o.at[i,c]-t.at[i,c])<1e-12))]
    print(n, 'changed cells:', len(diff))
    for i,c,a,b in diff: print(f'  row{i+2} {t.at[i,"Stage"]} {t.at[i,"Material"]} | {c}: {a} -> {b}')
