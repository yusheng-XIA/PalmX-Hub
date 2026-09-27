import fitz, numpy as np, pandas as pd
d=fitz.open('../../deliver/Main_Figures_revised/Figure2.pdf'); p=d[0]
dr=p.get_drawings()
cb=[x for x in dr if x['rect'].x0>480 and x['rect'].x1<500 and 219<x['rect'].y0<251 and x.get('fill')]
cb=sorted(cb,key=lambda x:x['rect'].y0)
cbcol=np.array([x['fill'] for x in cb]); n=len(cbcol)
cbval=np.linspace(2,-2,n)
na=(0.7249866724014282, 0.7449912428855896, 0.764995813369751)
cells=[x for x in dr if x['rect'].x0>320 and x['rect'].x1<485 and 213<x['rect'].y0<262 and x.get('fill') and (x['rect'].x1-x['rect'].x0)<10]
xs=sorted(set(round(c['rect'].x0,0) for c in cells)); ys=sorted(set(round(c['rect'].y0,0) for c in cells))
print(len(cells),len(xs),len(ys))
names=['TG 42:0','PA 36:2','FA 18:1+3O','FA 18:3','LPE 18:0','Azelaic acid','Hyperoside','P-Coumaric acid']
st=["0d","15d","35d","50d","65d","80d","95d","110d","125d","140d","155d","170d","185d","12h","24h","36h","48h","60h","72h"]
rows=[]
for c in cells:
    i=ys.index(round(c['rect'].y0,0)); j=xs.index(round(c['rect'].x0,0))
    f=np.array(c['fill'])
    if np.abs(f-na).max()<0.01: v=np.nan
    else:
        k=np.argmin(((cbcol-f)**2).sum(1)); v=cbval[k]; dist=np.sqrt(((cbcol[k]-f)**2).sum())
    rows.append((names[i],st[j],v))
fig=pd.DataFrame(rows,columns=['display_name','stage','fig_delta'])
S='../../deliver/Source_Data_split/Source_Data_Fig2.xlsx'
sd=pd.read_excel(S,'Fig2c_metabolites')
print(sd.display_name.unique())
pv=sd.pivot_table(index=['display_name','stage'],columns='genotype',values='mean_consensus_z',dropna=False)
pv['sd_delta']=(pv.FL-pv.TN).clip(-2,2)
m=fig.merge(pv.reset_index(),on=['display_name','stage'],how='left')
m['diff']=m.fig_delta-m.sd_delta
print(m['diff'].abs().describe())
print(m[m['diff'].abs()>0.1])
print('NA cells fig',m.fig_delta.isna().sum(),'NA sd',m.sd_delta.isna().sum())
m.to_csv('fig2c_pdf_vs_sd.tsv',sep='\t',index=False)
