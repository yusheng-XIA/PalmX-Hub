import pandas as pd, numpy as np
w=pd.read_csv('fig2_words.txt',sep='\t',header=None,names=['x','y','t'],dtype={'t':str},keep_default_na=False)
e=w[(w.y>405)&(w.y<560)&(w.x>90)&(w.x<250)]
# column headers
hdr=w[(w.y>395)&(w.y<405)&(w.x>90)]
print(hdr.to_string())
labs=w[(w.y>405)&(w.y<560)&(w.x>70)&(w.x<92)].sort_values('y')
vals=e[e.t.str.match(r'^[0-9.]+$|^NA$')].copy()
ys=sorted(vals.y.unique()); rows=[]
for y in ys:
    if not rows or y-rows[-1][-1]>3: rows.append([y])
    else: rows[-1].append(y)
print(len(rows), len(vals))
markers=['ACCase','KASIII','ENR','FATA/B','LACS','GPAT','FAD2-like','DGAT','OLE16']
markers=['ACCase','KASIII','ENR','FATA/B','LACS','GPAT','FAD2-like','DGAT','OLE16']
cols=sorted(vals.x.round(-1).unique()); print(cols)
st=['125d','170d','185d','12h','36h']
fig={}
for i,y in enumerate(rows):
    mk=markers[i//2]; lay='RNA' if i%2==0 else 'Protein'
    r=vals[vals.y.isin(y)].sort_values('x')
    for j,(xx,t) in enumerate(zip(r.x,r.t)): fig[(mk,lay,st[j])]=t
rna=pd.read_csv('out/fig2e_rna_recomputed.tsv',sep='\t')
pr=pd.read_csv('out/fig2e_protein_recomputed.tsv',sep='\t')
pmap={'ENR':'FabI (ENR)','FAD2-like':'FAD2'}
out=[]
for (mk,lay,s),t in fig.items():
    if lay=='RNA':
        v=rna[(rna.marker==mk)&(rna.stage==s)].FL_over_TN.iloc[0]
        det=''
    else:
        f=pmap.get(mk,mk); q=pr[(pr.family==f)&(pr.stage==s)].set_index('genotype')
        det=f"{q.loc['FL','detected_reps']}/{q.loc['TN','detected_reps']}"
        v=q.loc['FL','mean_detected']/q.loc['TN','mean_detected']
    rv='NA' if not np.isfinite(v) else (f'{v:.1f}' if v>=10 else f'{v:.2f}')
    out.append(dict(marker=mk,layer=lay,stage=s,figure=t,recomputed=v,recomputed_rounded=rv,prot_detected_reps_FL_TN=det,match=(t==rv)))
o=pd.DataFrame(out); o.to_csv('fig2e_figure_vs_recomputed.tsv',sep='\t',index=False)
print(o.match.value_counts()); print(o[~o.match])
# extra text checks
f=rna[rna.marker=='FAD2-like']; print('FAD2 RNA FL>TN stages:',(f.FL>f.TN).sum(), f[f.FL<=f.TN][['stage','FL_over_TN']].to_string())
d=pr.pivot_table(index=['family','genotype'],values='detected_reps',aggfunc=lambda x:(x>=2).sum()); print(d.unstack())
o16=pr[pr.family=='OLE16']; print(o16.loc[o16.groupby('genotype').mean_detected.idxmax()][['genotype','stage']])
print(pr[(pr.stage=='185d')&pr.family.isin(['OLE16','DGAT'])])
