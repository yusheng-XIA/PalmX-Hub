import fitz, numpy as np, pandas as pd
d=fitz.open('../../deliver/Main_Figures_revised/Figure2.pdf'); p=d[0]
dr=p.get_drawings()
S='../../deliver/Source_Data_split/Source_Data_Fig2.xlsx'
sd=pd.read_excel(S,'Fig2d_axis_scores')
teal=(0.09799344092607498, 0.49799343943595886, 0.5569848418235779)
lines=[x for x in dr if x['type']=='s' and x['rect'].x0>280 and x['rect'].y0>285 and x['rect'].y1<375 and len(x['items'])==18]
bands=[x for x in dr if x['type']=='f' and x['rect'].x0>280 and x['rect'].y0>285 and x['rect'].y1<375 and len(x['items'])==38]
panel={(296,297):'P02',(410,297):'P01',(296,343):'P03',(410,343):'P04'}
def pid(r):
    return 'P02' if r.x0<300 and r.y0<330 else 'P01' if r.y0<330 else 'P03' if r.x0<300 else 'P04'
out=[]
for ax in ['P02','P01','P03','P04']:
    L={('FL' if l['color']==teal else 'TN'):l for l in lines if pid(l['rect'])==ax}
    B={('FL' if b['fill']==teal else 'TN'):b for b in bands if pid(b['rect'])==ax}
    pts={}
    for g,l in L.items():
        pp=[l['items'][0][1]]+[it[2] for it in l['items']]
        pts[g]=np.array([(q.x,q.y) for q in pp])
    sub=sd[sd.axis_id==ax].sort_values('stage_index')
    fl=sub[sub.genotype=='FL'].set_index('stage_index'); tn=sub[sub.genotype=='TN'].set_index('stage_index')
    yall=np.r_[pts['FL'][:,1],pts['TN'][:,1]]; vall=np.r_[fl['mean'].values,tn['mean'].values]
    A=np.c_[vall,np.ones_like(vall)]; coef,res,_,_=np.linalg.lstsq(A,yall,rcond=None)
    pred=A@coef; resid_pt=np.abs(pred-yall)
    # convert residual to data units
    resid=resid_pt/abs(coef[0])
    # swapped-genotype check
    yswap=np.r_[pts['TN'][:,1],pts['FL'][:,1]]; c2,_,_,_=np.linalg.lstsq(A,yswap,rcond=None); rs=np.abs(A@c2-yswap).max()/abs(c2[0])
    # band check: se half-width
    bw={}
    for g,b in B.items():
        pp=[b['items'][0][1]]+[it[2] for it in b['items'] if it[0]=='l']
        ys=np.array([q.y for q in pp])
        bw[g]=ys
    print(ax,'scale pt/unit',round(coef[0],2),'max resid (data units)',round(resid.max(),4),'swapped max resid',round(rs,3), 'x monotone', np.all(np.diff(pts['FL'][:,0])>0))
    for g,df in [('FL',fl),('TN',tn)]:
        ys=bw[g]; upper=ys[:19]; lower=ys[19:38][::-1] if len(ys)>=38 else None
        if lower is not None:
            hw=np.abs(upper-lower)/2/abs(coef[0])
            print('  ',g,'band half-width vs SE max diff',round(np.nanmax(np.abs(hw-df['se'].values)),4))
