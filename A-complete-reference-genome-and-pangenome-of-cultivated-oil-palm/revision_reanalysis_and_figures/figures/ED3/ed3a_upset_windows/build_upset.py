"""ED3a UpSet from gene-stage robust ASE calls; phase windows as argument. Validates the published windows first."""
import pandas as pd, numpy as np, sys
u=pd.read_csv('../../ase_genogroup/in/gene_stage_ASE_unified.tsv.gz',sep='\t',usecols=['analysis','gene_id','stage','eligible','robust_ase','ase_call'],low_memory=False)
for c in ['eligible','robust_ase']: u[c]=u[c].astype(str).str.lower().eq('true')
ov=pd.read_csv('in/ASE_overlap_one_to_one_current.tsv',sep='\t')
OLD={'Early':['0d','15d','35d','50d'],'Mid':['65d','80d','95d','110d'],'Late':['125d','140d','155d','170d','185d'],'Postharvest':['12h','24h','36h','48h','60h','72h']}
NEW={'Early':['0d','15d','35d','50d','65d'],'Mid':['80d','95d','110d','125d','140d'],'Late':['155d','170d','185d'],'Postharvest':['12h','24h','36h','48h','60h','72h']}
def phase_sets(win):
    out={}
    for an in ('FL','TN'):
        x=u[u.analysis==an]
        genes=set(x.gene_id)
        d={}
        for p,st in win.items():
            d[p]=set(x[x.stage.isin(st)&x.robust_ase].gene_id)
        out[an]=(genes,d)
    return out
def upset(win):
    S=phase_sets(win)
    z=ov.copy()
    z=z[z.gene_africa.isin(S['FL'][0]) & z.gene_dura.isin(S['TN'][0])]
    cols=[]
    for an,key in (('FL','gene_africa'),('TN','gene_dura')):
        for p in win:
            c=f'{an}_{p}'; cols.append(c); z[c]=z[key].isin(S[an][1][p])
    z['pattern']=z[cols].astype(int).astype(str).agg(''.join,axis=1)
    pat=z.groupby('pattern').size().rename('orthogroups').reset_index().sort_values('orthogroups',ascending=False)
    ss=pd.DataFrame({'set':cols,'orthogroups':[int(z[c].sum()) for c in cols]})
    return pat,ss,z
if __name__=='__main__':
    pat,ss,z=upset(OLD)
    pub_p=pd.read_csv('../../../redraw/ED3/data/ED3a_upset_intersections.tsv',sep='\t',dtype={'pattern':str})
    pub_s=pd.read_csv('../../../redraw/ED3/data/ED3a_upset_set_sizes.tsv',sep='\t')
    print('rows in overlap used',len(z))
    m=ss.merge(pub_s,on='set',suffixes=('_re','_pub')); print(m)
    pp=pat.merge(pub_p,on='pattern',how='outer',suffixes=('_re','_pub')).fillna(0)
    pp['diff']=pp.orthogroups_re-pp.orthogroups_pub
    print('pattern diffs', (pp['diff']!=0).sum(), 'of', len(pp)); print(pp.sort_values('orthogroups_pub',ascending=False).head(8))
    pat2,ss2,z2=upset(NEW)
    pat2.to_csv('ED3a_upset_intersections_new.tsv',sep='\t',index=False); ss2.to_csv('ED3a_upset_set_sizes_new.tsv',sep='\t',index=False)
    print(ss2); print(pat2.head(18)); print('n intersections',len(pat2))
