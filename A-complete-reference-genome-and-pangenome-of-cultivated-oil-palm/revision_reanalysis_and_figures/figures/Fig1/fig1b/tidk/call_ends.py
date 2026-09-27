import pandas as pd, itertools
def call(g, nwin=5, thr=3, strand=False):
    d=pd.read_csv(f'{g}_telomeric_repeat_windows.tsv',sep='\t')
    out=[]
    for i,(cid,s) in enumerate(d.groupby('id',sort=False),1):
        s=s.sort_values('window')
        L=s.head(nwin); R=s.tail(nwin)
        if strand: lv=L.reverse_repeat_number; rv=R.forward_repeat_number
        else: lv=L[['forward_repeat_number','reverse_repeat_number']].max(axis=1); rv=R[['forward_repeat_number','reverse_repeat_number']].max(axis=1)
        out.append((i,cid,int(lv.max()),int(lv.max()>=thr),int(rv.max()),int(rv.max()>=thr)))
    return pd.DataFrame(out,columns=['chr','id','L_max','L_pos','R_max','R_pos'])
if __name__=='__main__':
    for g in ['EG11','EO12']:
        t=call(g); t.insert(0,'assembly',g)
        # add strand info of the max
        print(t.to_string(index=False)); print(g,'positive ends =',t.L_pos.sum()+t.R_pos.sum(),'/32\n')
        t.to_csv(f'{g}_ends_primary_50kb_ge3.tsv',sep='\t',index=False)
    print('sensitivity (EG11, EO12):')
    for strand in [False,True]:
        for nwin in [1,2,5,10]:
            row=[]
            for thr in [3,10,50,100]:
                r=[ (lambda t:t.L_pos.sum()+t.R_pos.sum())(call(g,nwin,thr,strand)) for g in ['EG11','EO12']]
                row.append(f"thr>={thr}: {r[0]}/{r[1]}")
            print(f"strand-aware={strand} terminal={nwin*10}kb  ", ' | '.join(row))
