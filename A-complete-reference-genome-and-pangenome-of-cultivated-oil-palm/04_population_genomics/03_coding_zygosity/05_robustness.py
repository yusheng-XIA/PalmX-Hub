"""Per-accession derived load (missense / synonymous / stop-gained), date-palm polarized, 308 accessions."""
import gzip, re, collections, glob, os, numpy as np
O=os.environ.get('WORK', 'coding_zygosity')
B=os.environ.get('DSV_DIR', 'dsv_analysis')
PAF=B+'/results_hap38/06_dSNP_phoenix_polarity/alignments/Phoenix_vs_Africa_hap2.asm20.cs.primary.paf'
# genome
fai={l.split()[0]:list(map(int,l.split()[1:4])) for l in open(B+'/input/Africa_hap2.fa.fai')}
fh=open(B+'/input/Africa_hap2.fa','rb')
def chrom_seq(c):
    L,off,lb=fai[c][0],fai[c][1],fai[c][2]; lB=lb+1
    fh.seek(off); raw=fh.read(L+L//lb+1).decode().replace('\n','')
    return raw[:L].upper()
code={}; bases='TCAG'; aa='FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'; i=0
for a in bases:
    for b in bases:
        for d in bases: code[a+b+d]=aa[i]; i+=1
comp={'A':'T','C':'G','G':'C','T':'A'}
# models
models=collections.defaultdict(list)
for l in open(O+'/cds_models.tsv'):
    m,c,s,ex=l.rstrip().split('\t'); ex=[tuple(map(int,e.split('-'))) for e in ex.split(',')]; models[c].append((m,s,ex))
# site -> (model,strand,cds_index)
samples=None; rows=[]
per=None
cats=['missense','synonymous','stop_gained']
for gf in sorted(glob.glob(O+'/geno/chr*.tsv.gz')):
    c=os.path.basename(gf).split('.')[0]
    seq=chrom_seq(c)
    pos2=collections.defaultdict(list)
    for m,s,ex in models[c]:
        k=0; tot=sum(b-a+1 for a,b in ex)
        coords=[p for a,b in ex for p in range(a,b+1)]
        if s=='-': coords=coords[::-1]
        for idx,p in enumerate(coords): pos2[p].append((s,coords,idx))
    sites=[]
    with gzip.open(gf,'rt') as g:
        hdr=next(g).rstrip('\n').split('\t')
        if samples is None:
            samples=[h.split(']')[1].split(':')[0] for h in hdr[4:]]
        for l in g:
            f=l.rstrip('\n').split('\t'); p=int(f[1]); ref,alt=f[2],f[3]
            if p not in pos2 or seq[p-1]!=ref: continue
            s,coords,idx=pos2[p][0]
            cs=idx-idx%3
            if cs+3>len(coords): continue
            cod=''.join(seq[q-1] for q in coords[cs:cs+3])
            if s=='-': cod=''.join(comp.get(x,'N') for x in cod); a2=comp[alt]
            else: a2=alt
            cod2=cod[:idx%3]+a2+cod[idx%3+1:]
            A1,A2=code.get(cod,'X'),code.get(cod2,'X')
            if 'X' in (A1,A2): continue
            cat='synonymous' if A1==A2 else ('stop_gained' if A2=='*' else ('stop_lost' if A1=='*' else 'missense'))
            gt=np.array([(-1 if '.' in x else int(x[0])+int(x[-1])) for x in f[4:]],dtype=np.int8)
            sites.append((c,p,ref,alt,cat,gt))
    rows.extend(sites); print(c,len(sites),flush=True)
# Phoenix base at sites
need=collections.defaultdict(dict)
for i,(c,p,ref,alt,cat,gt) in enumerate(rows): need[c][p]=i
phx={}
import bisect
SP={c:np.array(sorted(d)) for c,d in need.items()}
tok=re.compile(r'(:\d+|\*[a-z][a-z]|\+[a-z]+|-[a-z]+|~[a-z]{2}\d+[a-z]{2})')
for l in open(PAF):
    f=l.rstrip('\n').split('\t')
    if int(f[11])<20: continue
    t=f[5]
    if t not in need: continue
    ts=int(f[7]); te=int(f[8]); arr=SP[t]; nd=need[t]
    lo=np.searchsorted(arr,ts+1); hi=np.searchsorted(arr,te,side='right')
    if hi<=lo: continue
    cs=[x for x in f[12:] if x.startswith('cs:Z:')]
    if not cs: continue
    pos=ts
    for m in tok.findall(cs[0][5:]):
        if m[0]==':':
            n=int(m[1:]); a1=np.searchsorted(arr,pos+1); a2=np.searchsorted(arr,pos+n,side='right')
            for q in arr[a1:a2]: phx.setdefault(nd[int(q)],set()).add('R')
            pos+=n
        elif m[0]=='*':
            q=pos+1
            if q in nd: phx.setdefault(nd[q],set()).add(m[2].upper())
            pos+=1
        elif m[0]=='-':
            pos+=len(m)-1
print('phoenix-covered sites',len(phx),'of',len(rows),flush=True)
# polarize
grp={}
for l in open(os.environ.get('GROUPS', 'groups308.txt')):   # accession -> origin group
    a,b=l.split(); grp[a]=b.split(',')
S=len(samples)
acc={k:np.zeros((S,4)) for k in cats}   # derived alleles, hom-derived, het, called sites
site_out=gzip.open(O+'/polarized_sites.tsv.gz','wt'); site_out.write('Chrom\tPos\tRef\tAlt\tCat\tDerived\tDAF_all\n')
G={}
for i,(c,p,ref,alt,cat,gt) in enumerate(rows):
    if i not in phx or cat not in cats: continue
    st=phx[i]
    if len(st)!=1: continue
    b=next(iter(st)); b=ref if b=='R' else b
    if b==ref: d=gt.astype(float); der='ALT'
    elif b==alt: d=np.where(gt<0,-1,2-gt).astype(float); der='REF'
    else: continue
    ok=d>=0
    if ok.sum()<0.8*S: continue
    daf=d[ok].sum()/(2*ok.sum())
    if daf==0 or daf==1: pass
    site_out.write(f'{c}\t{p}\t{ref}\t{alt}\t{cat}\t{der}\t{daf:.4f}\n')
    A=acc[cat]; A[ok,0]+=d[ok]; A[ok,1]+=(d[ok]==2); A[ok,2]+=(d[ok]==1); A[ok,3]+=1
    G.setdefault(cat,[]).append((c,np.where(ok,d,np.nan)))
site_out.close()
# robustness: per chromosome x sample, all sites and DAF<=0.5 sites
import pickle
res={}
for cat in cats:
    ch=np.array([c for c,_ in G[cat]]); M=np.vstack([d for _,d in G[cat]])
    daf=np.nanmean(M,1)/2
    for lab,msk in [('all',np.ones(len(ch),bool)),('daf50',daf<=0.5)]:
        for c in sorted(set(ch)):
            m=msk&(ch==c); X=M[m]
            res[(cat,lab,c)]=(np.nansum(X,0),np.nansum(X==2,0),np.sum(~np.isnan(X),0))
pickle.dump((samples,res),open(O+'/robust.pkl','wb'))
print('robust saved',flush=True)
with open(O+'/per_sample_load.tsv','w') as o:
    o.write('Sample\tGroups\t'+'\t'.join(f'{k}_{x}' for k in cats for x in ['der','hom','het','n'])+'\n')
    for j,sm in enumerate(samples):
        o.write(sm+'\t'+','.join(grp.get(sm,['NA']))+'\t'+'\t'.join(f'{acc[k][j,x]:.0f}' for k in cats for x in range(4))+'\n')
# group Rxy (Do et al. 2015) with chromosome jackknife
names=['COMM','NONC','AFR','HHG','IDB','SAEG','SEAA','SEAB','K4P1','K4P2','K4P3','K4P4']
idx={g:np.array([j for j,sm in enumerate(samples) if g in grp.get(sm,[])]) for g in names}
def freqs(cat):
    ch=np.array([c for c,_ in G[cat]]); M=np.vstack([d for _,d in G[cat]])
    return ch,M
F={cat:freqs(cat) for cat in cats}
def L(cat,x,y,mask=None):
    ch,M=F[cat]; fx=np.nanmean(M[:,idx[x]],1)/2; fy=np.nanmean(M[:,idx[y]],1)/2
    m=np.isfinite(fx)&np.isfinite(fy) if mask is None else (np.isfinite(fx)&np.isfinite(fy)&(ch!=mask))
    return np.nansum((fx*(1-fy))[m]), np.nansum((fy*(1-fx))[m])
chs=sorted(set(F['synonymous'][0]))
with open(O+'/rxy.tsv','w') as o:
    o.write('X\tY\tcat\tRxy_norm_syn\tjk_se\n')
    for x,y in [('COMM','NONC'),('SEAA','NONC'),('SAEG','NONC'),('SEAA','AFR'),('SAEG','AFR'),('HHG','AFR'),('IDB','AFR'),('SEAB','AFR'),('K4P2','K4P3'),('K4P1','K4P3'),('K4P4','K4P3')]:
        for cat in ['missense','stop_gained']:
            def r(mask=None):
                a,b=L(cat,x,y,mask); s1,s2=L('synonymous',x,y,mask); return (a/b)/(s1/s2)
            full=r(); jk=np.array([r(c) for c in chs]); n=len(jk); se=np.sqrt((n-1)/n*np.sum((jk-jk.mean())**2))
            o.write(f'{x}\t{y}\t{cat}\t{full:.4f}\t{se:.4f}\n')
print('done',flush=True)
