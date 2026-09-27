"""G6 checks: per-site table for all classified CDS SNPs (308 accessions) with
 - own codon class, CDS index / length, exon index / n exons (last-exon flag)
 - Phoenix (date palm) state (asm20 PAF, MAPQ>=20; as the formal dSNP rule)
 - E. oleifera states (MZ4 hap1/hap2 SyRI vs FL-Hap2, enh_B pickles)
 - Nigerian (EG_niriliya) alignment coverage (paftools call R regions)
 - SnpEff first annotation from the SnpEff-annotated population VCFs
Genotype matrix saved as int8 (-1 missing) in npz."""
import gzip, re, collections, glob, os, pickle, subprocess, numpy as np
O=os.environ.get('WORK', 'coding_zygosity')
B=os.environ.get('DSV_DIR', 'dsv_analysis')
PAF=B+'/results_hap38/06_dSNP_phoenix_polarity/alignments/Phoenix_vs_Africa_hap2.asm20.cs.primary.paf'
MZ=os.environ.get('OLEIFERA_DIR', 'work/outgroup_sensitivity/mz4')   # 12_*/06_outgroup_sensitivity/01a output
NR=os.environ.get('NIGERIAN_VAR', 'snv_calling/minimap2_results/EG_niriliya/EG_niriliya_vs_Africa_hap2.var.txt')
SE=os.environ.get('SNPEFF_DIR', 'snpeff')
BCF='bcftools'
fai={l.split()[0]:list(map(int,l.split()[1:4])) for l in open(B+'/input/Africa_hap2.fa.fai')}
fh=open(B+'/input/Africa_hap2.fa','rb')
def chrom_seq(c):
    L,off,lb=fai[c]; fh.seek(off); raw=fh.read(L+L//lb+1).decode().replace('\n',''); return raw[:L].upper()
code={}; bases='TCAG'; aa='FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'; i=0
for a in bases:
    for b in bases:
        for d in bases: code[a+b+d]=aa[i]; i+=1
comp={'A':'T','C':'G','G':'C','T':'A'}
models=collections.defaultdict(list)
for l in open(O+'/cds_models.tsv'):
    m,c,s,ex=l.rstrip().split('\t'); ex=[tuple(map(int,e.split('-'))) for e in ex.split(',')]; models[c].append((m,s,ex))
samples=None; meta=[]; GT=[]
for gf in sorted(glob.glob(O+'/geno/chr*.tsv.gz')):
    c=os.path.basename(gf).split('.')[0]; seq=chrom_seq(c)
    pos2={}
    for m,s,ex in models[c]:
        exo=ex if s=='+' else ex[::-1]; coords=[]; exi=[]
        for k,(a,b) in enumerate(exo):
            r=range(a,b+1) if s=='+' else range(b,a-1,-1)
            coords.extend(r); exi.extend([k]*(b-a+1))
        for idx,p in enumerate(coords):
            if p not in pos2: pos2[p]=(m,s,coords,idx,exi[idx],len(exo))
    with gzip.open(gf,'rt') as g:
        hdr=next(g).rstrip('\n').split('\t')
        if samples is None: samples=[h.split(']')[1].split(':')[0] for h in hdr[4:]]
        for l in g:
            f=l.rstrip('\n').split('\t'); p=int(f[1]); ref,alt=f[2],f[3]
            if p not in pos2 or seq[p-1]!=ref: continue
            m,s,coords,idx,ei,ne=pos2[p]; cs=idx-idx%3
            if cs+3>len(coords): continue
            cod=''.join(seq[q-1] for q in coords[cs:cs+3])
            if s=='-': cod=''.join(comp.get(x,'N') for x in cod); a2=comp[alt]
            else: a2=alt
            cod2=cod[:idx%3]+a2+cod[idx%3+1:]
            A1,A2=code.get(cod,'X'),code.get(cod2,'X')
            if 'X' in (A1,A2): continue
            cat='synonymous' if A1==A2 else ('stop_gained' if A2=='*' else ('stop_lost' if A1=='*' else 'missense'))
            meta.append([c,p,ref,alt,cat,m,s,idx,len(coords),ei,ne])
            GT.append(np.array([(-1 if '.' in x else int(x[0])+int(x[-1])) for x in f[4:]],dtype=np.int8))
    print(c,len(meta),flush=True)
GT=np.vstack(GT)
need=collections.defaultdict(dict)
dup=collections.Counter((r[0],r[1]) for r in meta)
for i,r in enumerate(meta): need[r[0]][r[1]]=i
print('duplicated positions',sum(v>1 for v in dup.values()),flush=True)
SP={c:np.array(sorted(d)) for c,d in need.items()}
# Phoenix
phx=collections.defaultdict(set)
tok=re.compile(r'(:\d+|\*[a-z][a-z]|\+[a-z]+|-[a-z]+|~[a-z]{2}\d+[a-z]{2})')
for l in open(PAF):
    f=l.rstrip('\n').split('\t')
    if int(f[11])<20 or f[5] not in need: continue
    t=f[5]; ts=int(f[7]); te=int(f[8]); arr=SP[t]; nd=need[t]
    if np.searchsorted(arr,te,side='right')<=np.searchsorted(arr,ts+1): continue
    cs=[x for x in f[12:] if x.startswith('cs:Z:')]
    if not cs: continue
    pos=ts
    for m in tok.findall(cs[0][5:]):
        if m[0]==':':
            n=int(m[1:]); a1=np.searchsorted(arr,pos+1); a2=np.searchsorted(arr,pos+n,side='right')
            for q in arr[a1:a2]: phx[nd[int(q)]].add('R')
            pos+=n
        elif m[0]=='*':
            q=pos+1
            if q in nd: phx[nd[q]].add(m[2].upper())
            pos+=1
        elif m[0]=='-': pos+=len(m)-1
def state(st,ref,alt):
    if len(st)!=1: return 'NA'
    b=next(iter(st)); b=ref if b=='R' else b
    return 'REF' if b==ref else ('ALT' if b==alt else 'NA')
# MZ4
mz={}
for h in ['meizhou4_hap1','meizhou4_hap2']:
    d=pickle.load(open(f'{MZ}/{h}.pkl','rb')); out={}
    for c,arr in SP.items():
        bl=d['blocks'].get(c,np.zeros((0,2),np.int64)); dl=d['del'].get(c,np.zeros((0,2),np.int64)); sn=d['snp'].get(c,{})
        j=np.searchsorted(bl[:,0],arr-1,side='right')-1
        inb=(j>=0)&(arr-1<bl[np.clip(j,0,None),1]) if len(bl) else np.zeros(len(arr),bool)
        k=np.searchsorted(dl[:,0],arr,side='right')-1
        ind=(k>=0)&(arr<=dl[np.clip(k,0,None),1]) if len(dl) else np.zeros(len(arr),bool)
        for q,a,b in zip(arr,inb,ind):
            i=need[c][int(q)]; r=meta[i]
            if not a or b: out[i]='NA'
            elif int(q) in sn: out[i]='ALT' if sn[int(q)]==r[3] else 'NA'
            else: out[i]='REF'
    mz[h]=out
    print(h,'done',flush=True)
# Nigerian aligned regions
nreg=collections.defaultdict(list)
for l in open(NR):
    if l[0]=='R':
        f=l.split('\t'); nreg[f[1]+'B' if not f[1].endswith('B') else f[1]].append((int(f[2]),int(f[3])))
ncov={}
for c,arr in SP.items():
    iv=np.array(sorted(nreg.get(c,[])),dtype=np.int64).reshape(-1,2)
    j=np.searchsorted(iv[:,0],arr-1,side='right')-1
    ok=(j>=0)&(arr-1<iv[np.clip(j,0,None),1]) if len(iv) else np.zeros(len(arr),bool)
    for q,v in zip(arr,ok): ncov[need[c][int(q)]]=bool(v)
# SnpEff first annotation (03_snpeff_first_annotation.sh; chr01B, 08B, 10B, 14B, 16B)
se={}
G5=O+'/snpeff'   # 03_snpeff_first_annotation.sh output
for c in ['chr01B','chr08B','chr10B','chr14B','chr16B']:
    nd=need[c]
    with gzip.open(f'{G5}/{c}.ann.tsv.gz','rt') as g:
        for l in g:
            f=l.rstrip('\n').split('\t'); i=nd.get(int(f[0]))
            if i is None or meta[i][2]!=f[1] or meta[i][3]!=f[2]: continue
            se[i]=(f[3],f[5])
print('snpeff matched',len(se),flush=True)
with gzip.open(O+'/g6_sites_full.tsv.gz','wt') as o:
    o.write('Chrom\tPos\tRef\tAlt\tCat\tModel\tStrand\tCDS_idx\tCDS_len\tExon_idx\tN_exons\tPhoenix\tMZ4_h1\tMZ4_h2\tNigerian_aligned\tSnpEff_effect\tSnpEff_gene\n')
    for i,r in enumerate(meta):
        o.write('\t'.join(map(str,r))+f"\t{state(phx.get(i,set()),r[2],r[3])}\t{mz['meizhou4_hap1'].get(i,'NA')}\t{mz['meizhou4_hap2'].get(i,'NA')}\t{ncov.get(i,False)}\t{se.get(i,('NA','NA'))[0]}\t{se.get(i,('NA','NA'))[1]}\n")
np.savez_compressed(O+'/g6_gt.npz',GT=GT,samples=np.array(samples))
print('done',flush=True)
