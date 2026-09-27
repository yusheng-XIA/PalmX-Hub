import json,glob,os,collections,csv
sites=json.load(open('prom/fixed_sites.json'))
tr=str.maketrans('ACGT','TGCA')
def rc(s): return s.translate(tr)[::-1]
K=8
eoK={};egK={}
for idx,s in enumerate(sites):
    if s['region']=='up' and s['col']< (len(s['eo'][0])+0) : pass
    l,b,r=s['eo']; k=l[-K:]+b+r[:K]
    if len(k)!=2*K+1: continue
    eoK[k]=idx; eoK[rc(k)]=idx
    for l2,b2,r2 in s['eg']:
        k2=l2[-K:]+b2+r2[:K]
        if len(k2)==2*K+1: egK[k2]=idx; egK[rc(k2)]=idx
# remove kmers present in both
both=set(eoK)&set(egK)
for k in both: eoK.pop(k,None); egK.pop(k,None)
L=2*K+1
out=[]
for f in sorted(glob.glob('dl/*.hits.tsv')):
    r=os.path.basename(f).split('.')[0]
    if not os.path.exists(f'dl/{r}.done'): continue
    samf=f'rq/{r}.sam'
    if not os.path.exists(samf): continue
    keep=set()
    for l in open(samf):
        q=l.split('\t',3)
        if q[2] in ('African_hap2__evm.TU.chr11B.1497','American_hap1__evm.TU.chr11A.1059'): keep.add(q[0])
    eo=eg=0; site_eo=collections.Counter(); site_eg=collections.Counter()
    for line in open(f):
        p=line.rstrip('\n').split('\t')
        if len(p)<2: continue
        if p[0] not in keep: continue
        s=p[1]; he=set(); hg=set()
        for i in range(len(s)-L+1):
            k=s[i:i+L]
            if k in eoK: he.add(eoK[k])
            elif k in egK: hg.add(egK[k])
        if he and not hg: eo+=1; site_eo.update(he)
        elif hg and not he: eg+=1; site_eg.update(hg)
    out.append((r,eo,eg))
w=open('allele_kmer_counts.tsv','w'); w.write('run\tOLE16a_Eo_allele_reads\tOLE16a_Eg_allele_reads\n')
for x in out: w.write('\t'.join(map(str,x))+'\n')
print(len(eoK),len(egK))
for x in out: print(*x)
