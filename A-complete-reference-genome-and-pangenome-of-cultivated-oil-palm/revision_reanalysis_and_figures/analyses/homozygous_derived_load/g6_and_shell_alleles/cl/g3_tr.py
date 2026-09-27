import re
A='${ANALYSIS_DIR}/08_hifi_chromosome/11_final_corrected_39genomes_annotations_20260812/attempt_20260812_01'
O='${CLUSTER_WORK}/headline_pop/out'
fai={l.split()[0]:l.split() for l in open(A+'/01_genomes/African_hap2.fa.fai')}
print(list(fai)[:20])
c=[k for k in fai if k in('chr01','chr01B')][0]
_,L,off,lb,lB=fai[c]; L,off,lb,lB=map(int,(L,off,lb,lB))
fh=open(A+'/01_genomes/African_hap2.fa','rb')
def fetch(s,e):  # 1-based inclusive
    s0=s-1; a=off+(s0//lb)*lB+s0%lb; b=off+((e-1)//lb)*lB+(e-1)%lb
    fh.seek(a); return fh.read(b-a+1).decode().replace('\n','').upper()
cds={}
for l in open(O+'/g3_region.gff3'):
    f=l.rstrip('\n').split('\t')
    if f[2]=='CDS':
        p=re.search('Parent=([^;]+)',f[8]).group(1); cds.setdefault(p,[]).append((int(f[3]),int(f[4]),f[6]))
code={}
bases='TCAG';aa='FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
i=0
for a in bases:
    for b in bases:
        for d in bases: code[a+b+d]=aa[i]; i+=1
comp=str.maketrans('ACGT','TGCA')
for p,ex in sorted(cds.items(), key=lambda x:min(e[0] for e in x[1])):
    st=ex[0][2]; ex=sorted(ex)
    s=''.join(fetch(a,b) for a,b,_ in ex)
    if st=='-': s=s.translate(comp)[::-1]
    pr=''.join(code.get(s[i:i+3],'X') for i in range(0,len(s)-2,3))
    tag='MADS' if re.search('GR[GK]K[IV]E[IL]K',pr) else ''
    print(p, st, min(e[0] for e in ex), max(e[1] for e in ex), len(pr), tag, pr[:70])
    if tag: print('  exons', ex)
print('lead SNP base', fetch(3153030,3153030), 'context', fetch(3153020,3153040))
