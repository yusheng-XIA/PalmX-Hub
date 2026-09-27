import sys, collections
hdr=open(sys.argv[1]).read().rstrip('\n').split('\t')
samples=hdr[9:]
grp={}
for p in (1,2,3,4):
    for s in open(f'{sys.argv[3]}/K4_Pop{p}.txt'):
        s=s.strip()
        if s: grp[s]=p
print('samples in vcf',len(samples),'assigned',sum(s in grp for s in samples))
hap=collections.Counter(); sites=collections.Counter(); site_ok={p:0 for p in (1,2,3,4)}; n=0
gtcount=collections.Counter()
for line in open(sys.argv[2]):
    if line.startswith('#'): continue
    f=line.rstrip('\n').split('\t')
    n+=1
    bad=set()
    for s,g in zip(samples,f[9:]):
        gt=g.split(':')[0]
        gtcount[gt]+=1
        if '/' not in gt and '|' not in gt:
            hap[s]+=1; bad.add(grp.get(s))
    for p in (1,2,3,4):
        if p not in bad: site_ok[p]+=1
print('lines',n,'chrom range',f[0],f[1])
print('GT types',gtcount.most_common(12))
print('fully-diploid sites per pop',site_ok)
print('samples with haploid GT (sample,pop,count):',[(s,grp.get(s),c) for s,c in hap.most_common()])
