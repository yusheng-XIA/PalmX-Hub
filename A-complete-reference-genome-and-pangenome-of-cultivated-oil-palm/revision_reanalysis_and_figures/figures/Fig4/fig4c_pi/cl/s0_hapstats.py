import glob, collections
hdr=open('probe/header_chrom.txt').read().rstrip('\n').split('\t'); S=hdr[9:]
grp={}
for p in (1,2,3,4):
    for s in open('groups/K4_Pop%d.txt'%p):
        if s.strip(): grp[s.strip()]=p
sites=collections.Counter(); hap=collections.defaultdict(collections.Counter)
for f in glob.glob('stats/stats*.tsv'):
    for l in open(f):
        x=l.rstrip('\n').split('\t')
        if x[0]=='SITES': sites[x[1]]+=int(x[2])
        elif x[0]=='HAP': hap[x[1]][int(x[2])]+=int(x[3])
print('total sites',sum(sites.values()))
print('chrom\tsites\ttotal_dotGT\tsamples_dot>=50%sites(sample:pop)')
for c in sites:
    tot=sum(hap[c].values())
    big=['%s:P%d(%.0f%%)'%(S[i],grp[S[i]],100*n/sites[c]) for i,n in hap[c].most_common() if n>=0.5*sites[c]]
    print('%s\t%d\t%d\t%s'%(c,sites[c],tot,' '.join(big)))
