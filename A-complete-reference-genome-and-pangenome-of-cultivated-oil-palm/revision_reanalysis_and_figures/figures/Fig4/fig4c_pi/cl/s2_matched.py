import glob, collections
POPS=['K4_Pop1','K4_Pop2','K4_Pop3','K4_Pop4']
def rd(p, col):
    d={}
    with open(p) as h:
        hd=h.readline().split()
        for l in h:
            f=dict(zip(hd,l.split())); d[(f['CHROM'],int(f['BIN_START']))]=(float(f[col]),int(f['N_VARIANTS']))
    return d
def fixed(nm,suf,col):
    d={}
    for p in glob.glob('out/fixed/*__'+nm+suf): d.update(rd(p,col))
    return d
print('A) window-matched (same windows present in both):')
for nm in POPS:
    o=rd('orig/%s_100kb.pi.windowed.pi'%nm,'PI'); x=fixed(nm,'_100kb.pi.windowed.pi','PI')
    k=[w for w in o if w in x]
    print(' pi %s n=%d orig=%.6g fixed=%.6g  delta=%+.2f%%  nvar orig=%d fixed=%d'%(nm,len(k),sum(o[w][0] for w in k)/len(k),sum(x[w][0] for w in k)/len(k),100*(sum(x[w][0] for w in k)/sum(o[w][0] for w in k)-1),sum(o[w][1] for w in k),sum(x[w][1] for w in k)))
print('B) per-chromosome pi (fixed), windows mean:')
for nm in POPS:
    x=fixed(nm,'_100kb.pi.windowed.pi','PI'); by=collections.defaultdict(list)
    for (c,s),v in x.items(): by[c].append(v[0])
    print(' '+nm+' '+' '.join('%s:%.2f'%(c,1e3*sum(v)/len(v)) for c,v in sorted(by.items())))
print('C) chr04B pos of first window where orig ends / fixed Pop1 pi before/after 82.5 Mb:')
x=fixed('K4_Pop1','_100kb.pi.windowed.pi','PI')
a=[v[0] for (c,s),v in x.items() if c=='chr04B' and s<82500000]; b=[v[0] for (c,s),v in x.items() if c=='chr04B' and s>=82500000]
print(' chr04B <82.5Mb %.2f (%d win), >=82.5Mb %.2f (%d win)'%(1e3*sum(a)/len(a),len(a),1e3*sum(b)/len(b),len(b)))
a=[v[0] for (c,s),v in x.items() if c=='chr05B' and s<19000000]; b=[v[0] for (c,s),v in x.items() if c=='chr05B' and s>=19000000]
print(' chr05B <19Mb %.2f (%d win), >=19Mb %.2f (%d win)'%(1e3*sum(a)/len(a),len(a),1e3*sum(b)/len(b),len(b)))
