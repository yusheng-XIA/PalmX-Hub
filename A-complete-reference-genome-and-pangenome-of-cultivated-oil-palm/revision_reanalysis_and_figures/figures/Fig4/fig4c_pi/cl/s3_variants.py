import glob
POPS=['K4_Pop1','K4_Pop2','K4_Pop3','K4_Pop4']
PAIRS=['K4_Pop1_K4_Pop2','K4_Pop1_K4_Pop3','K4_Pop1_K4_Pop4','K4_Pop2_K4_Pop3','K4_Pop2_K4_Pop4','K4_Pop3_K4_Pop4']
def rd(p):
    out=[]
    with open(p) as h:
        hd=h.readline().split()
        for l in h: out.append(dict(zip(hd,l.split())))
    return out
def fx(nm,suf):
    r=[]
    for p in sorted(glob.glob('out/fixed/*__'+nm+suf)): r+=rd(p)
    return r
def bad(r):  # regions where a block of samples is absent from the joint-genotyping output
    s=int(r['BIN_END'])
    return (r['CHROM']=='chr04B' and s>82400000) or (r['CHROM']=='chr05B' and s>19000000)
def pi(rows): 
    v=[(float(r['PI']),int(r['BIN_END'])-int(r['BIN_START'])+1) for r in rows]; return sum(a*b for a,b in v)/sum(b for a,b in v), len(v)
def fst(rows):
    v=[(float(r['MEAN_FST']),int(r['N_VARIANTS'])) for r in rows if r['MEAN_FST'] not in ('nan','-nan')]; return sum(a*b for a,b in v)/sum(b for a,b in v), len(v)
orig={'K4_Pop1':0.0056650237166,'K4_Pop2':0.00470491658359,'K4_Pop3':0.00577306769722,'K4_Pop4':0.00597151287293,
'K4_Pop1_K4_Pop2':0.126357843235,'K4_Pop1_K4_Pop3':0.0657473265601,'K4_Pop1_K4_Pop4':0.0436513270168,'K4_Pop2_K4_Pop3':0.0928941177597,'K4_Pop2_K4_Pop4':0.129326048602,'K4_Pop3_K4_Pop4':0.0819804232249}
print('metric\tgroup\tfigure(orig)\tfixed_all16\tfixed_all16_%\tfixed_excl_chr04B>82.4Mb_chr05B>19Mb\texcl_%\tnwin_excl')
for nm in POPS:
    r=fx(nm,'_100kb.pi.windowed.pi'); a,_=pi(r); b,n=pi([x for x in r if not bad(x)])
    print('pi\t%s\t%.4f\t%.4f\t%+.1f\t%.4f\t%+.1f\t%d'%(nm,1e3*orig[nm],1e3*a,100*(a/orig[nm]-1),1e3*b,100*(b/orig[nm]-1),n))
for nm in PAIRS:
    r=fx(nm,'_100kb_fst.windowed.weir.fst'); a,_=fst(r); b,n=fst([x for x in r if not bad(x)])
    print('fst\t%s\t%.4f\t%.4f\t%+.1f\t%.4f\t%+.1f\t%d'%(nm,orig[nm],a,100*(a/orig[nm]-1),b,100*(b/orig[nm]-1),n))
