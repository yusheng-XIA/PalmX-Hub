import sys,collections
B='${ANALYSIS_DIR}/20_results/10_database/05_cds_sequences/'
A='${ANALYSIS_DIR}/20_results/Figure1/10_new_figure/06_DEA/DEA_new/allele/'
def fa(f):
    d={};k=None;s=[]
    for l in open(f):
        l=l.rstrip()
        if l.startswith('>'):
            if k: d[k]=''.join(s).upper()
            k=l[1:].split()[0];s=[]
        else: s.append(l)
    if k: d[k]=''.join(s).upper()
    return d
def norm(d):
    o=dict(d)
    for k,v in d.items():
        o.setdefault(k.replace('evm.model.','').replace('evm.TU.',''),v)
        o.setdefault(k.replace('evm.TU.','evm.model.'),v)
    return o
sets={'bk':('bk_hap1_cds.fasta','bk_hap2_cds.fasta'),'American_Africa':('American_hap1_cds.fasta','Africa_hap2_cds.fasta')}
for m,(f1,f2) in sets.items():
    h1=norm(fa(B+f1));h2=norm(fa(B+f2))
    print(m,'ids',list(h1)[:2],list(h2)[:2])
    for fn in (('same_cds_pairs.tsv.bak','biallelic_pairs.tsv.bak') if m=='bk' else ('same_cds_pairs.tsv','biallelic_pairs.tsv')):
        c=collections.Counter()
        for i,l in enumerate(open(A+m+'/final_results/'+fn)):
            if i==0: continue
            a,b=l.split('\t')[:2]
            s1=h1.get(a);s2=h2.get(b)
            if s1 is None or s2 is None: c['missing']+=1;continue
            if s1==s2: c['identical']+=1
            elif len(s1)==len(s2): c['same_len_diff']+=1
            else: c['diff_len']+=1
        print(m,fn,dict(c))
