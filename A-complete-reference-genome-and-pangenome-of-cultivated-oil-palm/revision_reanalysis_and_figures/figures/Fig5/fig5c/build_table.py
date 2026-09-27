import pandas as pd
from collections import defaultdict
sv=pd.read_csv('sv.bed',sep='\t',header=None,names=['Chrom','ctx_s0','ctx_e','SV_ID','SVTYPE','SVLEN_bp','Start','End','key'])
def load(f):
    d=defaultdict(set)
    for l in open(f):
        a,b=l.rstrip('\n').split('\t'); d[a].add(b)
    return d
H={k:load(f'hit_{k}.txt') for k in ['cds','gene','up','down']}
near=defaultdict(lambda:(set(),None))
for l in open('nearest.txt'):
    a,g,d=l.rstrip('\n').split('\t')
    if g=='.': continue
    s,_=near[a]; s.add(g); near[a]=(s,int(d))
ctx=[];gid=[];ndist=[]
for i in sv.key:
    for k,lab in [('cds','CDS'),('gene','Intron'),('up','Upstream 2 kb'),('down','Downstream 2 kb')]:
        if H[k].get(i):
            ctx.append(lab); gid.append(','.join(sorted(H[k][i]))); ndist.append(0 if k in('cds','gene') else near[i][1]); break
    else:
        ctx.append('Intergenic'); s,d=near[i]; gid.append(','.join(sorted(s))); ndist.append(d)
sv['Genomic_context']=ctx; sv['Gene_ID']=gid; sv['Distance_to_gene_bp']=ndist
sv['Context_interval_start0']=sv.ctx_s0; sv['Context_interval_end']=sv.ctx_e
chrorder=sorted(sv.Chrom.unique())
sv['row']=sv.key.str[1:].astype(int); sv=sv.sort_values('row')
out=sv[['SV_ID','Chrom','Start','End','SVTYPE','SVLEN_bp','Context_interval_start0','Context_interval_end','Genomic_context','Gene_ID','Distance_to_gene_bp']]
out.to_csv('Fig5c_SV_gene_context_loci.tsv',sep='\t',index=False)
order=['CDS','Intron','Upstream 2 kb','Downstream 2 kb','Intergenic']
print(out.Genomic_context.value_counts().reindex(order))
ct=pd.crosstab(out.SVTYPE,out.Genomic_context,margins=True,margins_name='Total')[order+['Total']]
print(ct); ct.to_csv('Fig5c_context_by_SVTYPE.tsv',sep='\t')
print('rows',len(out),'unique SV_ID',out.SV_ID.nunique())
dup=out[out.SV_ID.duplicated(keep=False)]
print('dup rows',len(dup)); print(dup.head(6).to_string())
print('identical-coordinate dup rows', out.duplicated(['SV_ID','Chrom','Start','End','SVTYPE'],keep=False).sum())
print('empty gene id by ctx:'); print(out[out.Gene_ID==''].Genomic_context.value_counts())
print('multi-gene rows by ctx:'); print(out[out.Gene_ID.str.contains(',')].Genomic_context.value_counts())
