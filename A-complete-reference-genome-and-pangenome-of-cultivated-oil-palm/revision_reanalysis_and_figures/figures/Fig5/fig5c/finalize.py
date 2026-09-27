import pandas as pd
d=pd.read_csv('Fig5c_SV_gene_context_loci.tsv',sep='\t',keep_default_na=False)
if 'Catalogue_row' in d.columns: raise SystemExit('already finalized')
g=pd.read_csv('gene_id.bed',sep='\t',header=None,names=['c','s','e','id','x','st']).set_index('id')
def gap(r):
    if r.Genomic_context in ('CDS','Intron'): return 0
    if r.Genomic_context=='Intergenic': return int(r.Distance_to_gene_bp)-1   # bedtools closest -d: book-ended = 1
    best=None
    for gid in r.Gene_ID.split(','):
        s,e=g.at[gid,'s'],g.at[gid,'e']
        x=max(s-r.Context_interval_end, r.Context_interval_start0-e, 0)
        best=x if best is None else min(best,x)
    return best
d['Distance_to_gene_bp']=d.apply(gap,axis=1)
d.insert(0,'Catalogue_row',range(1,len(d)+1))
d=d.rename(columns={'Chrom':'Chrom','Start':'Start','End':'End'})
d.to_csv('Fig5c_SV_gene_context_loci.tsv',sep='\t',index=False)
print(d.groupby('Genomic_context').Distance_to_gene_bp.agg(['min','max']))
print(d.Genomic_context.value_counts()); print(d.shape)
