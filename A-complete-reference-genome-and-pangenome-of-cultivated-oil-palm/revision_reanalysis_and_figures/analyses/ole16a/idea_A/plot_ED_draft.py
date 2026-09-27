import pandas as pd, numpy as np, matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from Bio import Phylo
plt.rcParams.update({'font.size':6.5,'font.family':'Arial','axes.linewidth':0.5})
ST=["0d","15d","35d","50d","65d","80d","95d","110d","125d","140d","155d","170d","185d","12h","24h","36h","48h","60h","72h"]
CF,CT='#c0392b','#2c6fbb'
fig=plt.figure(figsize=(7.2,6.2)); gs=fig.add_gridspec(2,3,hspace=0.45,wspace=0.4)
# a tree (pruned labels)
ax=fig.add_subplot(gs[:,0]); t=Phylo.read('res/dom/oleosin_iq.contree','newick'); t.root_at_midpoint(); t.ladderize()
def lab(c):
    n=c.name or ''
    if n.startswith('REF'): p=n.split('|'); return f"{p[2]} [{p[1]}]"
    if n: p=n.split('|'); return f"{p[2].replace('evm.TU.','')} ({p[1]})"
    return None
for c in t.get_nonterminals():
    c.name=None
    if c.confidence is not None and c.confidence<70: c.confidence=None
Phylo.draw(t,axes=ax,do_show=False,label_func=lab,label_colors=lambda s: CF if ('chr11' in s or 'chr9.' in s) else ('#7a4' if 'chr04' in s or 'chr3.2082' in s else 'k'))
for tx in ax.texts: tx.set_fontsize(3.6)
ax.set_title('a  Oleosin phylogeny (IQ-TREE, UFBoot ≥70)',loc='left',fontsize=7); ax.set_xlabel(''); ax.set_ylabel(''); ax.set_yticks([]); ax.tick_params(labelsize=5)
# b RNA
L=pd.read_csv('res/rna/ld_gene_sample_long.tsv.gz',sep='\t')
ax=fig.add_subplot(gs[0,1]); x=np.arange(len(ST))
for g,ls,nm in [('evm.TU.chr11B.1497','-','OLE16a (chr11)'),('evm.TU.chr04B.697','--','OLE16b (chr04)'),('evm.TU.chr10B.1088',':','LDAP (chr10B.1088)')]:
    for m,c in [('FL',CF),('TN',CT)]:
        d=L[(L.gene==g)&(L.genotype==m)].groupby('stage').norm.mean().reindex(ST)
        ax.plot(x,np.log10(d+1),ls=ls,color=c,lw=0.9,label=f'{nm} {m}')
ax.axvspan(11.5,18.5,color='0.92',zorder=0); ax.set_xticks(x); ax.set_xticklabels(ST,rotation=90,fontsize=5)
ax.set_ylabel('log10(normalized count + 1)'); ax.legend(fontsize=4.2,frameon=False,ncol=1); ax.set_title('b  Bulk RNA (n = 3 per stage)',loc='left',fontsize=7)
# c precursor
S=pd.read_csv('res/ole16_group_sample_sums.tsv',sep='\t'); ax=fig.add_subplot(gs[0,2]); P=ST[10:]
for grp,ls,nm in [('OLE16a_chr11_shared','-','OLE16a shared peptides'),('OLE16b_chr04_shared','--','OLE16b shared peptide')]:
    for m,c in [('FL',CF),('TN',CT)]:
        d=S[(S.grp==grp)&(S.mat==m)]
        for s in P:
            v=d[d.stage==s]['Precursor.Quantity']; xi=P.index(s)
            ax.scatter([xi]*len(v),np.log10(v),s=4,color=c,marker='o' if ls=='-' else '^',lw=0)
        mm=d.groupby('stage')['Precursor.Quantity'].mean().reindex(P); ax.plot(range(len(P)),np.log10(mm),ls=ls,color=c,lw=0.8,label=f'{nm} {m}')
Q=pd.read_csv('res/run_precursor_quantiles.tsv',sep='\t',index_col=0); ax.axhline(np.log10(Q.p50.median()),color='0.5',lw=0.5,ls=':'); ax.text(0,np.log10(Q.p50.median())+0.1,'run median precursor',fontsize=4.5,color='0.4')
ax.set_xticks(range(len(P))); ax.set_xticklabels(P,rotation=90,fontsize=5); ax.set_ylabel('log10 precursor quantity (sum)'); ax.legend(fontsize=4.2,frameon=False)
ax.set_title('c  Peptide-level protein (DIA-NN)',loc='left',fontsize=7)
# d family deployment
D=pd.read_csv('res/prot/ld_family_stage_means_1e6.tsv',sep='\t'); ax=fig.add_subplot(gs[1,1]); P2=["185d","12h","24h","36h","48h","60h","72h"]
fams=[('oleosin','Oleosin'),('REF','LDAP (REF/SRPP)'),('caleosin','Caleosin')]
w=0.13
for i,(f,nm) in enumerate(fams):
    for j,(m,c) in enumerate([('FL',CF),('TN',CT)]):
        v=[D[D.family==f][f'{m}{s}'].sum() for s in P2]
        ax.bar(np.arange(len(P2))+(i*2+j-2.5)*w,np.log10(np.array(v)*1e6),w,color=c,alpha=[1,0.6,0.3][i],label=f'{nm} {m}')
ax.set_ylim(6,10); ax.set_xticks(range(len(P2))); ax.set_xticklabels(P2,fontsize=5); ax.set_ylabel('log10 summed directLFQ'); ax.legend(fontsize=4,frameon=False,ncol=2)
ax.set_title('d  LD-coat protein families, 185 d–72 h',loc='left',fontsize=7)
# e snRNA
R=pd.read_csv('res/sn/ld_genes_by_cluster.tsv',sep='\t'); R=R[R.gene=='evm.TU.chr11B.1497']; ax=fig.add_subplot(gs[1,2])
cls=[7,9,18,12,15,0,1,6,10,14,16]
for j,(lb,c) in enumerate([('FL_185',CF),('TN_185',CT)]):
    v=[R[(R.lib==lb)&(R.cl==k)].pct_pos.sum() if ((R.lib==lb)&(R.cl==k)).any() else np.nan for k in cls]
    ax.bar(np.arange(len(cls))+(j-0.5)*0.4,v,0.4,color=c,label=lb.replace('_',' ')+' d')
ax.set_xticks(range(len(cls))); ax.set_xticklabels([f'C{k}' for k in cls],fontsize=5); ax.set_ylabel('% nuclei with OLE16a UMI (clusters ≥30 nuclei)'); ax.legend(fontsize=5,frameon=False)
ax.set_title('e  snRNA-seq, 185 d',loc='left',fontsize=7)
fig.savefig('ED_draft/ED_oleosin_draft.png',dpi=300,bbox_inches='tight'); fig.savefig('ED_draft/ED_oleosin_draft.pdf',bbox_inches='tight')
