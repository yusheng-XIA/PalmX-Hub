"""Build CDS model table (one mRNA per gene: longest CDS) and per-chrom CDS BED from the dSV-pipeline FL-Hap2 annotation."""
import re, collections, sys
B='${ANALYSIS_DIR}/21_MS/06_result/dSVs/input'
O='${CLUSTER_WORK}/headline_pop/g6'
import os; os.makedirs(O+'/bed',exist_ok=True)
cds=collections.defaultdict(list); gene_of={}
for l in open(B+'/Africa_hap2.EVM.gff3'):
    if l.startswith('#'): continue
    f=l.rstrip('\n').split('\t')
    if len(f)<9: continue
    if f[2]=='mRNA':
        m=re.search('ID=([^;]+)',f[8]).group(1); p=re.search('Parent=([^;]+)',f[8]).group(1); gene_of[m]=p
    elif f[2]=='CDS':
        p=re.search('Parent=([^;]+)',f[8]).group(1); cds[p].append((f[0],int(f[3]),int(f[4]),f[6]))
best={}
for m,ex in cds.items():
    g=gene_of.get(m,m); L=sum(e[2]-e[1]+1 for e in ex)
    if g not in best or L>best[g][0]: best[g]=(L,m)
out=open(O+'/cds_models.tsv','w'); bed=collections.defaultdict(list)
for g,(L,m) in best.items():
    ex=sorted(cds[m],key=lambda e:e[1]); c=ex[0][0]; s=ex[0][3]
    out.write(f"{m}\t{c}\t{s}\t{','.join(f'{a}-{b}' for _,a,b,_ in ex)}\n")
    for _,a,b,_ in ex: bed[c].append((a-1,b))
for c,v in bed.items():
    v.sort(); mg=[]
    for a,b in v:
        if mg and a<=mg[-1][1]: mg[-1]=(mg[-1][0],max(b,mg[-1][1]))
        else: mg.append((a,b))
    with open(f'{O}/bed/{c}.bed','w') as fh:
        for a,b in mg: fh.write(f'{c}\t{a}\t{b}\n')
print(len(best),'genes', sum(len(v) for v in bed.values()))
