import subprocess, os, re, sys
W='${CLUSTER_WORK}/ole16_cis'; os.chdir(W)
MM='minimap2'
def readfa(p):
    s={};k=None
    for l in open(p):
        l=l.strip()
        if l.startswith('>'): k=l[1:]; s[k]=[]
        elif k: s[k].append(l)
    return {k:''.join(v) for k,v in s.items()}
def rc(x): return x[::-1].translate(str.maketrans('ACGTNacgtn','TGCANtgcan'))
def hits(target,query,preset):
    out=subprocess.run([MM,'-c','--cs']+preset+[target,query],stdout=subprocess.PIPE,stderr=subprocess.DEVNULL,universal_newlines=True).stdout
    H=[]
    for l in out.splitlines():
        f=l.split('\t'); H.append(dict(qs=int(f[2]),qe=int(f[3]),st=f[4],ts=int(f[7]),te=int(f[8]),nm=int(f[9]),al=int(f[10]),mq=int(f[11]),ql=int(f[1])))
    return sorted(H,key=lambda h:-h['nm'])
gs=[l.split()[0] for l in open('genomes.txt')]
os.makedirs('prom',exist_ok=True)
rows=[];P={};G={};D={}
short=['-x','asm20','-k','11','-w','5','-s','40']
for g in gs:
    lf=f'loc/{g}.locus.fa'
    if not os.path.exists(lf): rows.append([g,'NOLOCUS']);continue
    seq=list(readfa(lf).values())[0]
    ho=hits(lf,'q/ole16a_gene_fwdstrand.fa',short)
    h98=hits(lf,'q/g1498.fa',['-x','asm20'])
    h96=hits(lf,'q/g1496.fa',['-x','asm20'])
    if not ho: rows.append([g,'NO_OLE16A']);continue
    o=ho[0]
    # gene on target: query is FL-Hap2 + strand where gene is minus; hit '+' => gene minus on target
    gstrand='-' if o['st']=='+' else '+'
    ident=o['nm']/o['al'] if o['al'] else 0
    # neighbour G1498 hits near o (union of blocks, take span)
    def span(H,maxd=60000):
        H=[h for h in H if min(abs(h['ts']-o['te']),abs(h['te']-o['ts']))<maxd]
        if not H: return None
        return min(h['ts'] for h in H),max(h['te'] for h in H),sum(h['nm'] for h in H)
    s98=span(h98); s96=span(h96)
    if gstrand=='-':
        tss=o['te']  # 5' end at higher coord
        up_end = s98[0] if s98 and s98[0]>tss else None
        prom = rc(seq[tss:up_end]) if up_end else rc(seq[tss:tss+15000])
        down = rc(seq[max(0,o['ts']-2000):o['ts']])
        gene = rc(seq[o['ts']:o['te']])
    else:
        tss=o['ts']
        up_end = s98[1] if s98 and s98[1]<tss else None
        prom = seq[up_end:tss] if up_end else seq[max(0,tss-15000):tss]
        down = seq[o['te']:o['te']+2000]
        gene = seq[o['ts']:o['te']]
    # prom is oriented 5'->3' ending at the gene 5' end
    P[g]=prom; G[g]=gene; D[g]=down
    dist96 = (o['ts']-s96[1]) if (s96 and gstrand=='-') else ((s96[0]-o['te']) if s96 else None)
    rows.append([g,'OK',gstrand,o['al'],round(ident,4),len(prom),'bounded' if up_end else 'open15k',dist96,(s98[2] if s98 else 0),(s96[2] if s96 else 0)])
with open('prom/promoters.fa','w') as f:
    for g,s in P.items(): f.write(f'>{g}\n{s}\n')
with open('prom/genes.fa','w') as f:
    for g,s in G.items(): f.write(f'>{g}\n{s}\n')
with open('prom/down2k.fa','w') as f:
    for g,s in D.items(): f.write(f'>{g}\n{s}\n')
with open('prom/summary.tsv','w') as f:
    f.write('genome\tstatus\tgene_strand\tgene_aln_len\tgene_identity_vs_FLHap2\tupstream_intergenic_bp\tupstream_bound\tdownstream_gap_bp\tG1498_matches\tG1496_matches\n')
    for r in rows: f.write('\t'.join(map(str,r))+'\n')
print(open('prom/summary.tsv').read())
