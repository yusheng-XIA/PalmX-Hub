import re, pandas as pd
R="${ANALYSIS_DIR}/22_answer_reviews/00_ms/03_V3/03_figure3/07_proteomics_multiomics_20260724/01_reference_database"
O="out/dom/"
def readfa(p):
    d={};k=None
    for l in open(p):
        l=l.strip()
        if l.startswith('>'): k=l[1:].split()[0]; d[k]=[]
        elif k: d[k].append(l)
    return {k:''.join(v).rstrip('*') for k,v in d.items()}
seq=readfa(R+"/FL_TN_unified_exact_sequence_nr.fasta")
ref=readfa("ref_oleosins.fa")
mp=pd.read_csv(R+"/FL_TN_unified_exact_sequence_nr_map.tsv",sep="\t",dtype=str).set_index("unified_protein_id")
def dom(p):
    rows=[]
    for l in open(p):
        if l.startswith('#'): continue
        f=l.split()
        rows.append(dict(target=f[0],tlen=int(f[2]),pfam=f[3],acc=f[4],evalue=float(f[6]),score=float(f[7]),i_eval=float(f[12]),hmm_from=int(f[15]),hmm_to=int(f[16]),ali_from=int(f[17]),ali_to=int(f[18])))
    return pd.DataFrame(rows)
d=dom(O+"unified_ld3.domtbl")
best=d.sort_values("score",ascending=False).groupby("target").head(1)
KD=dict(A=1.8,R=-4.5,N=-3.5,D=-3.5,C=2.5,Q=-3.5,E=-3.5,G=-0.4,H=-3.2,I=4.5,L=3.8,K=-3.9,M=1.9,F=2.8,P=-1.6,S=-0.8,T=-0.7,W=-0.9,Y=-1.3,V=4.2)
def hairpin(s):
    # longest run where 9-aa Kyte-Doolittle window mean >= 1.6 (hydrophobic core)
    w=9; v=[sum(KD.get(c,0) for c in s[i:i+w])/w for i in range(max(0,len(s)-w+1))]
    best=(0,0,0);st=None
    for i,x in enumerate(v+[-9]):
        if x>=1.6 and st is None: st=i
        if x<1.6 and st is not None:
            L=i-st+w-1
            if L>best[0]: best=(L,st+1,i+w-1)
            st=None
    return best
PK=re.compile(r"P.{5}SP.{3}P")
recs=[]
for r in best.itertuples():
    s=seq[r.target]; L,a,b=hairpin(s); m=PK.search(s)
    info=mp.loc[r.target]
    recs.append(dict(unified_id=r.target,family={"Oleosin":"oleosin","Caleosin":"caleosin","SRP":"LDAP_REF"}.get(r.pfam,r.pfam),pfam=r.pfam,score=r.score,evalue=r.evalue,
        len=len(s),hydrophobic_core_len=L,core_from=a,core_to=b,proline_knot=m.group(0) if m else "",pk_pos=(m.start()+1) if m else "",
        pk_in_core=bool(m and a<=m.start()+1<=b),materials=info.materials,sources=info.sources,genes=info.source_gene_ids))
T=pd.DataFrame(recs).sort_values(["family","genes"]); T.to_csv(O+"ld_candidates_domains.tsv",sep="\t",index=False)
print(T.family.value_counts()); print(T[T.family=="oleosin"].to_string())
# references
dr=dom(O+"ref_ld3.domtbl")
rr=[]
for k,s in ref.items():
    L,a,b=hairpin(s); m=PK.search(s); rr.append(dict(ref=k,len=len(s),core=L,pk=m.group(0) if m else "",pfam_hit=k in set(dr.target)))
print(pd.DataFrame(rr).to_string())
pd.DataFrame(rr).to_csv(O+"ref_domains.tsv",sep="\t",index=False)
with open(O+"oleosin_tree_input.fa","w") as f:
    for r in T[T.family=="oleosin"].itertuples():
        f.write(f">{r.unified_id}|{r.materials}|{r.genes.split(';')[0].split(':')[1]}\n{seq[r.unified_id]}\n")
    for k,s in ref.items(): f.write(f">{k}\n{s}\n")
    for r in T[T.family=="oleosin"].itertuples(): pass
