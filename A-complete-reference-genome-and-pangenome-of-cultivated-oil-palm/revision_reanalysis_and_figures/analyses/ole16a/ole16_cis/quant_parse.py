import sys,os,re,collections,glob,csv
fam={};uid={}
for l in open('target_family.tsv'):
    k,u,f=l.rstrip('\n').split('\t'); fam[k]=f; uid[k]=u
L={}
for l in open('targets.bed'):
    f=l.split('\t'); L[f[0]]=int(f[3])-int(f[2])+1
tot={}
for f in glob.glob('dl/*.count'):
    r=os.path.basename(f)[:-6]
    try: tot[r]=int(open(f).read().strip() or 0)
    except: pass
meta={}
for fn,cols in [('prjeb11097.tsv',('run_accession','scientific_name','sample_alias')),('prjdb1773.tsv',('run_accession','scientific_name','sample_title'))]:
    for r in csv.DictReader(open(fn),delimiter='\t'): meta[r['run_accession']]=(r['scientific_name'],r[cols[2]])
OLEa={'African_hap2__evm.TU.chr11B.1497':'Eg','American_hap1__evm.TU.chr11A.1059':'Eo'}
rows=[]
for sam in sorted(glob.glob('rq/*.sam')):
    r=os.path.basename(sam)[:-4]
    best={}; alle=collections.defaultdict(dict)
    for l in open(sam):
        f=l.rstrip('\n').split('\t')
        q,flag,t=f[0],int(f[1]),f[2]
        if t=='*' or flag&4: continue
        cig=f[5]; alen=sum(int(n) for n,o in re.findall(r'(\d+)([MI=X])',cig))
        tags=dict((x[:2],x[5:]) for x in f[6:] if len(x)>5)
        nm=int(tags.get('NM',99)); AS=int(tags.get('AS',-999))
        if alen<60 or nm>0.03*alen: continue
        if t in OLEa: alle[q][OLEa[t]]=max(AS,alle[q].get(OLEa[t],-999))
        if flag&256 or flag&2048: continue
        best[q]=t
    cnt=collections.Counter(uid.get(t,t) for t in best.values())
    famc=collections.Counter(fam.get(t,'?') for t in best.values())
    eo=sum(1 for q,d in alle.items() if 'Eo' in d and d.get('Eo',-999)>d.get('Eg',-999))
    eg=sum(1 for q,d in alle.items() if 'Eg' in d and d.get('Eg',-999)>d.get('Eo',-999))
    T=tot.get(r,0) or 1
    sp,lab=meta.get(r,('?','?'))
    rows.append(dict(run=r,species=sp,sample=lab,total_reads=T,
        OLE16a_rpm=round(cnt['UFTN008330']/T*1e6,1),OLE16b_rpm=round(cnt['UFTN008332']/T*1e6,1),
        LDAP_chr10_rpm=round(cnt['UFTN008529']/T*1e6,1),ACTIN_rpm=round(cnt['ACTIN']/T*1e6,1),
        oleosin_all_rpm=round(famc['oleosin']/T*1e6,1),LDAP_all_rpm=round(famc['REF']/T*1e6,1),caleosin_all_rpm=round(famc['caleosin']/T*1e6,1),
        OLE16a_Eo_allele_reads=eo,OLE16a_Eg_allele_reads=eg))
w=csv.DictWriter(open('public_rna_quant.tsv','w'),fieldnames=list(rows[0].keys()),delimiter='\t'); w.writeheader(); [w.writerow(x) for x in rows]
for x in rows: print('\t'.join(str(v) for v in x.values()))
