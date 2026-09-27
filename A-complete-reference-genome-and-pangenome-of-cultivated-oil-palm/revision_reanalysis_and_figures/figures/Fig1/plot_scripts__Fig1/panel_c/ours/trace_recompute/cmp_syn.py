import csv, collections
W='${CLUSTER_WORK}/trace/Fig1c'
sd=[r for r in csv.reader(open(W+'/sd/Fig1c_synteny.tsv'),delimiter='\t')][1:]
sdk={(r[0],int(float(r[1])),int(float(r[2])),r[3],int(float(r[4])),int(float(r[5])),r[6],int(float(r[8]))):int(float(r[7])) for r in sd}
mine={}
for l in open(W+'/out/sam_blocks.tsv'):
    f=l.split('\t'); rec=(f[0],int(f[1]),int(f[2]),f[3],int(f[4]),int(f[5]),f[6],int(f[8])); mine.setdefault(rec,[]).append(int(f[7]))
ok=[k for k in sdk if k in mine]; print('coord match',len(ok),'of',len(sdk))
d=collections.Counter(sdk[k]-min(mine[k],key=lambda x:abs(x-sdk[k])) for k in ok); print('blocklen diff (SD-mine) top',d.most_common(8))
miss=[k for k in sdk if k not in mine]; print('missing',len(miss),miss[:3])
# blocks in mine passing filter but not in SD by coords
m2=[k for k,v in mine.items() if k[0]==k[3] and k[7]>=20 and max(v)>=50000 and k not in sdk]; print('extra',len(m2),m2[:5])
