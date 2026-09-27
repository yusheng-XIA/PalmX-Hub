#!/bin/bash
W=${CLUSTER_WORK}/ole16_cis; cd $W/prom; export TMPDIR=$W/tmp
MM=minimap2
MAFFT=mafft
python3 - <<'P'
def readfa(p):
    s={};k=None
    for l in open(p):
        l=l.strip()
        if l.startswith('>'): k=l[1:]; s[k]=[]
        elif k: s[k].append(l)
    return {k:''.join(v) for k,v in s.items()}
P=readfa('promoters.fa')
for ref in ['African_hap2','MZ4_hap1']:
    open(f'ref_{ref}.fa','w').write(f'>{ref}\n{P[ref]}\n')
with open('prox3k.fa','w') as f:
    for g,s in P.items(): f.write(f'>{g}\n{s[-3000:]}\n')
P
for ref in African_hap2 MZ4_hap1; do $MM -c --cs -x asm20 -t 8 ref_$ref.fa promoters.fa > vs_$ref.paf 2>/dev/null; done
$MM -X -c -x asm20 -t 8 promoters.fa promoters.fa > allvsall.paf 2>/dev/null
$MAFFT --auto --thread 8 --quiet prox3k.fa > prox3k.aln.fa
$MAFFT --auto --thread 8 --quiet genes.fa > genes.aln.fa
echo ok
