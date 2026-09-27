A=${ANALYSIS_DIR}/08_hifi_chromosome/11_final_corrected_39genomes_annotations_20260812/attempt_20260812_01
O=${CLUSTER_WORK}/headline_pop/out
ls $A; ls $A/* | head -40
P=$(ls $A/*/African_hap2*pep* $A/*/African_hap2*prot* $A/*/African_hap2*.pep.fa* 2>/dev/null | head -1)
echo PEP=$P
if [ -n "$P" ]; then
python3 - "$P" <<'PY' > $O/g3_pep.txt
import sys,re
seqs={};name=None
for l in open(sys.argv[1]):
    if l.startswith('>'): name=l[1:].split()[0]; seqs[name]=[]
    else: seqs[name].append(l.strip())
for n,s in seqs.items():
    m=re.search(r'chr01B?\.(\d+)$',n)
    if m and 140<=int(m.group(1))<=180:
        s=''.join(s); print(n,len(s),s[:80], 'MADS' if re.search('KIEIKRIE|RGKIEIK|GRGKIE',s) else '')
PY
cat $O/g3_pep.txt
fi
