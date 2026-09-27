#!/usr/bin/env python3
"""Contig Nx for the 8 RagTag haplotypes (run on ${COMPUTE_HOST}).
Contigs = scaffold sequences split at every run of >=1 N/n (same regex as
20_results/Figure1/09_new_figure/02_N50/plot_Nx_curves_full.py::read_fasta_as_contigs,
which also drops pieces <1,000 bp; we report how many such pieces exist)."""
import re, sys, csv
BASE = '${ANALYSIS_DIR}/22_answer_reviews/00_ms/03_V3/01_figure1/00_final_input_resources/04_ragtag_6haps/fasta'
S = ['dura_hap1','dura_hap2','pisifera_hap1','pisifera_hap2','nrly_hap1','nrly_hap2','meizhou4_hap1','meizhou4_hap2']
rn = re.compile(r'[Nn]+')
def nx(L):
    L = sorted(L, reverse=True); tot = sum(L); t = [tot*x/100 for x in range(10,101,10)]; out=[]; c=0; i=0
    for l in L:
        c += l
        while i < 10 and c >= t[i]: out.append(l); i += 1
    return out, tot, len(L), L[0]
w = csv.writer(open(sys.argv[1], 'w'), delimiter='\t')
w.writerow(['Sample','level']+[f'N{k}' for k in range(10,101,10)]+['Total','N_seqs','Largest','N_gap_runs','N_bases','min_gap_run','max_gap_run','pieces_lt_1kb'])
for s in S:
    scaf=[]; ctg=[]; gl=[]
    def flush(seq):
        if not seq: return
        x=''.join(seq); scaf.append(len(x))
        gl.extend(len(m.group()) for m in rn.finditer(x))
        ctg.extend(len(p) for p in rn.split(x) if p)
    seq=[]
    with open(f'{BASE}/{s}.ragtag.fasta') as f:
        for line in f:
            if line[0]=='>': flush(seq); seq=[]
            else: seq.append(line.strip())
    flush(seq)
    small = sum(1 for l in ctg if l < 1000)
    for lvl,L in (('scaffold',scaf),('contig_splitN',ctg)):
        o,tot,n,big = nx(L)
        w.writerow([s,lvl]+o+[tot,n,big,len(gl),sum(gl),min(gl) if gl else 0,max(gl) if gl else 0,small])
    print(s,'done',flush=True)
