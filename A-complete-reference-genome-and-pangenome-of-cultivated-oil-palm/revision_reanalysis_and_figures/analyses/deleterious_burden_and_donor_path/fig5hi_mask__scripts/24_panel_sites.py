#!/usr/bin/env python3
"""(1) Reproduction: Nigerian ALT-derived conserved dSNP site sets from the re-run on the 7/25 assemblies vs the
carriers recorded in results-8.9 alt_derived_conserved_snp.tsv (the set used for Fig. 5h,i).
(2) African35 (and All38) frequency-inclusive site counts after replacing the Nigerian carriers with the calls on
the final assemblies (Derived_Site_Count, Fixed_Unavoidable_Site_Count; 70_prepare_freqinclusive_dsnp.py rules)."""
import csv, sys
from collections import defaultdict
W = "${CLUSTER_WORK}/fig5hi_mask"
R = "${ANALYSIS_DIR}/21_MS/06_result/dSVs/results-8.9"
TAB = R + "/09_ideal_parent_haplotypes/data/phoenix_polarity_conserved/alt_derived_conserved_snp.tsv"
MAN = R + "/config/Sample_Manifest.tsv"
man = list(csv.DictReader(open(MAN), delimiter="\t"))
A35 = {r["Sample_ID"] for r in man if r["Include_African35"].strip().lower() in ("yes", "true", "1")}
A38 = {r["Sample_ID"] for r in man if r["Include_All38"].strip().lower() in ("yes", "true", "1")}
assert len(A35) == 35 and len(A38) == 38
NIG = ("nrly_hap1", "nrly_hap2")
def load_sites(lab):
    s = set()
    for r in csv.DictReader(open(f"{W}/dsnp/{lab}.sites.tsv"), delimiter="\t"):
        if r["Polarity"] == "ALT_Derived": s.add((r["Chrom"], int(r["Pos"]), r["Ref"], r["Alt"]))
    return s
carriers = {}; old = {h: set() for h in NIG}
with open(TAB) as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        k = (r["Chrom"], int(r["Pos"]), r["Ref"], r["Alt"]); c = set(r["Samples"].split(";"))
        carriers[k] = c
        for h in NIG:
            if h in c: old[h].add(k)
out = open(f"{W}/dsnp/panel_site_counts.tsv", "w")
out.write("Item\tValue\n")
for h in NIG:
    out.write(f"catalog_{h}_sites\t{len(old[h])}\n")
    for tag in ("nrly_old_" + h[-4:], "nrly_final_" + h[-4:]):
        try: s = load_sites(tag)
        except FileNotFoundError: continue
        out.write(f"{tag}_sites\t{len(s)}\n{tag}_shared_with_catalog_{h}\t{len(s & old[h])}\n"
                  f"{tag}_only\t{len(s - old[h])}\ncatalog_{h}_only_vs_{tag}\t{len(old[h] - s)}\n")
# replace Nigerian carriers with the final-assembly calls
new = {}
try:
    new = {h: load_sites("nrly_final_" + h[-4:]) for h in NIG}
except FileNotFoundError:
    pass
if len(new) == 2:
    c2 = {k: set(v) - set(NIG) for k, v in carriers.items()}
    for h in NIG:
        for k in new[h]: c2.setdefault(k, set()).add(h)
    for lab, cmap in (("original", carriers), ("final_nigerian", c2)):
        for pan, S in (("African35", A35), ("All38", A38)):
            der = fix = 0; priv = defaultdict(int)
            for k, v in cmap.items():
                sel = v & S
                if not sel: continue
                der += 1
                if sel == S: fix += 1
                if len(sel) == 1: priv[next(iter(sel))] += 1
            out.write(f"{lab}_{pan}_Derived_Site_Count\t{der}\n{lab}_{pan}_Fixed_Unavoidable_Site_Count\t{fix}\n")
            for h in NIG: out.write(f"{lab}_{pan}_private_{h}\t{priv[h]}\n")
out.close()
print(open(f"{W}/dsnp/panel_site_counts.tsv").read())
