"""Recompute Supplementary Fig. 2b (CAFE5 changes at the oil palm ancestral node, per lipid enzyme) from a CAFE5
gamma run. Linked orthogroups and enzyme labels are the authors' (unique_OG_audit.tsv); the node is the MRCA of
American_hap1, Dura and Pisifera. Validated by reproducing the current SF2b from the authors' run."""
import csv, re, sys
from collections import OrderedDict
D = "${ANALYSIS_DIR}/22_answer_reviews/00_ms/03_V3/02_figure/03_Figure2_panels_flat_20260727/"
AUD = D + "Fig2b_core_lipid_copy_number_rigorous_v2_unique_OG_audit.tsv"
REF = D + "Fig2b_core_lipid_copy_number_rigorous_v2_ancestral_changes.tsv"
run, out = sys.argv[1], sys.argv[2]
line = next(l for l in open(run + "/Gamma_asr.tre") if "TREE" in l and "=" in l)
m = re.search(r"\(American_hap1<\d+>[^()]*\([^()]*\)<\d+>[^()]*\)<(\d+)>", line)
node = m.group(1)
hdr, *rows = [l.rstrip("\n").split("\t") for l in open(run + "/Gamma_change.tab")]
col = next(i for i, h in enumerate(hdr) if h.endswith(f"<{node}>") and not re.match(r"[A-Za-z]", h))
chg = {r[0]: int(r[col]) for r in rows}
sig = {l.split("\t")[0]: l.split("\t")[2].strip() == "y" for l in open(run + "/Gamma_family_results.txt") if not l.startswith("#")}
ref = list(csv.DictReader(open(REF), delimiter="\t"))
aud = list(csv.DictReader(open(AUD), delimiter="\t"))
res = []
for r in ref:
    e = r["Enzyme"]
    ogs = [a["Orthogroup"] for a in aud if e in [x.strip() for x in a["Enzyme_labels"].split(";")]]
    c = [chg.get(o, 0) for o in ogs]
    res.append(OrderedDict(Section=r["Section"], Enzyme=e, N_linked_OGs=len(ogs), N_expanded_OGs=sum(x > 0 for x in c),
               Total_inferred_gains=sum(x for x in c if x > 0), N_contracted_OGs=sum(x < 0 for x in c),
               Total_inferred_losses=-sum(x for x in c if x < 0), N_unchanged_OGs=sum(x == 0 for x in c),
               Net_change=sum(c), N_familywide_significant_OGs=sum(sig.get(o, False) for o in ogs),
               Familywide_significance_is_branch_specific="False"))
with open(out, "w") as f:
    w = csv.DictWriter(f, fieldnames=list(res[0]), delimiter="\t", lineterminator="\n"); w.writeheader(); w.writerows(res)
print("node", node, "column", hdr[col], "OGs missing from run:", sum(1 for a in aud if a["Orthogroup"] not in chg))
