#!/usr/bin/env python3
"""Recompute favourable-locus capture on the current African35 W=0, P=15 path.
Definition copied from results/09_ideal_parent_haplotypes/scripts/09_mosaic_step9.py::captured():
  window = Pos // 500000 (FL-Hap2 = Africa_hap2 coordinates, same reference as the current path);
  selected donor = path[chrom][window]; ca = step9 cache carries_ALT(sv, donor);
  captured iff ca is known and (Target==ALT and ca) or (Target==REF and not ca); unknown -> not captured.
Donor names: current path uses phased dura_hap1/2, pisifera_hap1/2, nrly_hap1/2 that have no call in the
legacy cache (09b Sample_Crosswalk: Unknown_Pending_Assembly_Genotyping). Sensitivity: collapsed proxy
EG_dura / EG_pisifera / EG_niriliya (audit-only).
"""
import csv, os
from collections import defaultdict, Counter
H = os.path.dirname(os.path.abspath(__file__))
M9 = os.path.join(H, "..", "misc9", "data")
WINDOW = 500_000
PROXY = {"dura_hap1": "EG_dura", "dura_hap2": "EG_dura", "pisifera_hap1": "EG_pisifera",
         "pisifera_hap2": "EG_pisifera", "nrly_hap1": "EG_niriliya", "nrly_hap2": "EG_niriliya"}

def rd(p): return list(csv.DictReader(open(p), delimiter="\t"))
path = defaultdict(dict)
for r in rd(f"{M9}/ideal_loadonly_path.tsv"): path[r["Chrom"]][int(r["Window_Index"])] = r["Donor_ID"]
alle = defaultdict(dict)
for r in rd(f"{H}/src/step9_donor_alleles.tsv"): alle[r["sv_id"]][r["donor"]] = r["carries_ALT"] == "1"
# exact-sequence states (09b audit) as a second sensitivity
exact = defaultdict(dict)
for r in rd(f"{H}/src/Current_All38_Direct_Allele_Matrix.tsv"): exact[r["GWAS_SV_ID"]][r["Sample_ID"]] = r["Formal_Allele_State"]
reg = {r["GWAS_SV_ID"]: r for r in rd(f"{H}/src/GWAS_Target_Registry.tsv")}
loci = rd(f"{M9}/ideal_step9_ALL_loci.tsv")
chromlen = {c: max(path[c]) for c in path}

def state(sv, tgt, donor):
    ca = alle[sv].get(donor)
    if ca is None: return None
    return int((tgt == "ALT" and ca) or (tgt == "REF" and not ca))

out = []
for r in loci:
    sv, ch, pos, tgt = r["SV"], r["Chrom"], int(r["Pos"]), r["Target"]
    w = pos // WINDOW
    proj = ch in path and w in path[ch]
    don = path[ch][w] if proj else "NA"
    s = state(sv, tgt, don) if proj else None
    status = "Captured" if s == 1 else ("Missed" if s == 0 else "Unknown_no_donor_call")
    cap = 1 if s == 1 else 0
    pd_ = PROXY.get(don)
    sp = state(sv, tgt, pd_) if pd_ else s
    # exact-sequence (09b) state
    ex = exact.get(sv, {}).get(don, "NA")
    ta = reg.get(sv, {}).get("Core_Target_Allele", "NA")
    out.append(dict(SV=sv, Chrom=ch, Pos=pos, Window_Index=w, Window_Start_0based=w*WINDOW,
        Target=tgt, Projectable=int(proj), Current_W0_Selected_Donor=don,
        Donor_Call_Source=("Legacy_step9_cache" if don in {d for v in alle.values() for d in v} else "None_phased_assembly_not_genotyped"),
        Selected_Donor_Carries_ALT=("NA" if alle[sv].get(don) is None else int(alle[sv][don])),
        Captured_W0=cap, Capture_Status_W0=status,
        Proxy_Donor=pd_ or "", Captured_W0_with_collapsed_proxy=("NA" if sp is None else sp),
        Exact_Sequence_State_09b=ex,
        Legacy_W4_Selected_Donor=r["Selected_Donor"], Legacy_W4_Captured=int(r["Captured"])))
flds = list(out[0])
with open(f"{H}/capture_current_path.tsv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=flds, delimiter="\t", lineterminator="\n"); w.writeheader(); w.writerows(out)

chroms = sorted(path, key=lambda c: int(c[3:5]))
per = {c: Counter() for c in chroms}
for o in out:
    p = per[o["Chrom"]]; p["tot"] += 1; p["cap"] += o["Captured_W0"]; p["unk"] += o["Capture_Status_W0"].startswith("Unknown")
    p["miss"] += o["Capture_Status_W0"] == "Missed"; p["old"] += o["Legacy_W4_Captured"]
    p["prox"] += 1 if o["Captured_W0_with_collapsed_proxy"] == 1 else 0
    p["prox_unk"] += o["Captured_W0_with_collapsed_proxy"] == "NA"
with open(f"{H}/per_chrom_capture.tsv", "w") as f:
    f.write("Chrom\tFav_Total\tCaptured_W0\tMissed_W0\tUnknown_W0_no_donor_call\tCapture_Rate_W0_pct\tCaptured_W0_collapsed_proxy\tUnknown_after_proxy\tLegacy_W4_Captured\tDelta_W0_minus_W4\n")
    T = Counter()
    for c in chroms:
        p = per[c]; T.update(p)
        f.write(f"{c}\t{p['tot']}\t{p['cap']}\t{p['miss']}\t{p['unk']}\t{100*p['cap']/p['tot'] if p['tot'] else 0:.2f}\t{p['prox']}\t{p['prox_unk']}\t{p['old']}\t{p['cap']-p['old']}\n")
    f.write(f"TOTAL\t{T['tot']}\t{T['cap']}\t{T['miss']}\t{T['unk']}\t{100*T['cap']/T['tot']:.2f}\t{T['prox']}\t{T['prox_unk']}\t{T['old']}\t{T['cap']-T['old']}\n")
# position-level marks (as the plotting script: max over loci sharing a coordinate)
mark = {}
for o in out:
    k = (o["Chrom"], o["Pos"]); mark[k] = max(mark.get(k, 0), o["Captured_W0"])
oldmark = {}
for o in out:
    k = (o["Chrom"], o["Pos"]); oldmark[k] = max(oldmark.get(k, 0), o["Legacy_W4_Captured"])
with open(f"{H}/position_marks.tsv", "w") as f:
    f.write("Chrom\tPos\tMark_W0\tMark_Legacy_W4\tChanged\n")
    for k in sorted(mark, key=lambda k: (int(k[0][3:5]), k[1])):
        f.write(f"{k[0]}\t{k[1]}\t{mark[k]}\t{oldmark[k]}\t{int(mark[k]!=oldmark[k])}\n")
print("positions", len(mark), "gold", sum(mark.values()), "grey", len(mark)-sum(mark.values()),
      "| old gold", sum(oldmark.values()), "changed", sum(mark[k]!=oldmark[k] for k in mark))
print("TOTAL", T)
print("unprojectable", sum(1 for o in out if not o["Projectable"]))
print(Counter(o["Current_W0_Selected_Donor"] for o in out if o["Capture_Status_W0"].startswith("Unknown")))
