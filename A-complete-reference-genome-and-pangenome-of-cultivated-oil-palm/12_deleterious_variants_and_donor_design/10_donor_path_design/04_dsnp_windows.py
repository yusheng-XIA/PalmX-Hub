#!/usr/bin/env python3
"""Per-sample dSNP extraction replicating stages 05 (40_build_dsnp_hap38_catalog.py::normalize_records),
09 Job A (70_prepare_freqinclusive_dsnp.py: Filter=PASS, Conserved_Region_Proxy, Phoenix ALT_Derived, MAPQ>=20)
and 09 Job B (71_build_loadonly_iph.py: window = (Pos-1)//500000, capped at last window).
usage: 21_dsnp_windows.py LABEL SNP_TXT OUT_PREFIX
writes OUT_PREFIX.sites.tsv (Chrom Pos Ref Alt) and OUT_PREFIX.windows.tsv (Chrom Window_Index DSNP_Count)
"""
import sys, csv, bisect, importlib.util
from collections import Counter, defaultdict
import os
D = os.environ.get("DSV_DIR", "dsv_analysis")
REF_FA = D + "/input/Africa_hap2.fa"; REF_FAI = D + "/input/Africa_hap2.fa.fai"
SNP_FAI = "FL_Hap2.fasta.fai"
BED = D + "/results/archive_old_versions/05_dsv_v3_multi_outgroup/multi_outgroup_conserved_regions/multi_outgroup_conserved_regions.bed"
PHX = D + "/results_hap38/06_dSNP_phoenix_polarity/alignments/Phoenix_vs_Africa_hap2.asm20.cs.primary.paf"
W = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # 12_deleterious_variants_and_donor_design/
spec = importlib.util.spec_from_file_location("s70", W + "/08_freq_inclusive_dsnp.py")
s70 = importlib.util.module_from_spec(spec); spec.loader.exec_module(s70)
spec = importlib.util.spec_from_file_location("s40", W + "/03_build_dsnp_catalog.py")
s40 = importlib.util.module_from_spec(spec); spec.loader.exec_module(s40)
WIN = 500_000
label, snp_path, outp = sys.argv[1:4]
lengths = s40.read_fai(REF_FAI); chrom_map = s40.build_chrom_map(s40.read_fai(SNP_FAI), lengths)
seqs = s40.read_fasta(REF_FA, set(lengths))
cons = s40.read_bed_index(BED, 100000, label="Conserved_Region")
cnt = Counter(); seen = set(); keep = []
with open(snp_path) as fh:
    for line in fh:
        if not line.strip() or line.startswith("#"): continue
        f = line.rstrip("\n").split("\t")
        cnt["raw"] += 1
        if len(f) < 8: cnt["malformed"] += 1; continue
        chrom = chrom_map.get(f[0])
        if not chrom: cnt["unmapped"] += 1; continue
        ref, alt = f[3].upper(), f[4].upper(); pos = s40.as_int(f[1], None)
        if pos is None or pos < 1 or pos > lengths[chrom]: cnt["oob"] += 1; continue
        if len(ref) != 1 or len(alt) != 1 or ref not in "ACGT" or alt not in "ACGT": cnt["nonsnv"] += 1; continue
        rb = seqs[chrom][pos - 1]
        if rb == "N": cnt["refN"] += 1; continue
        if rb != ref: cnt["ref_mismatch"] += 1; continue
        k = (chrom, pos, ref, alt)
        if k in seen: cnt["dup"] += 1; continue
        seen.add(k); cnt["valid"] += 1
        if f[6] != "PASS": cnt["nonpass"] += 1; continue
        if not cons.hits(chrom, pos - 1): continue
        cnt["conserved"] += 1
        keep.append(k)
del seqs
rows = [(f"S{i}", c, p, r, a, 1, "0", label, (label,)) for i, (c, p, r, a) in enumerate(keep)]
calls, pst = s70.parse_paf(PHX, rows, 20)
win = Counter(); nder = 0
with open(outp + ".sites.tsv", "w") as out:
    out.write("Chrom\tPos\tRef\tAlt\tPolarity\n")
    for i, row in enumerate(rows):
        pol = s70.summarize_call(row, calls.get(i, []))["ALT_Polarity_Phoenix"]
        out.write(f"{row[1]}\t{row[2]}\t{row[3]}\t{row[4]}\t{pol}\n")
        if pol != "ALT_Derived": continue
        nder += 1
        nwin = (lengths[row[1]] + WIN - 1) // WIN
        win[(row[1], min((row[2] - 1) // WIN, nwin - 1))] += 1
with open(outp + ".windows.tsv", "w") as out:
    out.write("Chrom\tWindow_Index\tDSNP_Count\n")
    for c in sorted(lengths, key=lambda c: int(c[3:5])):
        for w in range((lengths[c] + WIN - 1) // WIN):
            out.write(f"{c}\t{w}\t{win[(c, w)]}\n")
cnt["alt_derived"] = nder
with open(outp + ".audit.tsv", "w") as out:
    for k, v in cnt.items(): out.write(f"{k}\t{v}\n")
print(label, dict(cnt))
