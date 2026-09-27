#!/usr/bin/env python3
"""Reproduction of the original Nigerian dSNP input from the re-run on the 7/25 assemblies:
(1) the full reference-checked SNV set of each re-run equals the set of stage-05 catalogue rows carrying that haplotype;
(2) the 2,285 / 2,549 rare ALT-derived dSNPs (dsnp_v2_phoenix_alt_derived_candidates.tsv, Samples == nrly_hapX) are all
    recovered; (3) same comparison for the final assemblies (informative)."""
import csv, importlib.util, sys
W = "${CLUSTER_WORK}/fig5hi_mask"
R = "${ANALYSIS_DIR}/21_MS/06_result/dSVs/results-8.9"
D = "${ANALYSIS_DIR}/21_MS/06_result/dSVs"
spec = importlib.util.spec_from_file_location("s40", W + "/scripts/author/40_build_dsnp_hap38_catalog.py")
s40 = importlib.util.module_from_spec(spec); spec.loader.exec_module(s40)
lengths = s40.read_fai(D + "/input/Africa_hap2.fa.fai")
cmap = s40.build_chrom_map(s40.read_fai("${ANALYSIS_DIR}/14_pan_genome/04_SNP_calling/00_renamed_genomes/Africa_hap2.fasta.fai"), lengths)
seqs = s40.read_fasta(D + "/input/Africa_hap2.fa", set(lengths))
def keys(path):
    out = set()
    for line in open(path):
        f = line.rstrip("\n").split("\t"); c = cmap.get(f[0]); ref, alt = f[3].upper(), f[4].upper(); pos = int(f[1])
        if not c or not (1 <= pos <= lengths[c]) or len(ref) != 1 or len(alt) != 1 or ref not in "ACGT" or alt not in "ACGT": continue
        b = seqs[c][pos - 1]
        if b == "N" or b != ref: continue
        out.add((c, pos, ref, alt))
    return out
cat = {"nrly_hap1": set(), "nrly_hap2": set()}
with open(R + "/05_dSNP_minimap_hap38/snp_population_catalog.tsv") as fh:
    rd = csv.reader(fh, delimiter="\t"); hdr = next(rd); iS = hdr.index("Samples")
    for f in rd:
        s = f[iS]
        if "nrly_hap" not in s: continue
        for h in cat:
            if h in s.split(";"): cat[h].add((f[1], int(f[2]), f[3], f[4]))
rare = {"nrly_hap1": set(), "nrly_hap2": set()}
for r in csv.DictReader(open(R + "/06_dSNP_phoenix_polarity/polarity/dsnp_v2_phoenix_alt_derived_candidates.tsv"), delimiter="\t"):
    if r["Samples"] in rare: rare[r["Samples"]].add((r["Chrom"], int(r["Pos"]), r["Ref"], r["Alt"]))
with open(W + "/dsnp/repro_old_vs_catalog.tsv", "w") as out:
    out.write("Haplotype\tAssembly\tRerun_SNV\tCatalog_SNV\tShared\tRerun_only\tCatalog_only\tRare_dSNP_catalog\tRare_dSNP_recovered\n")
    for h in cat:
        for asm in ("old", "final"):
            p = f"{W}/calls/nrly_{asm}_{h[-4:]}/nrly_{asm}_{h[-4:]}_vs_Africa_hap2.snp.txt"
            try: k = keys(p)
            except FileNotFoundError: continue
            out.write(f"{h}\t{asm}\t{len(k)}\t{len(cat[h])}\t{len(k & cat[h])}\t{len(k - cat[h])}\t{len(cat[h] - k)}\t{len(rare[h])}\t{len(rare[h] & k)}\n")
print(open(W + "/dsnp/repro_old_vs_catalog.tsv").read())
