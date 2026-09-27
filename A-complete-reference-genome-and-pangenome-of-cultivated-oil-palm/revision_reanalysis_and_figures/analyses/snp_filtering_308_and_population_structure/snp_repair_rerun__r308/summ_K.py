"""Genome-wide pi and FST for the K = K groups exactly as fig4c/cl/s2_final.py ('present' set): pi = mean of window
PI_MISSING_AWARE; FST = window MEAN_FST weighted by N_VARIANTS (WEIGHTED_FST reported for reference)."""
import sys, glob, os, pandas as pd
K = int(sys.argv[1]); N = f"${CLUSTER_WORK}/snp_repair_r308/fig4c_K{K}" + (sys.argv[2] if len(sys.argv) > 2 else "")
rows = []
for p in range(1, K + 1):
    d = pd.concat([pd.read_csv(f, sep="\t") for f in glob.glob(f"{N}/out/pima/*__K{K}_Pop{p}_100kb.pima.tsv")])
    rows.append(dict(metric="pi", group_or_pair=f"K{K}_Pop{p}", value=d.PI_MISSING_AWARE.mean(), value_vcftools=d.PI_VCFTOOLS.mean(), n_windows=len(d)))
for a in range(1, K + 1):
    for b in range(a + 1, K + 1):
        d = pd.concat([pd.read_csv(f, sep="\t") for f in glob.glob(f"{N}/out/fst/*__K{K}_Pop{a}_K{K}_Pop{b}_100kb_fst.windowed.weir.fst")])
        d = d[pd.to_numeric(d.MEAN_FST, errors="coerce").notna()]
        w = d.N_VARIANTS
        rows.append(dict(metric="fst", group_or_pair=f"K{K}_Pop{a}_K{K}_Pop{b}", value=(d.MEAN_FST * w).sum() / w.sum(),
                         ref_weighted_FST=(d.WEIGHTED_FST * w).sum() / w.sum(), n_windows=len(d)))
out = pd.DataFrame(rows); out.to_csv(f"{N}/pi_fst_K{K}.tsv", sep="\t", index=False); print(out.to_string())
