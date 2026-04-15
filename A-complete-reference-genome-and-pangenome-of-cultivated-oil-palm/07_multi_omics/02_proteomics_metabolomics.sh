#!/bin/bash
# DIA proteomics and untargeted metabolomics processing

# ============================================================
# 1. DIA proteomics (DIA-NN)
# ============================================================
diann --f raw_files_list.txt \
    --lib "" \
    --threads ${threads} \
    --fasta ${reference}.protein.faa \
    --out proteomics_report.tsv \
    --out-lib spectral_library.tsv \
    --gen-spec-lib \
    --predictor \
    --qvalue 0.01 \
    --matrices \
    --smart-profiling \
    --peak-center \
    --no-ifs-removal

# MaxLFQ normalization with match-between-runs
# Filter: >=2 unique peptides, >=50% valid values across 76 samples

# ============================================================
# 2. Untargeted metabolomics (UPLC-QTOF-MS)
# ============================================================
# Peak detection, alignment, and annotation: Progenesis QI (Waters)
# Export: peak intensities, retention times, m/z values
# Log2-transform metabolite intensities for downstream analysis
