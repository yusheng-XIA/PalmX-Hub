#!/bin/bash
# DIA proteomics (Thermo Astral): FL and TN, 19 stages x 3 biological replicates = 114 raw files
# DIA-NN v2.2.0 and directLFQ v0.3.3
set -euo pipefail
threads=32
FASTA=FL_TN_unified_exact_sequence_nr.fasta     # non-redundant FL-Hap1/Hap2 + TN-Hap1/Hap2 proteins

# 1. In silico predicted spectral library
diann-linux --fasta ${FASTA} --fasta-search --predictor --gen-spec-lib \
    --out-lib FL_TN.predicted.speclib --threads 16

# 2. First pass over all 114 runs -> refined library
diann-linux $(sed 's/^/--f /' raw_files_114.txt) --lib FL_TN.predicted.speclib --fasta ${FASTA} \
    --out first_pass/report.parquet --out-lib first_pass/current_114_refined --gen-spec-lib \
    --temp first_pass/temp --threads ${threads} --qvalue 0.01 --matrices

# 3. Final pass with the refined library (precursor q-value 0.01)
diann-linux $(sed 's/^/--f /' raw_files_114.txt) --lib first_pass/current_114_refined.parquet --fasta ${FASTA} \
    --out final_pass/report.parquet --temp final_pass/temp --threads ${threads} --qvalue 0.01 --matrices

# 4. Protein-group quantification with directLFQ from the precursor table (global 1% FDR)
directlfq lfq --input_file current114_final_global1pct_directlfq_input.tsv \
    --input_type_to_use diann_precursors --min_nonan 1 --filename_suffix primary --num_cores ${threads}

# OLE16a / OLE16b: peptides unique to one locus across FL-Hap1, FL-Hap2, TN-Hap1 and TN-Hap2; FL vs TN compared
# on peptides shared by all alleles of a locus (precursor quantities, no imputation). Oleosin, LDAP and caleosin
# abundances are sums over the protein groups assigned to each family (oleosins: hmmsearch PF01277).
