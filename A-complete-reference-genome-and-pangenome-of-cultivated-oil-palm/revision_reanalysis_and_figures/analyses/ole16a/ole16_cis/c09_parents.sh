#!/bin/bash
F=${ANALYSIS_DIR}/22_answer_reviews/00_ms/03_V3/03_figure3/01_ASE/00_shared/legacy_expression_normalized_long.tsv
head -1 $F
grep -E "evm.TU.chr11B.1497|evm.TU.chr04B.697|evm.TU.chr10B.1088" $F
