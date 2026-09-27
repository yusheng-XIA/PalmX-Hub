#!/bin/bash
# Tandem repeats (Tandem Repeats Finder v4.09.1) and segmental duplications (BISER v1.4)
set -euo pipefail
threads=32

trf ${sample}.fa 2 6 6 80 10 50 2000 -d -h
biser --keep-contigs -t ${threads} -o ${sample}.biser.bedpe ${sample}.fa
