#!/bin/bash
# Single-nucleus RNA-seq processing with Cell Ranger
# FL and TN varieties at 95d, 125d, 185d

# ============================================================
# Cell Ranger alignment
# ============================================================
for sample in FL_95d FL_125d FL_185d TN_95d TN_125d TN_185d; do
    cellranger count \
        --id=${sample} \
        --transcriptome=${reference}_cellranger_ref \
        --fastqs=${sample}_fastqs/ \
        --sample=${sample} \
        --localcores=${threads} \
        --localmem=128 \
        --expect-cells=10000
done
