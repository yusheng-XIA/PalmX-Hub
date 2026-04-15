#!/bin/bash
# Synteny-based scaffolding using ALLMAPS (JCVI)
# Used when Hi-C data is not available

# Step 1: Generate synteny anchors
grep "gene" ${sample}.gff3 | awk -F ";" '{print $1}' \
    | awk '{print $1"\t"$4"\t"$5"\t"$9}' | sed 's/ID=//g' > ${sample}.bed
sed 's/evm.model/evm.TU/g' ${sample}.cds > ${sample}.cds.renamed

python -m jcvi.compara.catalog ortholog --nostdpf ref ${sample}

# Step 2: Transfer syntenic coordinates
python -m jcvi.assembly.syntenypath bed ref.${sample}.anchors \
    --switch -o ${sample}.transfer.bed

# Step 3: Merge and scaffold
python -m jcvi.assembly.allmaps mergebed ${sample}.transfer.bed \
    -o ${sample}.transfer.merge.bed
python -m jcvi.assembly.allmaps path ${sample}.transfer.merge.bed \
    ${sample}.contig.fasta

# Step 4: Lift over gene annotations
liftOver -gff ${sample}.gff3 \
    ${sample}.transfer.merge.chain \
    ${sample}.genome.gff3 \
    ${sample}.genome.unmapped
