#!/bin/bash
# Resistance gene analysis and gene family evolution
threads=64

# ============================================================
# 1. RGA prediction (RGAugury)
# ============================================================
for sample in $(cat genome_list.txt); do
    RGAugury.pl \
        -p ${sample}.protein.faa \
        -n ${sample}.cds.fna \
        -g ${sample}.fasta \
        -gff ${sample}.gff3 \
        -c ${threads} \
        -pfx ${sample}_RGA
done

# ============================================================
# 2. RGA tandem gene identification
# ============================================================
for class in NBS RLK RLP TMCC; do
    makeblastdb -in ${sample}_RGA.${class}.candidates.faa -dbtype prot
    blastp -db ${sample}_RGA.${class}.candidates.faa \
        -query ${sample}_RGA.${class}.candidates.faa \
        -evalue 1e-10 -outfmt 6 -num_threads ${threads} \
        -out ${sample}.${class}.self.blastp

    python -m jcvi.compara.catalog tandem \
        ${sample}.${class}.self.blastp \
        ${sample}_RGA.${class}.candidates.faa \
        ${sample}.bed \
        -o ${sample}.${class}.tandem.list
done

# ============================================================
# 3. Cross-species RGA synteny
# ============================================================
python -m jcvi.compara.catalog ortholog --dbtype=prot ${sample1} ${sample2}

cat ${sample1}.NBS.tandem.list ${sample1}.RLK.tandem.list \
    ${sample1}.RLP.tandem.list ${sample1}.TMCC.tandem.list \
    | sed 's/,/\t/g' > ${sample1}.RGA.tandem.list

python -m jcvi.compara.synteny mcscan ${sample1}.bed \
    ${sample1}.${sample2}.lifted.anchors \
    --mergetandem=${sample1}.RGA.tandem.list \
    --iter=1 -o ${sample1}.${sample2}.blocks

# ============================================================
# 4. Gene family evolution (CAFE5)
# ============================================================
# Build species tree with OrthoFinder single-copy orthologs
# Calibrate with TimeTree divergence times

cafe5 -i gene_family_counts.tsv \
    -t species_tree.nwk \
    -p \
    -o cafe5_results
