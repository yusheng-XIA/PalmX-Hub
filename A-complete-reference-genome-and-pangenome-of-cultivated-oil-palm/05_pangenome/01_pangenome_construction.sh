#!/bin/bash
# Pangenome: gene family clustering + graph construction
threads=64

# ============================================================
# 1. Gene family clustering (OrthoFinder)
# ============================================================
orthofinder -f proteomes_dir/ \
    -t ${threads} \
    -a ${threads} \
    -S diamond \
    -I 1.5 \
    -o orthofinder_results

# Classify by occupancy (34 genomes)
# core: 34/34, soft-core: >=33/34, shell: 2-32/34, cloud: 1/34

# ============================================================
# 2. Graph-based pangenome (minigraph-cactus)
# ============================================================
# Step 1: minigraph primary graph
minigraph -xggs -t ${threads} ${reference}.fasta \
    $(cat other_assemblies.list) > pangenome.gfa

# Step 2: Cactus refinement
cactus-pangenome jobstore \
    seqfile.txt \
    --outDir pangenome_out \
    --outName oil_palm_pangenome \
    --reference ${reference_name} \
    --vcf --giraffe \
    --maxCores ${threads}
