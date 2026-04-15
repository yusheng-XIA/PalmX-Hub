#!/bin/bash
# Functional annotation using InterProScan and eggNOG-mapper
threads=64
# eggNOG-mapper
emapper.py \
    -i ${sample}.protein.faa \
    --output ${sample}_eggnog \
    --cpu ${threads} \
    -m diamond

# BUSCO assessment of annotation quality
busco -m protein \
    -i ${sample}.protein.faa \
    -l embryophyta_odb10 \
    -o ${sample}_protein_busco \
    -c ${threads}
