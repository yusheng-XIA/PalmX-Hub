#!/bin/bash
# Reference-guided scaffolding using RagTag
# Used for dura, pisifera, and 28 pangenome accessions
threads=64

ragtag.py scaffold \
    ${reference}.fasta \
    ${sample}.contigs.fasta \
    -t ${threads} \
    -o ${sample}_ragtag_output
