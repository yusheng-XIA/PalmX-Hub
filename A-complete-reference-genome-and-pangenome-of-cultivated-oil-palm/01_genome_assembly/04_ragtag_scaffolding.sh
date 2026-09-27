#!/bin/bash
# Reference-guided scaffolding with RagTag
# (a) TN, TK, NS, Nigerian, E. oleifera: only where chromosome order or orientation needed additional support
# (b) 27 HiFi-only accessions: primary contigs scaffolded separately against FL-Hap2 and FL-Hap1; the two
#     results were integrated and, after auditing continuity, completeness, chromosome naming and HiFi-read
#     remapping, one representative chromosome-scale assembly was retained per accession.
set -euo pipefail
threads=32

for ref in FL_Hap2 FL_Hap1; do
    ragtag.py scaffold -t ${threads} -o ${sample}_ragtag_${ref} ${ref}.fa ${sample}.p_ctg.fa
done
# The FL-Hap2- and FL-Hap1-guided scaffolds (ragtag.scaffold.agp / ragtag.scaffold.fasta) were compared
# chromosome by chromosome and integrated into ${sample}.chr.fa (one representative assembly per accession).

# Remapping check of the retained assembly
minimap2 -ax map-hifi -t ${threads} ${sample}.chr.fa ${sample}.hifi.fq.gz \
    | samtools sort -@ ${threads} -o ${sample}.hifi_remap.bam
samtools flagstat ${sample}.hifi_remap.bam > ${sample}.hifi_remap.flagstat
