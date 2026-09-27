#!/bin/bash
# Allele pairing and gene composition within each material (Fig. 4g)
#   Haplotype-resolved materials (FL, TN, TK, NS, Nigerian, E. oleifera): Hap1 vs Hap2
#   27 HiFi-only materials: hifiasm primary (p_ctg) vs alternate (a_ctg) contigs
# JCVI v1.5.11, GeneTribe v1.2.1, BEDTools v2.31.1, BLAST+ v2.16.0, GMAP 2025-07-31
set -euo pipefail
threads=32
A=${sample}_p   # or Hap1
B=${sample}_a   # or Hap2

# 0. Split genome/GFF3/CDS/protein by contig set and add material prefixes to gene IDs
for x in ${A} ${B}; do
    python -m jcvi.formats.gff bed --type=mRNA --key=ID ${x}.gff3 -o ${x}.bed
done

# 1. Reciprocal best hits between the two protein sets (GeneTribe core workflow) = candidate allele pairs
GeneTribe core -l ${A} -f ${B} -n ${threads}

# 2. Full-length gene sequences re-evaluated with BLASTN (E <= 1e-5, dust off):
#    100% identity over the whole query = sequence-identical alleles; other valid alignments = sequence-different
#    biallelic genes; pairs without an interpretable full-length alignment = unresolved
bedtools getfasta -fi ${A}.fa -bed ${A}.bed -name -s > ${A}.gene.fa
bedtools getfasta -fi ${B}.fa -bed ${B}.bed -name -s > ${B}.gene.fa
makeblastdb -in ${B}.gene.fa -dbtype nucl
blastn -query ${A}.gene.fa -db ${B}.gene.fa -evalue 1e-5 -dust no -outfmt "6 std qlen slen" \
    -num_threads ${threads} > ${A}_vs_${B}.gene.blastn

# 3. Side-specific genes: no CDS-BLAST hit AND no GMAP placement on the opposite contig set
makeblastdb -in ${B}.fa -dbtype nucl
blastn -query ${A}.cds.fa -db ${B}.fa -evalue 1e-5 -dust no -outfmt 6 -num_threads ${threads} > ${A}_cds_vs_${B}.blastn
gmap_build -D gmapdb -d ${B} ${B}.fa
gmap -D gmapdb -d ${B} -t ${threads} -f samse ${A}.cds.fa > ${A}_cds_to_${B}.gmap.sam
# (repeat steps 2-3 in the B -> A direction)

# Fig. 3b (FL and TN) uses one-to-one allele pairs from chromosomal synteny supplemented by GMAP placement:
python -m jcvi.compara.catalog ortholog ${sample}_hap1 ${sample}_hap2 --cscore=.99 --no_strip_names
