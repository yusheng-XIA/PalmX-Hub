#!/bin/bash
# Allele-aware Hi-C scaffolding of the phased contigs (TN, TK, NS, Nigerian, E. oleifera)
# with HapHiC, followed by manual curation in Juicebox (v1.9.8)
set -euo pipefail
threads=64
nchrom=16

# 1. Hi-C read mapping and filtering (HapHiC recommended workflow)
cat ${sample}.hap1.p_ctg.fa ${sample}.hap2.p_ctg.fa > ${sample}.contigs.fa
bwa index ${sample}.contigs.fa
bwa mem -5SP -t ${threads} ${sample}.contigs.fa ${sample}_HiC_R1.fq.gz ${sample}_HiC_R2.fq.gz \
    | samblaster | samtools view - -@ ${threads} -S -h -b -F 3340 -o ${sample}.HiC.bam
filter_bam ${sample}.HiC.bam 1 --nm 3 --threads ${threads} | samtools view - -b -@ ${threads} -o ${sample}.HiC.filtered.bam

# 2. Allele-aware scaffolding (2 x 16 pseudochromosomes for a phased diploid)
haphic pipeline ${sample}.contigs.fa ${sample}.HiC.filtered.bam $((nchrom * 2)) --threads ${threads}

# 3. Juicebox review files; after manual curation the reviewed assembly is converted back to FASTA
cd 04.build
bash juicebox.sh
# ... curate out_JBAT.hic / out_JBAT.assembly in Juicebox, export out_JBAT.review.assembly ...
juicer post -o out_JBAT out_JBAT.review.assembly out_JBAT.liftover.agp ../${sample}.contigs.fa
