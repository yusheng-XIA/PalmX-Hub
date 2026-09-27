#!/bin/bash
# FL allele-specific mapping control (Extended Data Fig. 5)
# FL whole-genome short reads (MGI, paired-end) are aligned separately to FL-Hap2 (allele A) and FL-Hap1 (allele B);
# allele counts are taken at the ASE diagnostic sites (FL-Hap1 coordinates from 05a_lift_sites_to_hap1.py).
set -euo pipefail
threads=32

# 1. Separate alignments to the two haplotypes
for hap in FL_Hap2 FL_Hap1; do
    bwa mem -t ${threads} ${hap}.fa FL_WGS_R1.fq.gz FL_WGS_R2.fq.gz \
        | samtools sort -@ ${threads} -o FL_WGS.${hap}.bam
    samtools index FL_WGS.${hap}.bam
done

# 2. Site lists (FL-Hap2 coordinates; FL-Hap1 coordinates lifted from 201-bp flanks)
python 05a_lift_sites_to_hap1.py FL.marker_qc.tsv.gz FL_Hap2.fa FL_Hap1.fa sites
awk 'BEGIN{OFS="\t"}{print $1,$2-1,$2}' sites/FL_A_positions.tsv > sites/FL_A.bed
awk 'BEGIN{OFS="\t"}{print $1,$2-1,$2}' sites/FL_B_positions.tsv > sites/FL_B.bed

# 3. Base counts at the sites (mapping and base quality >= 20)
count_bases() {  # bam side
    samtools view -u -M -L sites/FL_$2.bed $1 \
    | samtools mpileup -B -d 100000 -q 20 -Q 20 --ff UNMAP,SECONDARY,QCFAIL,DUP,SUPPLEMENTARY \
        --no-output-ins --no-output-del --no-output-ends -l sites/FL_$2_positions.tsv - \
    | python parse_mpileup.py | gzip > pileup/$2.short_read.counts.tsv.gz
}
mkdir -p pileup
count_bases FL_WGS.FL_Hap2.bam A
count_bases FL_WGS.FL_Hap1.bam B

# 4. Gene-level DNA log2(A/B) and DNA-corrected ASE
python 05c_dna_corrected_ase.py --sites sites/FL_sites_liftB.tsv.gz --pileup-dir pileup \
    --ase FL.gene_stage_ASE.tsv.gz --out-prefix FL_DNA_control
