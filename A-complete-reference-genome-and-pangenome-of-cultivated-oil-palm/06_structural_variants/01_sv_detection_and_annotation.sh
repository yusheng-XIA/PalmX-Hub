#!/bin/bash
# Structural variant detection using multiple callers + SURVIVOR merging
threads=64

# ============================================================
# 1. Read-based SV calling: pbsv (HiFi reads)
# ============================================================
pbmm2 align ${reference} ${sample}.hifi.fq.gz ${sample}.pbmm2.bam \
    --sort --preset CCS --sample ${sample} -j 5 -J 5
pbsv discover ${sample}.pbmm2.bam ${sample}.svsig.gz \
    -s ${sample} --tandem-repeats ${reference}.trf.bed
pbsv call ${reference} ${sample}.svsig.gz ${sample}_pbsv.vcf \
    --min-sv-length 30 --max-ins-length 100K --max-dup-length 100K \
    -j ${threads} --ccs

awk -F '[;\t]' '{if ($1~/^#/) print; else if ($8=="SVTYPE=INS") print}' ${sample}_pbsv.vcf > ${sample}_pbsv_INS.vcf
awk -F '[;\t]' '{if ($1~/^#/) print; else if ($8=="SVTYPE=DEL") print}' ${sample}_pbsv.vcf > ${sample}_pbsv_DEL.vcf

# ============================================================
# 2. Read-based SV calling: cuteSV
# ============================================================
# HiFi reads
minimap2 -a -x map-hifi -t ${threads} ${reference} ${sample}.hifi.fq.gz --MD \
    | samtools sort -@ ${threads} -o ${sample}.hifi.sorted.bam
samtools index ${sample}.hifi.sorted.bam
cuteSV ${sample}.hifi.sorted.bam ${reference} ${sample}.cuteSV.vcf ./ \
    --sample ${sample} --min_size 30 --max_size 100000 \
    --max_cluster_bias_INS 1000 --diff_ratio_merging_INS 0.9 \
    --max_cluster_bias_DEL 1000 --diff_ratio_merging_DEL 0.5 \
    --min_support 3 -t ${threads}

# ONT reads
minimap2 -a -x map-ont -t ${threads} ${reference} ${sample}.ont.fq.gz --MD \
    | samtools sort -@ ${threads} -o ${sample}.ont.sorted.bam
samtools index ${sample}.ont.sorted.bam
cuteSV ${sample}.ont.sorted.bam ${reference} ${sample}.ont.cuteSV.vcf ./ \
    --sample ${sample} --min_size 30 --max_size 100000 \
    --max_cluster_bias_INS 100 --diff_ratio_merging_INS 0.3 \
    --max_cluster_bias_DEL 100 --diff_ratio_merging_DEL 0.3 \
    --min_support 10 -t ${threads}

awk -F '[;\t]' '{if ($1~/^#/) print; else if ($9=="SVTYPE=INS") {gsub(/END=[0-9]+/,"END="$2); print $1"\t"$2-1"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10}}' ${sample}.cuteSV.vcf > ${sample}.cuteSV_INS.vcf
awk -F '[;\t]' '{if ($1~/^#/) print; else if ($9=="SVTYPE=DEL") {gsub(/END=[0-9]+/,"END="$2); print $1"\t"$2-1"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10}}' ${sample}.cuteSV.vcf > ${sample}.cuteSV_DEL.vcf

# ============================================================
# 3. Assembly-based SV calling: SVIM-asm
# ============================================================
minimap2 -a -x asm5 --cs -r 2k -t ${threads} ${reference} ${sample}.contigs.fasta \
    | samtools sort -@ ${threads} -o ${sample}.asm.sorted.bam
samtools index ${sample}.asm.sorted.bam
svim-asm haploid ./ ${sample}.asm.sorted.bam ${reference} \
    --min_sv_size 30 --max_sv_size 100000 --sample ${sample}

awk -F '[;\t]' '{if ($1~/^#/) print; else if ($8=="SVTYPE=INS") print}' variants.vcf > ${sample}.SVIM_asm_INS.vcf
awk -F '[;\t]' '{if ($1~/^#/) print; else if ($8=="SVTYPE=DEL") print}' variants.vcf > ${sample}.SVIM_asm_DEL.vcf

# ============================================================
# 4. Assembly-based SV calling: SyRI
# ============================================================
nucmer -p ${sample}_vs_ref -t ${threads} ${reference} ${sample}.fasta
delta-filter -m ${sample}_vs_ref.delta > ${sample}_vs_ref.filtered.delta
show-coords -THrd ${sample}_vs_ref.filtered.delta > ${sample}_vs_ref.filtered.coords
syri -c ${sample}_vs_ref.filtered.coords \
    -d ${sample}_vs_ref.filtered.delta \
    -r ${reference} -q ${sample}.fasta \
    --prefix ${sample}.syri -s show-snps

# Inversions
grep "<INV>" ${sample}.syri.vcf | awk '{print $0"\tGT\t1/1"}' | bgzip > ${sample}.syri.INV.vcf.gz
tabix -p vcf ${sample}.syri.INV.vcf.gz

# Translocations
grep -P "<TRANS>|<INVTR>" ${sample}.syri.vcf | awk '{print $0"\tGT\t1/1"}' | bgzip > ${sample}.syri.TRANS.vcf.gz
tabix -p vcf ${sample}.syri.TRANS.vcf.gz

# Indels from SyRI
grep -v "#" ${sample}.syri.vcf | grep -P "INS|DEL" | grep -v -P "<INS|<DEL>" \
    | awk '{if ((length($5)-length($4))*(length($5)-length($4))<900) print $1"\t"$2"\t"tolower(substr($3,1,3))"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\tGT\t1/1"}' \
    > ${sample}.syri.InDel.vcf

# Visualization
plotsr --sr ${sample}.syri.out --genomes genomes.txt -H 8 -W 5 -o ${sample}_synteny.png

# ============================================================
# 5. SURVIVOR merge
# ============================================================
# INS: HiFi 3 callers, require 2/3 support, merge distance 50bp
echo -e "${sample}_pbsv_INS.vcf\n${sample}.cuteSV_INS.vcf\n${sample}.SVIM_asm_INS.vcf" > ${sample}_INS.files
SURVIVOR merge ${sample}_INS.files 50 2 1 1 0 30 ${sample}_INS.vcf

# DEL: HiFi 3 callers, merge at 50% reciprocal overlap
echo -e "${sample}_pbsv_DEL.vcf\n${sample}.cuteSV_DEL.vcf\n${sample}.SVIM_asm_DEL.vcf" > ${sample}_DEL.files
SURVIVOR merge ${sample}_DEL.files 0.5 2 1 1 0 30 ${sample}_DEL.vcf

# Normalize and sort
awk '/^##/{print} /^#CHROM/{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10} !/^#/{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\tGT\t1/1"}' \
    ${sample}_INS.vcf | bcftools norm -c x -f ${reference} -d all -Ov --threads ${threads} \
    | bcftools sort -Ov -o ${sample}_INS_norm.vcf

# ============================================================
# 6. SNP calling: longshot (long reads) + show-snps (assembly)
# ============================================================
# longshot from HiFi
longshot --bam ${sample}.hifi.sorted.bam -c 3 -D 3:10:50 \
    --ref ${reference} --out ${sample}.hifi.longshot.vcf

# longshot from ONT
longshot --bam ${sample}.ont.sorted.bam -c 10 -D 3:10:50 \
    --ref ${reference} --out ${sample}.ont.longshot.vcf

# show-snps from assembly alignment
nucmer -p ${sample}_vs_ref -t ${threads} ${reference} ${sample}.fasta
delta-filter -1 ${sample}_vs_ref.delta > ${sample}_vs_ref.1to1.delta
show-snps -C -T ${sample}_vs_ref.1to1.delta > ${sample}.snps.txt
show-snps -C -I -T ${sample}_vs_ref.1to1.delta > ${sample}.snps_only.txt

# ============================================================
# 7. TE-mediated SV analysis
# ============================================================
for svtype in INS DEL INV DUP; do
    for flank in up dn; do
        for repeat_type in TE LTR TIR HEL SD TR ALL; do
            bedtools coverage \
                -a ${sample}.${svtype}.${flank}.100.bed \
                -b ${sample}.${repeat_type}.bed \
                -wo > ${sample}.${svtype}.${flank}.${repeat_type}.coverage.bed
        done
    done
    paste ${sample}.${svtype}.up.*.coverage.bed ${sample}.${svtype}.dn.*.coverage.bed \
        > ${sample}.${svtype}.all.coverage.bed
done

# Simulation control
bedtools shuffle -i ${sample}.${svtype}.bed \
    -g ${sample}.genome.sizes \
    -excl ${sample}.svall.bed \
    > ${sample}.${svtype}.shuffled.bed

# ============================================================
# 8. SV genomic context and evolutionary age
# ============================================================
# Genomic context annotation
bedtools intersect -a merged_svs.bed -b ${reference}.genes.bed -wo > svs_genic.bed

# TE overlap (>=50% length)
bedtools intersect -a merged_svs.bed -b ${sample}.te.bed -wo -f 0.5 > svs_te_overlap.bed

# Evolutionary age: ancient (all 3 groups) / derived non-specific (2 groups) / derived specific (1 group)
# SV hotspots: 500-kb windows > 95th percentile density
bedtools makewindows -g ${reference}.genome.sizes -w 500000 > windows_500kb.bed
bedtools intersect -a windows_500kb.bed -b merged_svs.bed -c > sv_density_500kb.bed