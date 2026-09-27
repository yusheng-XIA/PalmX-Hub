#!/bin/bash
# Short-read mapping and SNP calling (308 accessions, jointly genotyped with 470 further accessions; 778 in total)
# Reference: FL-Hap2. NVIDIA Clara Parabricks v4.3.0-1, GATK v4.2.0.0, bcftools
set -euo pipefail
ref=FL_Hap2.fa

# 1. GPU-accelerated BWA-MEM mapping (fq2bamfast) and per-sample GVCF (HaplotypeCaller)
while read sample; do
    pbrun fq2bamfast --ref ${ref} --in-fq ${sample}_R1.fq.gz ${sample}_R2.fq.gz \
        --out-bam ${sample}.bam
    pbrun haplotypecaller --ref ${ref} --in-bam ${sample}.bam --gvcf --out-variants ${sample}.g.vcf.gz
done < samples_778.txt

# 2. Chromosome-wise joint genotyping
awk '{print $1"\t"$1".g.vcf.gz"}' samples_778.txt > sample_map.txt
for chr in $(seq -f "chr%02gB" 1 16); do
    gatk GenomicsDBImport --sample-name-map sample_map.txt --genomicsdb-workspace-path gdb_${chr} -L ${chr}
    gatk GenotypeGVCFs -R ${ref} -V gendb://gdb_${chr} -O raw_${chr}.vcf.gz
done
gatk GatherVcfs $(for c in $(seq -f "chr%02gB" 1 16); do echo -I raw_${c}.vcf.gz; done) -O raw_778.vcf.gz
tabix -p vcf raw_778.vcf.gz

# 3. Per-chromosome filter chain: biallelic SNPs, GATK hard filter, genotypes with GQ < 20, DP < 5 or DP > 100
#    set to missing, sites with > 25% missing genotypes removed, polymorphic sites kept
for chr in $(seq -f "chr%02gB" 1 16); do
    bcftools norm -r ${chr} -f ${ref} -m -any -Ou raw_${chr}.vcf.gz \
    | bcftools view -v snps -m2 -M2 -Ou \
    | bcftools annotate --set-id %CHROM:%POS:%REF:%FIRST_ALT -Ou \
    | bcftools filter -s SNP_HARD_FILTER \
        -e "INFO/QD<2.0 || INFO/MQ<40.0 || INFO/FS>60.0 || INFO/SOR>3.0 || INFO/MQRankSum<-12.5 || INFO/ReadPosRankSum<-8.0" -Ou \
    | bcftools view -f PASS,. -Ou \
    | bcftools filter -S . -e "FMT/GQ<20 | FMT/DP<5 | FMT/DP>100" -Ou \
    | bcftools +fill-tags -Ou -- -t AC,AN,AF,MAF,F_MISSING,NS \
    | bcftools view -i "INFO/F_MISSING<=0.25 && INFO/AC>=1 && INFO/AC<INFO/AN" -Oz -o ${chr}.diversity_snps.vcf.gz
    bcftools index -t ${chr}.diversity_snps.vcf.gz
done

# 4. 308-accession subset, polymorphic sites only (89,810,082 SNPs), and PLINK sets
for chr in $(seq -f "chr%02gB" 1 16); do
    bcftools view -S samples_308.txt -Ou ${chr}.diversity_snps.vcf.gz \
    | bcftools +fill-tags -Ou -- -t AN,AC,F_MISSING,MAF \
    | bcftools view -i "INFO/AC>0 && INFO/AC<INFO/AN" -Ob -o ${chr}.div308.bcf
    bcftools index ${chr}.div308.bcf
    # structure set (MAF >= 0.01, missing <= 10%) and GWAS set (MAF >= 0.05, missing <= 10%)
    plink --bcf ${chr}.div308.bcf --geno 0.1 --maf 0.01 --biallelic-only strict --allow-extra-chr \
        --set-missing-var-ids @:# --keep-allele-order --make-bed --out ${chr}.struct
    plink --bcf ${chr}.div308.bcf --geno 0.1 --maf 0.05 --allow-extra-chr \
        --set-missing-var-ids @:# --keep-allele-order --make-bed --out ${chr}.gwas
done
bcftools concat -Oz -o snp_308.final.vcf.gz $(seq -f "chr%02gB.div308.bcf" 1 16) && tabix -p vcf snp_308.final.vcf.gz
