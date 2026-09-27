#!/bin/bash
# SNP- and SV-GWAS with EMMAX (release 2012-02-10; emmax-intel64 -v -d 10)
# SNPs : 308 accessions, biallelic SNPs with MAF >= 0.05 and missing rate <= 10%; SNP kinship
# SVs  : 370,136 PanGenie-genotyped SVs (MAF >= 0.05, missing <= 0.1); SV kinship + intercept + SV PC1-5
set -euo pipefail
G=gwas/genotype

# ---- genotype preparation --------------------------------------------------------------
# per-chromosome GWAS sets from 04_population_genomics/01_mapping_and_variant_calling.sh (MAF >= 0.05, missing <= 10%)
for i in $(seq -w 1 16); do echo chr${i}B.gwas; done > gwas_merge.txt
plink --merge-list gwas_merge.txt --allow-extra-chr --keep-allele-order --make-bed --out ${G}/snp
for c in $(seq -f "chr%02gB" 1 16); do
    plink --bfile ${G}/snp --chr ${c} --allow-extra-chr --recode12 transpose --output-missing-genotype 0 \
        --out ${G}/snp_${c}
done
plink --bfile ${G}/snp --allow-extra-chr --recode12 transpose --output-missing-genotype 0 --out ${G}/snp_all
emmax-kin-intel64 -v -h -d 10 ${G}/snp_all          # -> snp_all.hBN.kinf (Balding-Nichols kinship)
mv ${G}/snp_all.hBN.kinf ${G}/snp.hBN.kinf

plink --vcf pangenie_sv_308.maf05_miss10.vcf.gz --allow-extra-chr --double-id --make-bed --out ${G}/sv
plink --bfile ${G}/sv --allow-extra-chr --recode12 transpose --output-missing-genotype 0 --out ${G}/sv_assoc
emmax-kin-intel64 -v -h -d 10 ${G}/sv_assoc && mv ${G}/sv_assoc.hBN.kinf ${G}/sv.kinf
plink --bfile ${G}/sv --allow-extra-chr --pca 5 --out ${G}/sv_pca
awk '{printf "%s %s 1", $1, $2; for (i = 3; i <= 7; i++) printf " %s", $i; print ""}' ${G}/sv_pca.eigenvec \
    > ${G}/sv_covariates_5PC_with_intercept.txt

# ---- association scans (one run per retained trait) --------------------------------------
tail -n +2 gwas/manifests/sv_tasks.tsv | while IFS=$'\t' read -r id category trait pheno; do
    out=gwas/snp/results/${category}/${trait}; mkdir -p ${out}
    for c in $(seq -f "chr%02gB" 1 16); do
        emmax-intel64 -v -d 10 -t ${G}/snp_${c} -p ${pheno} -k ${G}/snp.hBN.kinf -o ${out}/emmax_${c}
    done
    out=gwas/sv/results/${category}/${trait}; mkdir -p ${out}
    emmax-intel64 -v -d 10 -t ${G}/sv_assoc -p ${pheno} -k ${G}/sv.kinf \
        -c ${G}/sv_covariates_5PC_with_intercept.txt -o ${out}/${trait}
done

# Bonferroni: SNP 0.05 / number of SNP tests; SV 0.05 / 370,136 = 1.3509e-7.
# Significant variants are merged into reporting intervals with a maximum span of 250 kb (03d_summarize_scans.py).
