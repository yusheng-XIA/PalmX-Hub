#!/bin/bash
# Population structure, diversity and LD (308 accessions; FL-Hap2 coordinates)
set -euo pipefail
threads=32
VCF=snp_308.final.vcf.gz          # 89,810,082 biallelic SNPs polymorphic among the 308 accessions

# ---- 1. Structure SNP set: per-chromosome *.struct sets (MAF >= 0.01, missing <= 10%; 01_*.sh, step 4)
#         merged (40,158,388 SNPs) and LD-pruned (50 SNPs, step 10, r2 0.2) to 3,987,046 SNPs
for i in $(seq -w 1 16); do echo chr${i}B.struct; done > merge.txt
plink --merge-list merge.txt --allow-extra-chr --keep-allele-order --make-bed --out all.missing_maf
plink --bfile all.missing_maf --indep-pairwise 50 10 0.2 --allow-extra-chr --out tmp.ld
plink --bfile all.missing_maf --extract tmp.ld.prune.in --allow-extra-chr --keep-allele-order --make-bed --out all.LDfilter

# ---- 2. PCA (PLINK v1.9, top 10 PCs) -------------------------------------------------------------------------
plink --bfile all.LDfilter --pca 10 --allow-extra-chr --out PCA_10

# ---- 3. ADMIXTURE v1.3.0, K = 2-8, 10-fold cross-validation (K = 4 used as working classification) ----------
# ADMIXTURE needs integer chromosome codes
awk 'BEGIN{OFS="\t"}{c=$1; sub(/^chr/,"",c); sub(/B$/,"",c); $1=c+0; print}' all.LDfilter.bim > all.bim
ln -sf all.LDfilter.bed all.bed; ln -sf all.LDfilter.fam all.fam
for K in $(seq 2 8); do
    admixture --cv=10 -j${threads} all.bed ${K} > admix.${K}.log 2>&1
done
grep -h "CV error" admix.*.log

# ---- 4. pi and FST: MAF >= 0.05, missing <= 20%, QUAL >= 30 (29,856,610 SNPs); 100-kb windows, 50-kb step ------
bcftools view -i 'QUAL>=30' ${VCF} -Ou | bcftools view -q 0.05:minor -Ou \
    | bcftools view -i 'F_MISSING<=0.2' -Oz -o div.vcf.gz
for pop in $(cat populations.txt); do        # one sample list per group: ${pop}.txt
    vcftools --gzvcf div.vcf.gz --keep ${pop}.txt --window-pi 100000 --window-pi-step 50000 --out pi_${pop}
done
# pi per site from the called alleles, summed per window and divided by the window length;
# genome-wide pi = unweighted mean of window values
for pair in $(cat population_pairs.txt); do        # e.g. POP1,POP2
    p1=${pair%,*}; p2=${pair#*,}
    vcftools --gzvcf div.vcf.gz --weir-fst-pop ${p1}.txt --weir-fst-pop ${p2}.txt \
        --fst-window-size 100000 --fst-window-step 50000 --out fst_${p1}_${p2}
done

# ---- 5. LD decay: div.vcf.gz thinned to >= 2 kb spacing (676,023 SNPs), PopLDdecay v3.42 ---------------------
vcftools --gzvcf div.vcf.gz --thin 2000 --recode --stdout | bgzip > ld.thin2kb.vcf.gz
for pop in $(cat populations.txt); do        # one sample list per group: ${pop}.txt
    PopLDdecay -InVCF ld.thin2kb.vcf.gz -SubPop ${pop}.txt -MaxDist 500 -MAF 0.005 -Het 0.9 -Miss 0.25 \
        -OutStat ld_${pop}
done
# Mean r2 per distance bin weighted by SNP-pair number (2-kb bins to 100 kb, 5-kb bins to 250 kb, 10-kb beyond);
# LD50 = distance (linear interpolation) at which binned mean r2 first falls below half of the first-bin value.
