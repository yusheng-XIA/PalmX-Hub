#!/bin/bash
# Population structure, diversity, and LD analysis
threads=64

# ============================================================
# 1. LD pruning for structure analysis
# ============================================================
plink2 --vcf final_snps.vcf.gz \
    --indep-pairwise 50 10 0.3 \
    --out ld_pruned \
    --allow-extra-chr

plink2 --vcf final_snps.vcf.gz \
    --extract ld_pruned.prune.in \
    --make-bed --out ld_pruned \
    --allow-extra-chr

# ============================================================
# 2. ADMIXTURE (K = 2 to 10)
# ============================================================
for K in $(seq 2 10); do
    admixture --cv ld_pruned.bed ${K} -j${threads} \
        | tee admixture_K${K}.log
done

# Determine optimal K by cross-validation error
grep -h "CV" admixture_K*.log > cv_errors.txt

# ============================================================
# 3. PCA
# ============================================================
plink2 --bfile ld_pruned --pca 20 --out pca_results --allow-extra-chr

# ============================================================
# 4. Nucleotide diversity (pi) and Fst
# ============================================================
for pop in pop1 pop2 pop3 pop4 pop5 pop6; do
    vcftools --gzvcf final_snps.vcf.gz \
        --keep ${pop}.samples.txt \
        --window-pi 100000 \
        --out ${pop}_pi
done

# Pairwise Fst
vcftools --gzvcf final_snps.vcf.gz \
    --weir-fst-pop pop1.samples.txt \
    --weir-fst-pop pop2.samples.txt \
    --fst-window-size 100000 \
    --out pop1_vs_pop2_fst

# ============================================================
# 5. Heterozygosity and inbreeding
# ============================================================
vcftools --gzvcf final_snps.vcf.gz --het --out heterozygosity

# ============================================================
# 6. LD decay
# ============================================================
PopLDdecay -InVCF final_snps.vcf.gz \
    -OutStat ld_decay \
    -SubPop pop1.samples.txt \
    -bin1 100 -bin2 1000 -break 500
