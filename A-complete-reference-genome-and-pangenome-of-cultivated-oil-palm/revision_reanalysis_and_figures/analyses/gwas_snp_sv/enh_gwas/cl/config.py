"""Paths shared by the enh_gwas scripts (read-only inputs belong to user)."""
from pathlib import Path
W = Path("${CLUSTER_WORK}/enh_gwas")
RUN = Path("${ANALYSIS_DIR}/22_answer_reviews/00_ms/03_V3/05_figure/07_GWAS_zero_extreme_excluded_20260723")
GENO = Path("${ANALYSIS_DIR}/05_GWAS/00_analysis/06_GWAS/shared_genotype")
BED, BIM, FAM = GENO / "GP_maf_allChr.bed", GENO / "GP_maf_allChr.bim", GENO / "GP_maf_allChr.fam"
KIN_SNP = GENO / "GP_maf_allChr.hBN.kinf"
SVK = Path("${ANALYSIS_DIR}/14_pan_genome/06_Minigraph/Pangenie/02_sv_combined/步骤七_SV_GWAS/kinship")
SV_TPED, SV_TFAM, KIN_SV = SVK / "sv_assoc.tped", SVK / "sv_assoc.tfam", SVK / "sv.kinf"
SV_COV = Path("${ANALYSIS_DIR}/14_pan_genome/06_Minigraph/Pangenie/02_sv_combined/步骤五_群体结构/pca/sv_covariates_5PC_with_intercept.txt")
SNP_PCA = W / "in/PCA_out.eigenvec"      # PLINK 1.9 --pca 10 on the 5,883,108 LD-pruned SNPs (Fig. 4b)
Q4, QFAM = W / "in/all.4.Q", W / "in/all.fam"
N_SNP_TESTS = 28192651
BONF_SNP = 0.05 / N_SNP_TESTS
BONF_SV = 0.05 / 370706
CHROMS = [f"chr{i:02d}B" for i in range(1, 17)]
FOCAL = ["C12_0_Lauric_acid", "C14_0_Myristic_acid", "Flesh_thickness_mm", "C18_3n3_Alpha_linolenic_acid"]
VALID = ["Nut_weight_g", "C12_0_Lauric_acid", "Flesh_thickness_mm"]
