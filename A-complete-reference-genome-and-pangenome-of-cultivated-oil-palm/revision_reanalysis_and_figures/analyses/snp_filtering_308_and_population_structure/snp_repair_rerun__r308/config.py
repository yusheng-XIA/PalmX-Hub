"""Rerun of the SNP-GWAS on the repaired joint-called SNP call set (308-accession subset; --geno 0.1 --maf 0.05)."""
from pathlib import Path
W = Path("${CLUSTER_WORK}/snp_repair_r308/gwas")
RUN = Path("${ANALYSIS_DIR}/22_answer_reviews/00_ms/03_V3/05_figure/07_GWAS_zero_extreme_excluded_20260723")
GENO = Path("${CLUSTER_WORK}/snp_repair_r308/gwas/geno")
BED, BIM, FAM = GENO / "GP_new_o.bed", GENO / "GP_new_o.bim", GENO / "GP_new_o.fam"
KIN_SNP = GENO / "GP_new_o.aBN.kinf"
SVK = Path("${ANALYSIS_DIR}/14_pan_genome/06_Minigraph/Pangenie/02_sv_combined/步骤七_SV_GWAS/kinship")
SV_TPED, SV_TFAM, KIN_SV = SVK / "sv_assoc.tped", SVK / "sv_assoc.tfam", SVK / "sv.kinf"
SV_COV = Path("${ANALYSIS_DIR}/14_pan_genome/06_Minigraph/Pangenie/02_sv_combined/步骤五_群体结构/pca/sv_covariates_5PC_with_intercept.txt")
SNP_PCA = Path("${CLUSTER_WORK}/snp_repair_r308/struct/PCA_10.eigenvec")
Q4, QFAM = Path("${CLUSTER_WORK}/enh_gwas/in/all.4.Q"), Path("${CLUSTER_WORK}/enh_gwas/in/all.fam")
N_SNP_TESTS = int(open(GENO / "n_tests.txt").read()) if (GENO / "n_tests.txt").exists() else 0
BONF_SNP = 0.05 / N_SNP_TESTS if N_SNP_TESTS else None
BONF_SV = 0.05 / 370706
CHROMS = [f"chr{i:02d}B" for i in range(1, 17)]
FOCAL = ["C12_0_Lauric_acid", "C14_0_Myristic_acid", "Flesh_thickness_mm", "C18_3n3_Alpha_linolenic_acid"]
VALID = ["Nut_weight_g", "C12_0_Lauric_acid", "Flesh_thickness_mm"]
