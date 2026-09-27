"""Input locations shared by the GWAS scripts. Set GWAS_DIR (analysis directory) and GWAS_WORK (output
directory) or edit the paths below.

GWAS_DIR layout
  manifests/sv_tasks.tsv                  category, trait (the 60 retained traits)
  phenotypes/<category>/<trait>.txt       FID IID value (zeros and 2.5%/97.5% tails set to NA; 01_phenotype_filter.py)
  tables/all_trait_gwas_loci.tsv          Bonferroni reporting intervals of the SNP and SV scans
  snp/results/<category>/<trait>/emmax_<chrom>.ps / .reml   EMMAX output (02_run_emmax.sh)
  genotype/snp.{bed,bim,fam}, genotype/snp.hBN.kinf        SNP genotypes (MAF >= 0.05, missing <= 10%) and kinship
  genotype/sv_assoc.{tped,tfam}, genotype/sv.kinf          PanGenie SV genotypes and SV kinship
  genotype/sv_covariates_5PC_with_intercept.txt            intercept + SV PC1-5
  genotype/snp_pca.eigenvec, genotype/all.4.Q, genotype/all.fam   PLINK PCA of LD-pruned SNPs; ADMIXTURE K = 4
"""
import os
from pathlib import Path

RUN = Path(os.environ.get("GWAS_DIR", "gwas"))
W = Path(os.environ.get("GWAS_WORK", "gwas_work"))
GENO = RUN / "genotype"
BED, BIM, FAM = GENO / "snp.bed", GENO / "snp.bim", GENO / "snp.fam"
KIN_SNP = GENO / "snp.hBN.kinf"
SV_TPED, SV_TFAM, KIN_SV = GENO / "sv_assoc.tped", GENO / "sv_assoc.tfam", GENO / "sv.kinf"
SV_COV = GENO / "sv_covariates_5PC_with_intercept.txt"
SV_META = GENO / "sv_qc.per_sv.tsv"          # per-SV id, svtype, svlen, ref_len, alt_len
SNP_PCA = GENO / "snp_pca.eigenvec"          # PLINK 1.9 --pca 10 on the LD-pruned SNP set
Q4, QFAM = GENO / "all.4.Q", GENO / "all.fam"
N_SNP_TESTS = sum(1 for _ in open(BIM)) if BIM.exists() else None   # markers in the bed (row-count checks)
N_SNP_BONF = 25928923                        # SNPs with MAF >= 0.05 and missing rate <= 10% among the 308 accessions
BONF_SNP = 0.05 / N_SNP_BONF                 # 1.9283e-9
N_SV_UNIQUE = 370136                         # 370,706 tests on 370,136 unique SVs
BONF_SV = 0.05 / N_SV_UNIQUE                 # 1.3509e-7
CHROMS = [f"chr{i:02d}B" for i in range(1, 17)]
FOCAL = ["C12_0_Lauric_acid", "C14_0_Myristic_acid", "Flesh_thickness_mm", "C18_3n3_Alpha_linolenic_acid"]
VALID = ["Nut_weight_g", "C12_0_Lauric_acid", "Flesh_thickness_mm"]
