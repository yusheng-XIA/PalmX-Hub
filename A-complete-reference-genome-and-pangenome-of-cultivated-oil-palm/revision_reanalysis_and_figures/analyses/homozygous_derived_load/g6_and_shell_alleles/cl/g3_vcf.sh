E=${CLUSTER_HOME}/miniconda3/envs/cs213/bin
IN=${JOINT_SNP_DIR}/06_filter/oil_palm_joint.population_snps.corrected.vcf.gz
O=${CLUSTER_WORK}/headline_pop/out
$E/bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%QUAL[\t%GT]\n' -r chr01B:3259080-3259300,chr01B:3153030 $IN > $O/g3_shell_exon1.tsv
$E/bcftools query -l $IN > $O/g3_samples_joint.txt
cut -f1-5 $O/g3_shell_exon1.tsv
grep -c . $O/g3_samples_joint.txt
awk '$1=="chr01B" && $4>=3237000 && $4<=3260000' ${ANALYSIS_DIR}/05_GWAS/00_analysis/06_GWAS/shared_genotype/GP_maf_allChr.bim | wc -l
awk '$1=="chr01B" && $4>=3259080 && $4<=3259300' ${ANALYSIS_DIR}/05_GWAS/00_analysis/06_GWAS/shared_genotype/GP_maf_allChr.bim
head -2 ${ANALYSIS_DIR}/05_GWAS/00_analysis/06_GWAS/shared_genotype/GP_maf_allChr.bim
