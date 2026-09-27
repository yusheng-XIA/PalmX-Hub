set -e
IMG=${ANALYSIS_DIR}/05_GWAS/GWAS_tur/software/software/Reseq_genek.sif
cd ${CLUSTER_WORK}/snp_repair_rerun/gwas/geno
awk '{print $1, $2}' ${ANALYSIS_DIR}/05_GWAS/00_analysis/06_GWAS/shared_genotype/GP_maf_allChr.fam > order.txt
singularity exec -B ${DATA_ROOT} $IMG plink --bfile GP_new --indiv-sort f order.txt --allow-extra-chr --keep-allele-order --make-bed --out GP_new_o > /dev/null
python - <<'P'
import numpy as np
old=[l.split()[1] for l in open('GP_new.tfam')]; new=[l.split()[1] for l in open('GP_new_o.fam')]
ref=[l.split()[1] for l in open('order.txt')]; assert new==ref, 'order'
K=np.loadtxt('GP_new.aBN.kinf'); idx=[old.index(s) for s in new]
np.savetxt('GP_new_o.aBN.kinf', K[np.ix_(idx,idx)], fmt='%.10f', delimiter='\t')
print('ok', K.shape)
P
echo DONE > reorder.done
