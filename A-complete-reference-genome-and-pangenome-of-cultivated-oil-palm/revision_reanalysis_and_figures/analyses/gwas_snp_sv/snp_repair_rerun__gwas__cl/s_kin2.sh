set -e
IMG=${ANALYSIS_DIR}/05_GWAS/GWAS_tur/software/software/Reseq_genek.sif
B=${ANALYSIS_DIR}/05_GWAS/00_analysis/06_GWAS/shared_genotype
W=${CLUSTER_WORK}/snp_repair_rerun/gwas
cd $W/kincal
# calibration: container emmax-kin (-x: chromosome names are chr01B.. and would otherwise be skipped as non-autosomal)
singularity exec -B ${DATA_ROOT} $IMG /opt/emmax-20120210/emmax-kin-intel64 -v -x -d 10 -o old_aBN_x.kinf $B/GP_maf_allChr > kin_old_x.log 2>&1
python - <<'P' > calib.txt
import numpy as np
a=np.loadtxt('old_aBN_x.kinf'); b=np.loadtxt('${ANALYSIS_DIR}/05_GWAS/00_analysis/06_GWAS/shared_genotype/GP_maf_allChr.hBN.kinf')
print('maxabsdiff',np.abs(a-b).max(),'corr',np.corrcoef(a.ravel(),b.ravel())[0,1])
P
cd $W/geno
singularity exec -B ${DATA_ROOT} $IMG plink --bfile GP_new_o --allow-extra-chr --keep-allele-order --recode 12 transpose --output-missing-genotype 0 --out GP_new_o > /dev/null
singularity exec -B ${DATA_ROOT} $IMG /opt/emmax-20120210/emmax-kin-intel64 -v -x -d 10 -o GP_new_o.aBN.kinf GP_new_o > kin_new_x.log 2>&1
rm -f GP_new_o.tped
echo DONE > kin2.done
