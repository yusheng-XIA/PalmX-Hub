IMG=${ANALYSIS_DIR}/05_GWAS/GWAS_tur/software/software/Reseq_genek.sif
B=${ANALYSIS_DIR}/05_GWAS/00_analysis/06_GWAS/shared_genotype
W=${CLUSTER_WORK}/snp_repair_rerun/gwas/kincal; mkdir -p $W; cd $W
singularity exec -B ${DATA_ROOT} $IMG /opt/emmax-20120210/emmax-kin-intel64 -v -d 10 -o old_aBN.kinf $B/GP_maf_allChr > kin_old.log 2>&1
python3 - <<'P'
import numpy as np
a=np.loadtxt('old_aBN.kinf'); b=np.loadtxt('${ANALYSIS_DIR}/05_GWAS/00_analysis/06_GWAS/shared_genotype/GP_maf_allChr.hBN.kinf')
print('maxabsdiff',np.abs(a-b).max(),'corr',np.corrcoef(a.ravel(),b.ravel())[0,1],'diag',a.diagonal()[:3],b.diagonal()[:3])
P
echo DONE > kin_old.done
