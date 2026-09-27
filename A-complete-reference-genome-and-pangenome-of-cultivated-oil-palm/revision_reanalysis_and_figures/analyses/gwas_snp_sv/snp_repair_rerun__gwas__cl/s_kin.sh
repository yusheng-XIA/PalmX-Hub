IMG=${ANALYSIS_DIR}/05_GWAS/GWAS_tur/software/software/Reseq_genek.sif
cd ${CLUSTER_WORK}/snp_repair_rerun/gwas/geno
singularity exec -B ${DATA_ROOT} $IMG emmax-kin-intel64 -v -d 10 GP_new > kin.log 2>&1 && rm -f GP_new.tped && echo DONE > geno.done
