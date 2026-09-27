IMG=${ANALYSIS_DIR}/05_GWAS/GWAS_tur/software/software/Reseq_genek.sif
cd ${CLUSTER_WORK}/snp_repair_r308/struct
singularity exec -B ${DATA_ROOT} $IMG plink --bfile all.LDfilter --allow-extra-chr --make-rel square --out grm > /dev/null
awk '{s+=$NR} END{print s}' grm.rel > grm_trace.txt
# same on the published LD-pruned set for calibration
singularity exec -B ${DATA_ROOT} $IMG plink --bfile ${ANALYSIS_DIR}/05_GWAS/00_analysis/02_genomeDB/05_snp/LD_pruned_bed --allow-extra-chr --make-rel square --out grm_old > /dev/null
awk '{s+=$NR} END{print s}' grm_old.rel > grm_old_trace.txt
