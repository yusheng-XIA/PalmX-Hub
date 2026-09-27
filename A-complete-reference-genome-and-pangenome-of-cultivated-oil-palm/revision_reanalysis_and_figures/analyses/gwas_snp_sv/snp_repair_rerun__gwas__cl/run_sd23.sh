cd ${CLUSTER_WORK}/snp_repair_rerun/gwas/cl
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4
PY=python
$PY s23_sd23.py > ../logs/s23.log 2>&1
$PY s22_far_finemap.py > ../logs/s22.log 2>&1
echo DONE >> ../logs/s23.log
