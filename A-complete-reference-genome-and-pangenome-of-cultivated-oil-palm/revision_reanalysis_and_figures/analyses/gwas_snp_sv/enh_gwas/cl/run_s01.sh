cd ${CLUSTER_WORK}/enh_gwas/cl
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
ulimit -v 40000000
python s01_lambda_orig.py 8 > ../logs/s01.log 2>&1; echo DONE >> ../logs/s01.log
