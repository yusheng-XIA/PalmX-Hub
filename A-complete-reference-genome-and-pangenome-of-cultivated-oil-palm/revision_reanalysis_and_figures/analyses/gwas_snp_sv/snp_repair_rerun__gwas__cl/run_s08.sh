cd ${CLUSTER_WORK}/enh_gwas/cl
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8
ulimit -v 40000000
python s08_sv_int.py > ../logs/s08.log 2>&1; echo DONE >> ../logs/s08.log
