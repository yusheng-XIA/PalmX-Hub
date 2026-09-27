cd ${CLUSTER_WORK}/enh_gwas/cl
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4
PY=python
$PY s00_prep.py > ../logs/s00.log 2>&1 && $PY t01_validate.py > ../logs/t01.log 2>&1; echo DONE >> ../logs/t01.log
