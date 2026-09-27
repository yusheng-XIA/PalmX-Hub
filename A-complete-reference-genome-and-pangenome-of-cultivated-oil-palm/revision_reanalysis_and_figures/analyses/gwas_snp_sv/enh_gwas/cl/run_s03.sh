cd ${CLUSTER_WORK}/enh_gwas/cl
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4
ulimit -v 60000000
PY=python
$PY s03_summarize.py > ../logs/s03.log 2>&1; echo DONE >> ../logs/s03.log
$PY s04_sv_complement.py > ../logs/s04.log 2>&1; echo DONE >> ../logs/s04.log
