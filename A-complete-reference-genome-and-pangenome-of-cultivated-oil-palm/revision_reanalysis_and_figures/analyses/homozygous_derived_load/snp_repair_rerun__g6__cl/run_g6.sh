cd ${CLUSTER_WORK}/snp_repair_rerun/g6
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2
ulimit -v 80000000
PY=python
$PY cl/g6_e_sites.py > e.log 2>&1 && $PY cl/g6_checks.py > g6_checks.out 2>&1
echo DONE $? > g6.done
