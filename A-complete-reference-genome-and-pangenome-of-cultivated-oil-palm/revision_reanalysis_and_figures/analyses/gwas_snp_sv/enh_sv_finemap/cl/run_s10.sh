cd ${CLUSTER_WORK}/enh_sv_finemap/cl
mkdir -p ../logs ../out ../tmp
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8 TMPDIR=${CLUSTER_WORK}/enh_sv_finemap/tmp
ulimit -v 60000000
python s10_finemap.py > ../logs/s10.log 2>&1; echo DONE >> ../logs/s10.log
