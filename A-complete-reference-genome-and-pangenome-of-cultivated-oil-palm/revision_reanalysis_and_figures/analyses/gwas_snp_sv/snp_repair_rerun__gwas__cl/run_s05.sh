cd ${CLUSTER_WORK}/enh_gwas/cl
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
PY=python
ls ../in/pos_*.npy | sed 's/.*pos_//; s/.npy//' | xargs -P 8 -I{} bash -c "ulimit -v 8000000; $PY s05_mac.py {} >> ../logs/s05.log 2>&1"
echo DONE >> ../logs/s05.log
