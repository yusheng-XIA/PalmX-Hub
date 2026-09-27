#!/bin/bash
# run on ${COMPUTE_HOST}: KR-normalise the FL Pore-C matrix (ED1e redraw2)
W=${CLUSTER_HOME}/OilPalm_final_0923/ed_redraw2
cd $W
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8
nice -n 10 ${DATA_DIR2}/anaconda3/envs/haphic/bin/python $W/ed1e_kr.py $W/ED1e_KR_500kb.npz > $W/ed1e_kr.log 2>&1
ls -la $W >> $W/ed1e_kr.log
