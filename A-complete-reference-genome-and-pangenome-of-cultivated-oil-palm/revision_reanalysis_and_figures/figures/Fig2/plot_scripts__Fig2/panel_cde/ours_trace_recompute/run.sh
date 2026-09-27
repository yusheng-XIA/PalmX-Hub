#!/bin/bash
cd ${CLUSTER_WORK}/trace/Fig2; mkdir -p out
PY=python
for p in met rna prot; do nohup $PY recompute_fig2cde.py $p > out/log_$p.txt 2>&1 & done
wait
echo ALLDONE > out/done.txt
