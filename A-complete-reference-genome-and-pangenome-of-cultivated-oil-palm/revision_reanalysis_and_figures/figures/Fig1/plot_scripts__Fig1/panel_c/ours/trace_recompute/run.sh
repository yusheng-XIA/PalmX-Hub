#!/bin/bash
cd ${CLUSTER_WORK}/trace/Fig1c
for s in gene te gap sd trf syn; do
  nohup python3 recompute.py $s > log_$s.txt 2>&1 &
done
wait
