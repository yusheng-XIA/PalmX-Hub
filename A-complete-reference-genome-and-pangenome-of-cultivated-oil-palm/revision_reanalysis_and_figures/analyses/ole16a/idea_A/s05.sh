#!/bin/bash
cd ${CLUSTER_WORK}/idea_A; ulimit -v 100000000; export OMP_NUM_THREADS=4
python s05_rna.py > s05.log 2>&1; echo EXIT $? >> s05.log
