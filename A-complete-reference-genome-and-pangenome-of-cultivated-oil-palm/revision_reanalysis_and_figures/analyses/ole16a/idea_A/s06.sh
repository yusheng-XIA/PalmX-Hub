#!/bin/bash
cd ${CLUSTER_WORK}/idea_A; ulimit -v 120000000; export OMP_NUM_THREADS=4
python s06_sn.py > s06.log 2>&1; echo EXIT $? >> s06.log
