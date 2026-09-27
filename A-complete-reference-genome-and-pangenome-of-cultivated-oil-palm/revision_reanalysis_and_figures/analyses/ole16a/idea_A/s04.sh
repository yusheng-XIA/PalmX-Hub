#!/bin/bash
W=${CLUSTER_WORK}/idea_A; cd $W
PY=python
[ -d pylib/pyarrow ] || $PY -m pip install --quiet --no-index --no-deps --target pylib pyarrow-14.0.2-cp38-cp38-manylinux_2_17_x86_64.manylinux2014_x86_64.whl
export PYTHONPATH=$W/pylib OMP_NUM_THREADS=4
ulimit -v 150000000
$PY s04_pep.py > s04.log 2>&1; echo EXIT $? >> s04.log
