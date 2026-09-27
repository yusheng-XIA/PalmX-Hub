#!/bin/bash
cd ${CLUSTER_WORK}/idea_A; ulimit -v 100000000
python s07_protfam.py > s07.log 2>&1; echo EXIT $? >> s07.log
