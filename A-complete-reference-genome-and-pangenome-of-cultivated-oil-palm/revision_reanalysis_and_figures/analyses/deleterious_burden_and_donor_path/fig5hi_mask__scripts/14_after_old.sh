#!/bin/bash
W=${CLUSTER_WORK}/fig5hi_mask
until grep -q ALLPOST $W/logs/post.log 2>/dev/null; do sleep 300; done
ulimit -v 60000000
python3 $W/scripts/25_repro_old.py > $W/logs/repro_old.log 2>&1
echo REPRODONE >> $W/logs/repro_old.log
