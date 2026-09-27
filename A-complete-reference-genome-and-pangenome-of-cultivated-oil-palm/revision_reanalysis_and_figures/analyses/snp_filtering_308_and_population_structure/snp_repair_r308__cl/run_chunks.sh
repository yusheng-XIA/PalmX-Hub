#!/bin/bash
# usage: run_chunks.sh LIST PARALLEL
W=${CLUSTER_WORK}/snp_repair_r308
cat $1 | xargs -P $2 -L 1 bash -c "c=\${1%%:*}; o=\$2; [ -s \$o.done ] && exit 0; bash $W/cl/filter_chr.sh \$0 \${1%%:*} \$2 \$1 2 > \$2.log 2>&1 || echo FAIL \$1 >> $W/logs/chunk_fail.txt"
