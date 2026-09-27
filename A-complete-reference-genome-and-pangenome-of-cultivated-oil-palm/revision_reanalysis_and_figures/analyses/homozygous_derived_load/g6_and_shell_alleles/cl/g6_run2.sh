#!/bin/bash
O=${CLUSTER_WORK}/headline_pop/g6
mkdir -p $O/snpeff; rm -f $O/snpeff/done.txt
ls $O/bed | sed 's/.bed//' | xargs -P 8 -I{} bash ${CLUSTER_WORK}/headline_pop/cl/g6_f_snpeff.sh {}
ulimit -v 80000000
cd $O && python ../cl/g6_e_sites.py > e.log 2>&1
