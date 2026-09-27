#!/bin/bash
O=${CLUSTER_WORK}/headline_pop/g6
mkdir -p $O/geno $O/tmp; export TMPDIR=$O/tmp
cd $O
python3 ${CLUSTER_WORK}/headline_pop/cl/g6_a_cds.py > $O/a.log 2>&1
ls $O/bed | sed 's/.bed//' | nice -n 10 xargs -P 8 -I{} bash ${CLUSTER_WORK}/headline_pop/cl/g6_b_geno.sh {}
echo ALLDONE >> $O/geno/done.txt
