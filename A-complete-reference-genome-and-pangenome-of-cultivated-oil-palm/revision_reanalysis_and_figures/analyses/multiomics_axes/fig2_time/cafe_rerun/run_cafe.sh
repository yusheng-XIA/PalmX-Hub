#!/bin/bash
# CAFE5 re-runs for Fig. 2b (${COMPUTE_HOST}; <=15 threads total). Inputs copied read-only from user layer2_palm.
set -uo pipefail
W=${CLUSTER_WORK}/fig2_time
SRC=${ANALYSIS_DIR}/20_results/Figure2/07_new_figure/02_comparative_orthofinder/phylo_divtime_cafe/03_cafe5/layer2_palm
cd $W
export LD_LIBRARY_PATH=${CLUSTER_HOME}/miniconda3/lib
CAFE=$W/cafe_pkg/bin/cafe5
mkdir -p inputs runs
cp -n $SRC/gene_families.tsv inputs/gene_families.tsv
cp -n $SRC/cafe_tree.nwk inputs/cafe_tree_orig.nwk
md5sum $SRC/gene_families.tsv inputs/gene_families.tsv $SRC/cafe_tree.nwk inputs/cafe_tree_orig.nwk > inputs/md5.txt
$CAFE --help 2>&1 | head -3 > runs/cafe_version.txt
for T in orig oldMCMC newMCMC; do
  ( cd runs; /usr/bin/time -v $CAFE --infile ../inputs/gene_families.tsv --tree ../inputs/cafe_tree_$T.nwk --cores 5 -k 3 \
      --output_prefix gamma_$T > gamma_$T.log 2> gamma_$T.time ) &
done
wait
for T in orig oldMCMC newMCMC; do
  ( cd runs; /usr/bin/time -v $CAFE --infile ../inputs/gene_families.tsv --tree ../inputs/cafe_tree_$T.nwk --cores 5 \
      --output_prefix base_$T > base_$T.log 2> base_$T.time ) &
done
wait
echo ALLDONE > runs/DONE
