#!/bin/bash
W=${CLUSTER_WORK}/idea_A; cd $W
cat out/ann_hits_africahap2.tsv
E=${CLUSTER_HOME}/miniconda3/envs
for e in phylo interproscan dram bio biotools; do echo "== $e"; ls $E/$e/bin 2>/dev/null | grep -E "^(hmmscan|hmmsearch|mafft|iqtree2?|FastTree|blastp|diamond|muscle|trimal|interproscan.sh|clustalo)$"; done
find $E/dram $E/interproscan -maxdepth 4 -iname "Pfam-A.hmm*" 2>/dev/null | head
find ${CLUSTER_HOME} -maxdepth 4 -iname "Pfam-A.hmm" 2>/dev/null | head
