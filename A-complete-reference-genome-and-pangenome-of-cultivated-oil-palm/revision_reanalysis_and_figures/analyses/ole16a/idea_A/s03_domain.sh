#!/bin/bash
# Gate 1: domain architecture + phylogeny of all oleosin / caleosin / LDAP(REF) proteins in the FL+TN search DB
set -e
W=${CLUSTER_WORK}/idea_A; cd $W; mkdir -p out/dom tmp
R=${ANALYSIS_DIR}/22_answer_reviews/00_ms/03_V3/03_figure3/07_proteomics_multiomics_20260724/01_reference_database
E=${CLUSTER_HOME}/miniconda3/envs
HS=$E/interproscan/bin/hmmsearch; [ -x $HS ] || HS=$E/dram/bin/hmmsearch
for pf in PF01277 PF05042 PF05755; do gunzip -c $pf.hmm.gz > tmp/$pf.hmm; done
cat tmp/PF01277.hmm tmp/PF05042.hmm tmp/PF05755.hmm > tmp/ld3.hmm
$HS --cpu 8 --domtblout out/dom/unified_ld3.domtbl -E 1e-3 tmp/ld3.hmm $R/FL_TN_unified_exact_sequence_nr.fasta > /dev/null
$HS --cpu 8 --domtblout out/dom/ref_ld3.domtbl -E 1e-3 tmp/ld3.hmm ref_oleosins.fa > /dev/null
PY=$E/rnaseq_analysis/bin/python
$PY s03_domain.py
$E/phylo/bin/mafft --auto --quiet out/dom/oleosin_tree_input.fa > out/dom/oleosin_aln.fa
$E/phylo/bin/iqtree -s out/dom/oleosin_aln.fa -m MFP -B 1000 -T 8 --prefix out/dom/oleosin_iq -redo > out/dom/iqtree.log 2>&1 || $E/phylo/bin/iqtree2 -s out/dom/oleosin_aln.fa -m MFP -B 1000 -T 8 --prefix out/dom/oleosin_iq -redo > out/dom/iqtree.log 2>&1
echo DONE > out/dom/done.txt
