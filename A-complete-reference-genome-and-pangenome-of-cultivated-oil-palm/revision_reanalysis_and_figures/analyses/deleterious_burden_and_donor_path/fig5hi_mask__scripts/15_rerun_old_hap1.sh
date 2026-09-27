#!/bin/bash
# clean re-run of nrly_old_hap1 (the first run's outputs were clobbered by a duplicate queue launch; moved to calls/_clobbered_nrly_old_hap1)
W=${CLUSTER_WORK}/fig5hi_mask
S=${ANALYSIS_DIR}/22_answer_reviews/00_ms/03_V3/01_figure1/00_final_input_resources/05_syri_plotsr_chain12/results/00_standardized_genomes
PY=python3
ulimit -v 199000000
mv $W/calls/nrly_old_hap1 $W/calls/_clobbered_nrly_old_hap1_$(date +%s)
bash $W/scripts/10_call_dsnp_asm.sh nrly_old_hap1 $S/nrly_hap1.chr16.fa > $W/logs/nrly_old_hap1.rerun.out 2> $W/logs/nrly_old_hap1.rerun.err
P=$W/calls/nrly_old_hap1/nrly_old_hap1_vs_Africa_hap2
$PY $W/scripts/20_coverage.py nrly_old_hap1 $P.sorted.paf $W/coverage/nrly_old_hap1.tsv >> $W/logs/post.log 2>&1
$PY $W/scripts/21_dsnp_windows.py nrly_old_hap1 $P.snp.txt $W/dsnp/nrly_old_hap1 > $W/logs/dsnp_nrly_old_hap1.log 2>&1
$PY $W/scripts/24_panel_sites.py > $W/logs/panel_sites.log 2>&1
$PY $W/scripts/25_repro_old.py > $W/logs/repro_old.log 2>&1
echo REPRODONE2 >> $W/logs/repro_old.log
