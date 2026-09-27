#!/bin/bash
# usage: 12_launch.sh NODE   (run on ${LOGIN_HOST}; only launches)
W=${CLUSTER_WORK}/fig5hi_mask
case $1 in
 ${COMPUTE_HOST}) Q="nrly_final_hap1:${ANALYSIS_DIR}/08_hifi_chromosome/11_final_corrected_39genomes_annotations_20260812/attempt_20260812_01/01_genomes/nrly_hap1.fa dura_hap2:${ANALYSIS_DIR}/22_answer_reviews/00_ms/03_V3/01_figure1/00_final_input_resources/05_syri_plotsr_chain12/results/00_standardized_genomes/dura_hap2.chr16.fa nrly_old_hap1:${ANALYSIS_DIR}/22_answer_reviews/00_ms/03_V3/01_figure1/00_final_input_resources/05_syri_plotsr_chain12/results/00_standardized_genomes/nrly_hap1.chr16.fa";;
 ${COMPUTE_HOST}) Q="pisifera_hap1:${ANALYSIS_DIR}/22_answer_reviews/00_ms/03_V3/01_figure1/00_final_input_resources/05_syri_plotsr_chain12/results/00_standardized_genomes/pisifera_hap1.chr16.fa pisifera_hap2:${ANALYSIS_DIR}/22_answer_reviews/00_ms/03_V3/01_figure1/00_final_input_resources/05_syri_plotsr_chain12/results/00_standardized_genomes/pisifera_hap2.chr16.fa";;
esac
ssh -n $1 "cd $W && setsid nohup bash scripts/11_queue.sh '$Q' > logs/q_$1.log 2>&1 < /dev/null &"
