#!/bin/bash
W=${CLUSTER_WORK}/ole16_cis; cd $W; mkdir -p sv
V=${GWAS_DIR}/02_analysis/10_snp_sv_gwas/01_pangenie/02_sv_combined/03_combined_matrix/pangenie_joint.biallelic.vcf.gz
BC=$(ls ${CLUSTER_HOME}/miniconda3/envs/*/bin/bcftools | head -1); echo $BC
cut -f1 samples308.tsv > sv/s308.txt
$BC view -r chr11B:112184000-112200000 -S sv/s308.txt --force-samples $V 2>/dev/null | $BC query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT[\t%GT]\n' | awk -v OFS='\t' '{print $1,$2,$3,length($4),length($5),$0}' | cut -f1-5,11- > sv/pangenie_region.tsv
$BC query -l -S sv/s308.txt --force-samples $V 2>/dev/null > /dev/null
$BC view -h -r chr11B:1-1 -S sv/s308.txt --force-samples $V 2>/dev/null | $BC query -l > sv/sample_order.txt
wc -l sv/pangenie_region.tsv
C=${ANALYSIS_DIR}/22_answer_reviews/00_ms/03_V3/05_figure/08_SV_hap39_Figure5abcd_20260808/results/final_catalog/oilpalm_hap38.final_highconfidence_sv.tsv
head -1 $C > sv/linear_catalog_region.tsv
awk -F'\t' 'NR>1' $C | awk -F'\t' '{for(i=1;i<=NF;i++) if($i ~ /^chr11/){c=i;break}} {print}' | grep -P "chr11B?\t11218[4-9][0-9]{3}|chr11B?\t11219[0-9]{4}" >> sv/linear_catalog_region.tsv
wc -l sv/linear_catalog_region.tsv
