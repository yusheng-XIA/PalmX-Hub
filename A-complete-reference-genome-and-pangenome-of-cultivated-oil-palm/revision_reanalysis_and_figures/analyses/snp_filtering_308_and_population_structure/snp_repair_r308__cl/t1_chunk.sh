set -e
W=${CLUSTER_WORK}/snp_repair_r308; E=${CLUSTER_HOME}/miniconda3/envs/genomics_a2/bin
RAW=${JOINT_SNP_DIR}/04_chr_vcf/chr12B.raw.vcf.gz
bash $W/cl/filter_chr.sh $RAW chr12B $W/test/c1 chr12B:1-1000000 4 2>/dev/null &
bash $W/cl/filter_chr.sh $RAW chr12B $W/test/c2 chr12B:1000001-2000000 4 2>/dev/null &
wait
for n in diversity population; do $E/bcftools concat -a -Oz -o $W/test/cc.$n.vcf.gz $W/test/c1.$n.vcf.gz $W/test/c2.$n.vcf.gz; $E/bcftools query -f "%POS[\t%GT]\n" $W/test/cc.$n.vcf.gz | md5sum; $E/bcftools query -f "%POS[\t%GT]\n" $W/test/chr12B_2Mb.$n.vcf.gz | md5sum; done
