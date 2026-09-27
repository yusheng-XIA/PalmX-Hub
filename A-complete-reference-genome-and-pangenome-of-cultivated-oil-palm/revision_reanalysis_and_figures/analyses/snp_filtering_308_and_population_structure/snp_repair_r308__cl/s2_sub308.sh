#!/bin/bash
# 308-accession subset of the repaired joint-called diversity call set, then PLINK sets mirroring the published 308 pipeline:
#  structure: --geno 0.1 --maf 0.01 --biallelic-only strict (run_evolution_pipeline.sh)
#  GWAS:      --geno 0.1 --maf 0.05 (shared_genotype/all_308.log)
set -euo pipefail
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
c=$1
E=${CLUSTER_HOME}/miniconda3/envs/genomics_a2/bin
W=${CLUSTER_WORK}/snp_repair_r308
F=${JOINT_SNP_DIR}/06_filter
IMG=${ANALYSIS_DIR}/05_GWAS/GWAS_tur/software/software/Reseq_genek.sif
S308=${CLUSTER_WORK}/fig4c_final/meta/s308.txt
IN=$W/vcf/$c.div308.bcf   # 308-based filter chain (r308/filter308.sh)
O=$W/s308; mkdir -p $O
$E/bcftools view --threads 4 -S $S308 -Ou $IN | $E/bcftools +fill-tags -Ou -- -t AN,AC,F_MISSING,MAF \
  | $E/bcftools view --threads 4 -i "INFO/AC>0 && INFO/AC<INFO/AN" -Ob -o $O/$c.div308.bcf
$E/bcftools index --threads 4 $O/$c.div308.bcf
echo -e "$c\t$($E/bcftools index -n $O/$c.div308.bcf)" > $O/$c.div308.count
cd $O
singularity exec -B ${DATA_ROOT} $IMG plink --bcf $c.div308.bcf --geno 0.1 --maf 0.01 --biallelic-only strict --allow-extra-chr --set-missing-var-ids @:# --keep-allele-order --make-bed --out $c.struct > /dev/null
singularity exec -B ${DATA_ROOT} $IMG plink --bcf $c.div308.bcf --geno 0.1 --maf 0.05 --allow-extra-chr --set-missing-var-ids @:# --keep-allele-order --make-bed --out $c.gwas > /dev/null
echo DONE > $c.s2.done
