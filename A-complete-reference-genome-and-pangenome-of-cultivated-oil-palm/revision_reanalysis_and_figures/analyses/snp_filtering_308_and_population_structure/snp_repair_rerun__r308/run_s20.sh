cd ${CLUSTER_WORK}/snp_repair_r308/gwas/cl
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4
PY=python
$PY s20_shell_finemap.py > ../logs/s20.log 2>&1
L=$(awk -F"\t" "NR>1 && \$2==\"Nut_weight_g\"{split(\$0,a,\"\t\")} END{}" ../shell/out/finemap_shell_new.tsv)
LEAD=$($PY -c "import pandas as pd; d=pd.read_csv(\"../shell/out/finemap_shell_new.tsv\",sep=\"\t\"); print(int(d[d.trait==\"Nut_weight_g\"].snp_lead_pos.iloc[0]))")
echo LEAD=$LEAD >> ../logs/s20.log
LEAD=$LEAD $PY s21_shell_alleles.py > ../logs/s21.log 2>&1
echo DONE >> ../logs/s21.log
