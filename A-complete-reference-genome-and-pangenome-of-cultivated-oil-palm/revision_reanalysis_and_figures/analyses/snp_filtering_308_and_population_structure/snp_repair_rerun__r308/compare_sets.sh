#!/bin/bash
# sites added/removed by the 308-based filters, per set and chromosome (old = current applied results)
O=${CLUSTER_WORK}/snp_repair_rerun; N=${CLUSTER_WORK}/snp_repair_r308
E=${CLUSTER_HOME}/miniconda3/envs/genomics_a2/bin; T=$N/tmp/cmp; mkdir -p $T
out=$N/logs/compare_sets.tsv; echo -e "set\tchrom\told\tnew\tshared\tremoved\tadded" > $out
cmpf() { local s=$1 c=$2 a=$3 b=$4
  sh=$(comm -12 $a $b | wc -l); o=$(wc -l < $a); n=$(wc -l < $b)
  echo -e "$s\t$c\t$o\t$n\t$sh\t$((o-sh))\t$((n-sh))" >> $out; }
for i in $(seq -w 1 16); do c=chr${i}B
  for s in struct gwas; do
    awk '{print $4}' $O/s308/$c.$s.bim | sort > $T/o.$s.$c; awk '{print $4}' $N/s308/$c.$s.bim | sort > $T/n.$s.$c; cmpf $s $c $T/o.$s.$c $T/n.$s.$c
  done
  $E/bcftools query -f '%POS\n' $O/s308/$c.div308.bcf | sort > $T/o.div.$c; $E/bcftools query -f '%POS\n' $N/vcf/$c.div308.bcf | sort > $T/n.div.$c; cmpf div308 $c $T/o.div.$c $T/n.div.$c
  $E/bcftools query -f '%POS\n' $O/fig4c/vcf/$c.present.vcf.gz | sort > $T/o.pf.$c; $E/bcftools query -f '%POS\n' $N/fig4c/vcf/$c.present.vcf.gz | sort > $T/n.pf.$c; cmpf pi_fst $c $T/o.pf.$c $T/n.pf.$c
  zcat $O/g6/geno/$c.tsv.gz | awk 'NR>1{print $2}' | sort > $T/o.g6.$c; zcat $N/g6/geno/$c.tsv.gz | awk 'NR>1{print $2}' | sort > $T/n.g6.$c; cmpf g6_cds $c $T/o.g6.$c $T/n.g6.$c
done
python3 - $out <<'P'
import sys, csv, collections
t = collections.defaultdict(lambda: [0]*5)
for r in csv.DictReader(open(sys.argv[1]), delimiter="\t"):
    for i, k in enumerate(["old", "new", "shared", "removed", "added"]): t[r["set"]][i] += int(r[k])
for s, v in t.items(): print(s, "old", v[0], "new", v[1], "removed", v[3], f"({100*v[3]/v[0]:.3f}%)", "added", v[4], f"({100*v[4]/v[0]:.3f}%)")
P
rm -rf $T
