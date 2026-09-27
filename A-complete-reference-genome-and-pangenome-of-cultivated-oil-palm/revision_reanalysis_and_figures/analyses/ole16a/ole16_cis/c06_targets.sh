#!/bin/bash
D=${ANALYSIS_DIR}/08_hifi_chromosome/11_final_corrected_39genomes_annotations_20260812/attempt_20260812_01
W=${CLUSTER_WORK}/ole16_cis; cd $W
ST=samtools
> targets.fa; > targets.bed
while IFS=: read g id; do
  [ "$g" = "Africa_hap2" ] && G=African_hap2 || G=American_hap1
  line=$(grep -P "\tgene\t" $D/02_annotations/$G.gene.gff3 | grep -F "ID=$id;" | head -1)
  [ -z "$line" ] && { echo "miss $g $id"; continue; }
  c=$(echo "$line"|cut -f1); s=$(echo "$line"|cut -f4); e=$(echo "$line"|cut -f5); st=$(echo "$line"|cut -f7)
  a=$((s-500)); b=$((e+500))
  if [ "$st" = "-" ]; then $ST faidx -i $D/01_genomes/$G.fa $c:$a-$b | sed "s/^>.*/>${G}__${id}/" >> targets.fa; else $ST faidx $D/01_genomes/$G.fa $c:$a-$b | sed "s/^>.*/>${G}__${id}/" >> targets.fa; fi
  echo -e "${G}__${id}\t$c\t$s\t$e\t$st" >> targets.bed
done < targets_ids.txt
grep -c ">" targets.fa
