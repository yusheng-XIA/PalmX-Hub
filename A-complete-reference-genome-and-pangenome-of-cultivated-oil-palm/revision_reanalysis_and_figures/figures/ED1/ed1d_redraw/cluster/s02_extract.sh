#!/bin/bash
# ED1d redraw step 2: region BAMs (reads overlapping the display windows) + per-read / depth tables (pysam)
# usage: bash s02_extract.sh hifi|ont
# v2 (junctions re-derived in s06_junction.sh: chr12A 75,515|75,516; chr07B 129,402,880|129,402,881)
set -euo pipefail
PL=$1
W=${CLUSTER_WORK}/ed1d_redraw
ST=samtools
PY=python3
P=${ANALYSIS_DIR}/01_seedless_results/12_evaluate
mkdir -p $W/region $W/tables; cd $W
# locus  sample  chrom  win_start  win_end  (1-based, inclusive)
while read L S C A B; do
  if [ $PL = hifi ]; then BAI=${CLUSTER_WORK}/enh_A/bai/$S.hifi_remap.bai
  else BAI=$W/bai/$S.ont_remap.bai; fi
  $ST view -@ 4 -b -X $P/$S/02_remapping/${PL}_remap.sorted.bam $BAI $C:$A-$B -o region/${L}.${PL}.bam
  $ST index region/${L}.${PL}.bam
  # depth over the whole terminal extension (MAPQ>=20 and all) for the summary
done <<'LOCI'
chr12A_left 02_american chr12A 1 130000
chr07B_right 01_africa chr07B 129352881 129452880
LOCI
# extension-wide depth (whole added sequence)
if [ $PL = hifi ]; then BA=${CLUSTER_WORK}/enh_A/bai; else BA=$W/bai; fi
for spec in "chr12A_ext 02_american chr12A 1 75515" "chr07B_ext 01_africa chr07B 129402881 130504653" "chr12A_anchor 02_american chr12A 75516 1075515" "chr07B_anchor 01_africa chr07B 128402881 129402880"; do
  set -- $spec
  for q in 0 20; do
    $ST depth -a -Q $q -r $3:$4-$5 -X $P/$2/02_remapping/${PL}_remap.sorted.bam $BA/$2.${PL}_remap.bai \
      | awk -v L=$1 -v q=$q -v pl=$PL '{s+=$3; n++; if($3==0) z++} END{printf "%s\t%s\tMAPQ>=%s\t%d\t%.2f\t%d\n", L, pl, q, n, s/n, z}'
  done
done > tables/extension_depth.${PL}.tsv
$PY $W/s03_tables.py $PL
echo DONE $PL $(date)
