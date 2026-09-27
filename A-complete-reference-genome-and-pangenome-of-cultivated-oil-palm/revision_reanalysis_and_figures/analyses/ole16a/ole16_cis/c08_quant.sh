#!/bin/bash
W=${CLUSTER_WORK}/ole16_cis; cd $W; mkdir -p rq; export TMPDIR=$W/tmp
MM=minimap2
for f in hits/*.hits.tsv; do r=$(basename $f .hits.tsv); [ -s rq/$r.sam ] && continue
  awk -F'\t' '{print ">"$1"\n"$2}' $f > rq/$r.fa
  $MM -ax sr -N 5 --secondary=yes -t 8 targets.fa rq/$r.fa 2>/dev/null | grep -v "^@" | cut -f1-6,12- > rq/$r.sam
  rm rq/$r.fa
done
echo ok
