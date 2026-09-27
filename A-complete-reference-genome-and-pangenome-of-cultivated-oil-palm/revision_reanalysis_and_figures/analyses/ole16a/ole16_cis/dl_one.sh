#!/bin/bash
# usage: dl_one.sh RUN URL NREADS
r=$1; u=$2; n=$3; cd "$(dirname "$0")/dl"
[ -s $r.done ] && exit 0
for try in 1 2 3; do
  curl -s -m 1800 "https://$u" | gunzip -c 2>/dev/null | head -n $((n*4)) | awk -v OFS='\t' 'NR%4==1{h=substr($1,2)} NR%4==2{print h,$0; c++} END{print c > "/dev/stderr"}' 2> $r.count | rg -F -f ../pats31.txt > $r.hits.tsv
  c=$(cat $r.count); if [ "${c:-0}" -ge $((n/2)) ]; then echo "$r $c" > $r.done; exit 0; fi
done
echo "$r FAIL ${c:-0}" > $r.fail
