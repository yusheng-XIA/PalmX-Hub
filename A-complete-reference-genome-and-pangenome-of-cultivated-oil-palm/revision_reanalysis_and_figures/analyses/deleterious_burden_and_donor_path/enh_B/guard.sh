#!/bin/bash
# kill my own minimap2 (PID $1) if node MemAvailable drops below 8 GB
P=$1; L=${CLUSTER_WORK}/enh_B/aln/guard.log
echo "guard start $(date) pid $P" >> $L
while kill -0 $P 2>/dev/null; do
  a=$(awk '/MemAvailable/{print int($2/1048576)}' /proc/meminfo)
  if [ "$a" -lt 8 ] && grep -q asm20 /proc/$P/cmdline; then kill $P; echo "killed_lowmem avail=${a}G $(date)" >> $L; fi
  sleep 10
done
echo "guard end $(date)" >> $L
