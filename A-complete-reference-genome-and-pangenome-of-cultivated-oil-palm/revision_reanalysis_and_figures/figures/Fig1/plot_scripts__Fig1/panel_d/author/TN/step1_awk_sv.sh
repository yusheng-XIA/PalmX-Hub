# 完全替代parse_syri_sv.py

grep -v "^#" SV.syri.out | \
  awk '
  BEGIN {
    split("INS|DEL|DUP|INVDP|TRANS|INVTR|CPL|CPG|INV|TDM", t, "|")
    for (i in t) valid[t[i]] = 1
  }
  {
    if (!($11 in valid)) next
    if ($2 == "." || $3 == ".") next

    # INS：长度取query序列长度（第5列），坐标用ref插入点
    if ($11 == "INS") {
      len = length($5)
      if (len < 50) next
      OFS="\t"
      print $1, $2-1, $2, $10, len, $11
    }
    # DEL及其他：长度取ref区间长度
    else {
      len = ($3 > $2) ? $3 - $2 : $2 - $3
      if (len < 50) next
      OFS="\t"
      print $1, $2-1, $3, $10, len, $11
    }
  }' | sort -k1,1 -k2,2n > sv_all.bed

# 单独提取倒位（图D用）
grep -v "^#" SV.syri.out | \
  awk '($11=="INV" || $11=="INVTR" || $11=="INVDP") && $2!="." && $3!="." {
    len = ($3 > $2) ? $3 - $2 : $2 - $3
    if (len < 50) next
    OFS="\t"
    print $1, $2-1, $3, $10, len, $11
  }' | sort -k1,1 -k2,2n > sv_inversions.bed

# 检查
echo "=== SV类型统计 ==="
cut -f6 sv_all.bed | sort | uniq -c | sort -rn
echo "Total: $(wc -l < sv_all.bed)"

echo ""
echo "=== 倒位统计 ==="
cut -f6 sv_inversions.bed | sort | uniq -c | sort -rn

echo ""
echo "=== 各类型长度分布（中位数bp）==="
cut -f5,6 sv_all.bed | sort -k2,2 | \
  awk '{a[$2][NR]=$1; c[$2]++} 
       END {for (t in c) {
         n=c[t]; mid=int(n/2)
         print t, a[t][mid], "bp (n="n")"
       }}' | sort -k1,1
