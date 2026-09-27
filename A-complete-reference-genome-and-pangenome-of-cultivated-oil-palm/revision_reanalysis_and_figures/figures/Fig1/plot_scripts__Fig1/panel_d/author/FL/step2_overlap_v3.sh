#!/bin/bash
# step2_overlap_v3.sh

SV=sv_all_hap2.bed
TE=Africa_hap2.fa.mod.EDTA.TEanno.gff3
GENE=Africa_hap2.EVM.gff3
GENOME=Africa_hap2.genome

set -euo pipefail
echo "[$(date '+%H:%M:%S')] 开始 step2_overlap 分析"

# ════════════════════════════════════════════
# 1. 准备 TE BED
# ════════════════════════════════════════════
echo "[$(date '+%H:%M:%S')] 1/6 提取TE区间..."

grep -v "^#" $TE | \
  awk '$3 != "long_terminal_repeat" && \
       $3 != "target_site_duplication" && \
       $3 != "repeat_region"' | \
  awk '{
    split($9, a, ";")
    class = "unknown"
    for (i in a) {
      if (a[i] ~ /^classification=/) {
        split(a[i], b, "=")
        class = b[2]
      }
    }
    OFS="\t"
    print $1, $4-1, $5, $3, class, $7
  }' | sort -k1,1 -k2,2n > te.bed

echo "  TE条目数: $(wc -l < te.bed)"

# ════════════════════════════════════════════
# 2. 准备基因各区域 BED
# ════════════════════════════════════════════
echo "[$(date '+%H:%M:%S')] 2/6 提取基因区域..."

# 基因体
grep -v "^#" $GENE | awk '$3=="gene" && $3!=""' | \
  awk '{OFS="\t"; print $1,$4-1,$5,$9,".",$7}' | \
  sort -k1,1 -k2,2n > gene_body.bed

# 外显子
grep -v "^#" $GENE | awk '$3=="exon" && $3!=""' | \
  awk '{OFS="\t"; print $1,$4-1,$5,$9,".",$7}' | \
  sort -k1,1 -k2,2n > exon.bed

# 内含子 = 基因体 - 外显子
bedtools subtract -a gene_body.bed -b exon.bed | \
  sort -k1,1 -k2,2n > intron.bed

# 上游2kb（不考虑链方向）
bedtools flank -i gene_body.bed -g $GENOME -l 2000 -r 0 | \
  sort -k1,1 -k2,2n > upstream2k.bed

# 下游2kb（不考虑链方向）
bedtools flank -i gene_body.bed -g $GENOME -l 0 -r 2000 | \
  sort -k1,1 -k2,2n > downstream2k.bed

# 合并genic区域（基因体 + 上下游2kb）
cat gene_body.bed upstream2k.bed downstream2k.bed | \
  sort -k1,1 -k2,2n | \
  bedtools merge > genic_region.bed

echo "  基因数:    $(wc -l < gene_body.bed)"
echo "  外显子数:  $(wc -l < exon.bed)"
echo "  内含子数:  $(wc -l < intron.bed)"

bedtools genomecov -i genic_region.bed -g $GENOME | \
  awk '$1=="genome" && $2==1 {print "  genic覆盖率:", $5}'
bedtools genomecov -i te.bed -g $GENOME | \
  awk '$1=="genome" && $2==1 {print "  TE覆盖率:   ", $5}'

# ════════════════════════════════════════════
# 3. 外圈分类：genic vs intergenic（互斥，覆盖全部SV）
#    判断标准：SV自身50%以上落在genic region即为genic
# ════════════════════════════════════════════
echo "[$(date '+%H:%M:%S')] 3/6 外圈：genic vs intergenic..."

bedtools intersect -a $SV -b genic_region.bed -u -f 0.5 > sv_genic.bed
bedtools intersect -a $SV -b genic_region.bed -v -f 0.5 > sv_intergenic.bed

echo "  genic:      $(wc -l < sv_genic.bed)"
echo "  intergenic: $(wc -l < sv_intergenic.bed)"

# ════════════════════════════════════════════
# 4. 内圈分类：TE overlap vs no TE overlap（互斥，覆盖全部SV）
#    判断标准：SV自身50%以上与TE重叠即为TE overlap
# ════════════════════════════════════════════
echo "[$(date '+%H:%M:%S')] 4/6 内圈：TE overlap vs no TE overlap..."

bedtools intersect -a $SV -b te.bed -u -f 0.5 > sv_with_te.bed
bedtools intersect -a $SV -b te.bed -v -f 0.5 > sv_no_te.bed

echo "  与TE重叠≥50%: $(wc -l < sv_with_te.bed)"
echo "  与TE重叠<50%: $(wc -l < sv_no_te.bed)"

# ════════════════════════════════════════════
# 5. 右侧条形图：仅针对genic SV的子区域细分
#    优先级（互斥）: exon > intron > upstream2k > downstream2k
# ════════════════════════════════════════════
echo "[$(date '+%H:%M:%S')] 5/6 genic SV子区域细分..."

# 先把genic SV与各子区域intersect
bedtools intersect -a sv_genic.bed -b exon.bed       -u -f 0.5 > sv_exon_raw.bed
bedtools intersect -a sv_genic.bed -b intron.bed     -u -f 0.5 > sv_intron_raw.bed
bedtools intersect -a sv_genic.bed -b upstream2k.bed -u -f 0.5 > sv_up2k_raw.bed
bedtools intersect -a sv_genic.bed -b downstream2k.bed -u -f 0.5 > sv_down2k_raw.bed

# 按优先级互斥分配
# exon优先级最高，直接取
cp sv_exon_raw.bed sv_exon.bed

# intron：去掉已归入exon的
bedtools intersect -a sv_intron_raw.bed \
  -b sv_exon.bed -v > sv_intron.bed

# upstream：去掉已归入exon和intron的
cat sv_exon.bed sv_intron.bed | sort -k1,1 -k2,2n | uniq > _assigned.bed
bedtools intersect -a sv_up2k_raw.bed \
  -b _assigned.bed -v > sv_upstream2k.bed

# downstream：去掉已归入exon、intron、upstream的
cat sv_exon.bed sv_intron.bed sv_upstream2k.bed | \
  sort -k1,1 -k2,2n | uniq > _assigned.bed
bedtools intersect -a sv_down2k_raw.bed \
  -b _assigned.bed -v > sv_downstream2k.bed

rm -f sv_exon_raw.bed sv_intron_raw.bed sv_up2k_raw.bed \
      sv_down2k_raw.bed _assigned.bed

echo "  exon:        $(wc -l < sv_exon.bed)"
echo "  intron:      $(wc -l < sv_intron.bed)"
echo "  upstream2k:  $(wc -l < sv_upstream2k.bed)"
echo "  downstream2k:$(wc -l < sv_downstream2k.bed)"

# ════════════════════════════════════════════
# 6. 汇总写入CSV
# ════════════════════════════════════════════
echo "[$(date '+%H:%M:%S')] 6/6 写入汇总文件..."

total=$(wc -l < $SV)
cat > sv_summary_for_plot.csv << EOF
total,${total}
genic,$(wc -l < sv_genic.bed)
intergenic,$(wc -l < sv_intergenic.bed)
te_overlap,$(wc -l < sv_with_te.bed)
no_te,$(wc -l < sv_no_te.bed)
exon,$(wc -l < sv_exon.bed)
intron,$(wc -l < sv_intron.bed)
up2k,$(wc -l < sv_upstream2k.bed)
down2k,$(wc -l < sv_downstream2k.bed)
EOF

echo ""
echo "════════════════════════════════"
echo "汇总结果："
cat sv_summary_for_plot.csv
echo "════════════════════════════════"
echo "[$(date '+%H:%M:%S')] step2 完成"
