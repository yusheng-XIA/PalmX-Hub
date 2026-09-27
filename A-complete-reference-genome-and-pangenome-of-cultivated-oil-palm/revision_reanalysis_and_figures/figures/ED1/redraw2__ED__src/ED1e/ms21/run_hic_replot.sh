#!/bin/bash
# ============================================================
# Boke Hi-C 重新过滤 + 绘图
# 去掉 --nm 3 限制, 输出到 03_figure
# ============================================================

source ~/.bashrc
conda activate haphic

SRC="${DATA_DIR2}/projects/1_oil_palm/7-manuscripts/1-SVs-validation/2-boken-sv"
OUT="${ANALYSIS_DIR}/21_MS/01_result/03_figure"
mkdir -p "$OUT"

cd "$SRC"

# ============================================================
# Step 1: 重新过滤 (去掉 --nm 3, 只保留 MAPQ >= 1)
# ============================================================
echo "[$(date)] Step 1: Re-filtering BAM (MAPQ>=1, no NM limit) ..."

~/tools/HapHiC/utils/filter_bam HiC.bam 1 --threads 14 \
  | samtools view - -b -@ 14 -o "$OUT/HiC.filtered.relaxed.bam"

echo "[$(date)] Filtering done."
echo "对比 reads 数:"
echo -n "  旧 (NM<3):  "; samtools view -c -@ 14 HiC.filtered.bam
echo -n "  新 (无NM):  "; samtools view -c -@ 14 "$OUT/HiC.filtered.relaxed.bam"

# ============================================================
# Step 2: 绘图 (haphic plot 输出到 03_figure)
# ============================================================
echo ""
echo "[$(date)] Step 2: Plotting ..."

cd "$OUT"

~/tools/HapHiC/haphic plot \
  "$SRC/chr_asm.agp" \
  "$OUT/HiC.filtered.relaxed.bam" \
  --bin_size 1000 --min_len 5

echo ""
echo "[$(date)] Done!"
echo "输出目录: $OUT"
ls -lh "$OUT"/contact_map* 2>/dev/null
ls -lh "$OUT"/contact_matrix* 2>/dev/null

# ============================================================
# Step 3: 清理大文件 (可选, 取消注释执行)
# ============================================================
# echo "清理中间 BAM ..."
# rm "$OUT/HiC.filtered.relaxed.bam"
