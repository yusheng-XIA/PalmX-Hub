#!/bin/bash
source ~/.bashrc
conda activate haphic

HAPHIC="${DATA_DIR2}/tools/HapHiC"
SRC="${DATA_DIR2}/projects/1_oil_palm/7-manuscripts/1-SVs-validation/2-boken-sv"
OUT="${ANALYSIS_DIR}/21_MS/01_result/03_figure"
mkdir -p "$OUT"

cd "$SRC"

echo "[$(date)] Step 1: Re-filtering BAM (MAPQ>=1, no NM limit) ..."
$HAPHIC/utils/filter_bam HiC.bam 1 --threads 14 \
  | samtools view - -b -@ 14 -o "$OUT/HiC.filtered.relaxed.bam"

echo "[$(date)] Filtering done."
echo -n "  旧 (NM<3):  "; samtools view -c -@ 14 HiC.filtered.bam
echo -n "  新 (无NM):  "; samtools view -c -@ 14 "$OUT/HiC.filtered.relaxed.bam"

echo ""
echo "[$(date)] Step 2: Plotting ..."
cd "$OUT"
$HAPHIC/haphic plot "$SRC/chr_asm.agp" "$OUT/HiC.filtered.relaxed.bam" \
  --bin_size 1000 --min_len 5

echo ""
echo "[$(date)] Done!"
ls -lh "$OUT"/contact_map* 2>/dev/null
