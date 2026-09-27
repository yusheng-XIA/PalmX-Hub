# 把sv_all.bed染色体名 chr01 → chr01B
sed 's/^chr\([0-9]*\)\t/chr\1B\t/' sv_all.bed > sv_all_hap2.bed

# 验证
echo "=== 修正后染色体名 ==="
cut -f1 sv_all_hap2.bed | sort -u

# 同时修正sv_inversions.bed
sed 's/^chr\([0-9]*\)\t/chr\1B\t/' sv_inversions.bed > sv_inversions_hap2.bed
