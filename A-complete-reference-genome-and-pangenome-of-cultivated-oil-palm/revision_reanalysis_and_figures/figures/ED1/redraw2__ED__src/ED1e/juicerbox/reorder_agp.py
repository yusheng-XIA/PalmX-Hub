#!/usr/bin/env python3

# 读取排序信息
order_map = {}
with open('group_chr_order_complete.txt', 'r') as f:  # 改为complete文件
    for line in f:
        group, chr_name, order = line.strip().split('\t')
        order_map[group] = (chr_name, int(order))

# 读取原始AGP
agp_lines = {}
with open('scaffolds.raw.agp', 'r') as f:
    for line in f:
        parts = line.strip().split('\t')
        group = parts[0]
        agp_lines[group] = line.strip()

# 按新顺序输出，替换group名为染色体名
sorted_groups = sorted(order_map.items(), key=lambda x: x[1][1])

with open('scaffolds_reordered.agp', 'w') as f:
    for group, (chr_name, order) in sorted_groups:
        if group in agp_lines:
            # 替换group名为染色体名
            line_parts = agp_lines[group].split('\t')
            line_parts[0] = chr_name
            f.write('\t'.join(line_parts) + '\n')

print(f"Successfully created scaffolds_reordered.agp with {len(sorted_groups)} chromosomes")
