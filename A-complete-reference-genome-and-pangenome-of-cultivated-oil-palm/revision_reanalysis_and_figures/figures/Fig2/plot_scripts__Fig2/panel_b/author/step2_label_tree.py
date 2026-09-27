"""
模仿教程: 将Gamma_clade_results.txt的扩缩数字写入树文件
所有节点都标注 +increase/-decrease
"""
import re, os

BASE = "${ANALYSIS_DIR}/20_results/Figure2/07_new_figure/02_comparative_orthofinder/phylo_divtime_cafe"
CAFE = f"{BASE}/03_cafe5"
OUT  = f"{BASE}/04_figures"

for layer_name in ["layer1_global", "layer2_palm", "layer3_oilcrop"]:
    print(f"\n--- {layer_name} ---")
    
    # 读取 clade_results
    gamma_data = {}
    clade_file = os.path.join(CAFE, layer_name, "gamma_results", "Gamma_clade_results.txt")
    with open(clade_file) as f:
        for line in f:
            if line.startswith('#'): continue
            parts = line.strip().split('\t')
            if len(parts) >= 3:
                taxon_id = parts[0].strip()
                increase = int(parts[1])
                decrease = int(parts[2])
                gamma_data[taxon_id] = (increase, decrease)
    
    print(f"  节点数: {len(gamma_data)}")
    
    # 读取 id_tree
    tree_file = os.path.join(OUT, f"{layer_name}_id_tree.nwk")
    with open(tree_file) as f:
        tree_content = f.read()
    
    # 替换所有节点: Name<ID> → Name<ID>+inc/-dec
    # 先处理有名字的节点 (tips): Dura<1> → Dura<1>+375/-576
    # 再处理内部节点: <4> → <4>+451/-649
    for taxon_id, (inc, dec) in gamma_data.items():
        pattern = re.escape(taxon_id)
        replacement = f"{taxon_id}+{inc}/-{dec}"
        tree_content = re.sub(pattern, replacement, tree_content)
    
    # 保存
    outfile = os.path.join(OUT, f"{layer_name}_labeled_tree.nwk")
    with open(outfile, 'w') as f:
        f.write(tree_content)
    print(f"  ✓ {layer_name}_labeled_tree.nwk")
    
    # 也生成显著性标注版本 (只标显著家族的扩缩)
    # 从Gamma_family_results.txt获取显著家族
    sig_file = os.path.join(CAFE, layer_name, "gamma_results", "Gamma_family_results.txt")
    change_file = os.path.join(CAFE, layer_name, "gamma_results", "Gamma_change.tab")
    
    # 读sig families
    sig_fams = set()
    with open(sig_file) as f:
        next(f)
        for line in f:
            p = line.strip().split('\t')
            if len(p) > 2 and p[2] == 'y':
                sig_fams.add(p[0])
    
    # 读change.tab获取每个节点在显著家族中的扩缩
    with open(change_file) as f:
        header = f.readline().strip().split('\t')
    
    col_nid = {}
    col_name = {}
    for i, col in enumerate(header):
        if i == 0: continue
        m = re.search(r'(.*)<(\d+)>', col)
        if m:
            name = m.group(1).strip()
            nid = int(m.group(2))
            col_nid[i] = nid
            col_name[nid] = name
    
    # 统计显著家族中每个节点的扩缩
    sig_stats = {}  # nid -> {expand: n, contract: n}
    for nid in col_nid.values():
        sig_stats[nid] = {'expand': 0, 'contract': 0}
    
    with open(change_file) as f:
        next(f)
        for line in f:
            parts = line.strip().split('\t')
            fam = parts[0]
            if fam not in sig_fams:
                continue
            for i, nid in col_nid.items():
                if i < len(parts):
                    val = int(parts[i].replace('+', ''))
                    if val > 0:
                        sig_stats[nid]['expand'] += 1
                    elif val < 0:
                        sig_stats[nid]['contract'] += 1
    
    # 生成数据表 (给R用)
    tsv_file = os.path.join(OUT, f"{layer_name}_node_data.tsv")
    with open(tsv_file, 'w') as f:
        f.write("node_id\tname\tis_tip\tall_expand\tall_contract\tall_net\tsig_expand\tsig_contract\tsig_net\n")
        for taxon_id, (inc, dec) in gamma_data.items():
            m = re.search(r'<(\d+)>', taxon_id)
            if m:
                nid = int(m.group(1))
                name = taxon_id.split('<')[0].strip()
                is_tip = name != ''
                se = sig_stats.get(nid, {}).get('expand', 0)
                sc = sig_stats.get(nid, {}).get('contract', 0)
                f.write(f"{nid}\t{name}\t{is_tip}\t{inc}\t{dec}\t{inc-dec}\t{se}\t{sc}\t{se-sc}\n")
    print(f"  ✓ {layer_name}_node_data.tsv")

print("\nStep 2 完成!")
