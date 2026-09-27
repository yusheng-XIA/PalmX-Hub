import re, sys, os

BASE = "${ANALYSIS_DIR}/20_results/Figure2/07_new_figure/02_comparative_orthofinder/phylo_divtime_cafe"
CAFE = f"{BASE}/03_cafe5"
MCMC = f"{BASE}/02_mcmctree"
OUT  = f"{BASE}/04_figures"

for layer_name, layer_dir in [("layer1_global", "layer1_global"), 
                               ("layer2_palm", "layer2_palm"),
                               ("layer3_oilcrop", "layer3_oilcrop")]:
    print(f"\n--- {layer_name} ---")
    
    asr_file = os.path.join(CAFE, layer_dir, "gamma_results", "Gamma_asr.tre")
    with open(asr_file) as f:
        for line in f:
            if line.strip().startswith("TREE "):
                tree_line = line.strip()
                break
    
    # 提取 = 后面的树
    nwk = tree_line.split("=", 1)[1].strip().rstrip(";")
    
    # 清理: 去掉 *标记 和 _count (保留<ID>和:brlen)
    # 原始: Dura<1>*_10:12  →  Dura<1>:12
    # 内部: )<4>_14:4  →  )<4>:4
    cleaned = re.sub(r'\*', '', nwk)
    cleaned = re.sub(r'(<\d+>)_\d+', r'\1', cleaned)
    cleaned = cleaned.strip() + ";"
    
    # 保存 id_tree (带<ID>的newick)
    with open(os.path.join(OUT, f"{layer_name}_id_tree.nwk"), 'w') as f:
        f.write(cleaned + "\n")
    
    print(f"  ✓ {layer_name}_id_tree.nwk")
    print(f"    前200字符: {cleaned[:200]}")

# 同时处理MCMCTree时间树 (用于ggtree带时间轴)
print(f"\n--- MCMCTree clean tree ---")
with open(os.path.join(MCMC, "mcmctree_clean.nwk")) as f:
    nwk = f.read().strip()
# 确保冒号后无空格
nwk = re.sub(r':\s+', ':', nwk)
with open(os.path.join(OUT, "mcmctree_clean.nwk"), 'w') as f:
    f.write(nwk + "\n")
print(f"  ✓ mcmctree_clean.nwk")

print("\nStep 1 完成!")
