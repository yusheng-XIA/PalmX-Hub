"""Inputs for Fig. 4a/4b candidates from the repaired call set (new ADMIXTURE Q, new PCA), in the author's layout."""
import numpy as np, pandas as pd, shutil
from pathlib import Path
from scipy.optimize import linear_sum_assignment
W = Path("${CLUSTER_WORK}/snp_repair_r308")
OLD = Path("${ANALYSIS_DIR}/05_GWAS/00_analysis/02_genomeDB/08_structure")
OR = Path("${ANALYSIS_DIR}/22_answer_reviews/00_ms/03_V3/04_figure4/K3_K4_pi_fst_pca_groups")
N = W / "fig4new"; (N / "structure").mkdir(parents=True, exist_ok=True); (N / "metadata").mkdir(exist_ok=True); (N / "summary").mkdir(exist_ok=True); (N / "figures").mkdir(exist_ok=True)
new_ids = [l.split()[1] for l in open(W / "struct/admix/all.fam")]
old_ids = [l.split()[1] for l in open(OLD / "all.fam")]
shutil.copy(W / "struct/admix/all.fam", N / "structure/all.fam")
shutil.copy(OLD / "integrated_rooted_tree_structure_sample_order.txt", N / "structure/")
perm = {}
for k in range(2, 9):
    if not (W / f"struct/admix/all.{k}.Q").exists(): continue
    qn = pd.read_csv(W / f"struct/admix/all.{k}.Q", sep=r"\s+", header=None); qn.index = new_ids
    try: qo = pd.read_csv(OLD / f"all.{k}.Q", sep=r"\s+", header=None); qo.index = old_ids
    except FileNotFoundError: continue
    qo = qo.loc[new_ids]
    C = np.array([[np.corrcoef(qo[i], qn[j])[0, 1] for j in range(k)] for i in range(k)])
    r, c = linear_sum_assignment(-C); perm[k] = list(c)
    qn.iloc[:, list(c)].to_csv(N / f"structure/all.{k}.Q", sep=" ", header=False, index=False, float_format="%.6f")
    print(k, "old->new cols", list(c), "matched r", np.round(C[r, c], 3))
q4 = pd.read_csv(N / "structure/all.4.Q", sep=r"\s+", header=None); q4.index = new_ids; q4.columns = ["Q1", "Q2", "Q3", "Q4"]
pc = pd.read_csv(W / "struct/PCA_10.eigenvec", sep=r"\s+", header=None).set_index(1)
po = pd.read_csv("${ANALYSIS_DIR}/05_GWAS/00_analysis/02_genomeDB/05_snp/PCA_10.eigenvec", sep=r"\s+", header=None).set_index(1)
for c in range(2, 12):   # PC signs are arbitrary: orient each new PC to correlate positively with the published PC of the same rank
    if np.corrcoef(pc.loc[new_ids, c], po.loc[new_ids, c])[0, 1] < 0: pc[c] = -pc[c]
pc.reset_index()[[0, 1] + list(range(2, 12))].to_csv(N / "structure/PCA_10_oriented.eigenvec", sep=" ", header=False, index=False)
arch = pd.read_csv(OR / "metadata/K4_dominant_assignments.tsv", sep="\t").set_index("sample")["old_6group"]
m = pd.DataFrame(index=new_ids)
m["sample"] = new_ids; m["K"] = "K4"; mx = q4.max(1); dom = q4.idxmax(1).str.replace("Q", "K4_Pop")
m["dominant_group"] = dom; m["threshold_0.7_group"] = np.where(mx >= 0.7, dom, "Admixed"); m["max_Q"] = mx
m["second_Q"] = q4.apply(lambda r: sorted(r)[-2], axis=1); m["old_6group"] = arch.reindex(new_ids).values
for i, p in enumerate(["PC1", "PC2", "PC3"]): m[p] = pc.loc[new_ids, i + 2].values
m = m.join(q4)
m.to_csv(N / "metadata/K4_dominant_assignments.tsv", sep="\t", index=False)
print("K4 sizes", m.dominant_group.value_counts().sort_index().to_dict(), "maxQ<0.7", int((mx < 0.7).sum()))
import json
v = json.load(open(W / "cl/fig4c_values_k4new.json"))
pops = ["K4_Pop1", "K4_Pop2", "K4_Pop3", "K4_Pop4"]
pd.DataFrame({"group": pops, "PI_weighted": [v["pi"][p[3:]] / 1000 for p in pops], "n_used_weighted": 0, "n_windows": 0}).to_csv(N / "summary/K4_PI_weighted.tsv", sep="\t", index=False)
M = pd.DataFrame(0.0, index=pops, columns=pops)
for k, x in v["fst"].items():
    a, b = ["K4_" + t for t in k.split("_")]; M.loc[a, b] = M.loc[b, a] = x
M.to_csv(N / "summary/K4_FST_weighted_matrix.tsv", sep="\t")
