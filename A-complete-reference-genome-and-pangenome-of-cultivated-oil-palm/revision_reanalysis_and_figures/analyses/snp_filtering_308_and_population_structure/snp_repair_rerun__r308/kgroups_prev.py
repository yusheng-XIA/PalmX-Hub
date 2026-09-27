"""K = 3 / 8 groups (argmax) from the existing ADMIXTURE Q of the previous structure set (0.16% fewer pruned SNPs than the
308-based set), used while the 308-based ADMIXTURE runs.  usage: kgroups_prev.py K"""
import sys, pandas as pd
from pathlib import Path
K = int(sys.argv[1]); W = Path("${CLUSTER_WORK}/snp_repair_r308"); OLD = Path("${CLUSTER_WORK}/snp_repair_rerun")
ids = [l.split()[1] for l in open(OLD / "fig4new/structure/all.fam")]
q = pd.read_csv(OLD / f"fig4new/structure/all.{K}.Q", sep=r"\s+", header=None); pop = q.to_numpy().argmax(1) + 1
G = W / f"groups_K{K}_prevQ"; G.mkdir(exist_ok=True)
for p in range(1, K + 1): open(G / f"K{K}_Pop{p}.txt", "w").write("".join(f"{s}\n" for s, x in zip(ids, pop) if x == p))
pd.DataFrame({"id": ids, "pop": pop, "maxQ": q.to_numpy().max(1)}).to_csv(G / "assignments.tsv", sep="\t", index=False)
print(K, [int((pop == p).sum()) for p in range(1, K + 1)])
