"""K = 3 / K = 8 groups from the 308-based ADMIXTURE Q (argmax ancestry, as for K = 4). Columns are matched to the
current reordered Q of the same K (fig4new/structure, colours of Fig. 4a) by maximal correlation.  usage: kgroups.py K"""
import sys, numpy as np, pandas as pd
from pathlib import Path
from scipy.optimize import linear_sum_assignment
K = int(sys.argv[1])
W = Path("${CLUSTER_WORK}/snp_repair_r308"); OLD = Path("${CLUSTER_WORK}/snp_repair_rerun")
ids = [l.split()[1] for l in open(W / f"struct/admixK{K}/all.fam")]
qn = pd.read_csv(W / f"struct/admixK{K}/all.{K}.Q", sep=r"\s+", header=None); qn.index = ids
oids = [l.split()[1] for l in open(OLD / "fig4new/structure/all.fam")]
qo = pd.read_csv(OLD / f"fig4new/structure/all.{K}.Q", sep=r"\s+", header=None); qo.index = oids; qo = qo.loc[ids]
C = np.array([[np.corrcoef(qo[i], qn[j])[0, 1] for j in range(K)] for i in range(K)])
r, c = linear_sum_assignment(-C); qn = qn.iloc[:, list(c)]; qn.columns = range(K)
qn.to_csv(W / f"struct/admixK{K}/all.{K}.Q.reordered", sep=" ", header=False, index=False, float_format="%.6f")
pop = qn.to_numpy().argmax(1) + 1
G = W / f"groups_K{K}"; G.mkdir(exist_ok=True)
for p in range(1, K + 1): open(G / f"K{K}_Pop{p}.txt", "w").write("".join(f"{s}\n" for s, q in zip(ids, pop) if q == p))
pd.DataFrame({"id": ids, "pop": pop, "maxQ": qn.to_numpy().max(1)}).to_csv(G / "assignments.tsv", sep="\t", index=False)
po = qo.to_numpy().argmax(1) + 1
print(f"K{K} matched r", np.round(C[r, c], 4), "sizes", [int((pop == p).sum()) for p in range(1, K + 1)],
      "changed vs current Q", int((po != pop).sum()), "max|dQ|", round(float(np.abs(qn.to_numpy() - qo.to_numpy()).max()), 4))
