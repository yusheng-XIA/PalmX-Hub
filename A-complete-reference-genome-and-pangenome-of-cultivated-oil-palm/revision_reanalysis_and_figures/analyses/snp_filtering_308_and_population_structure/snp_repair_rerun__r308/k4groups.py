"""Map the new K = 4 ADMIXTURE columns (308-based structure set) to the Pop1-4 labels of the current groups by maximal
Q-matrix correlation; write groups_k4new/ (previous copy kept as groups_k4new_prev/) and a comparison summary."""
import shutil, numpy as np, pandas as pd
from pathlib import Path
from scipy.optimize import linear_sum_assignment
W = Path("${CLUSTER_WORK}/snp_repair_r308"); OLD = Path("${CLUSTER_WORK}/snp_repair_rerun")
ids = [l.split()[1] for l in open(W / "struct/admix/all.fam")]
qn = pd.read_csv(W / "struct/admix/all.4.Q", sep=r"\s+", header=None); qn.index = ids
oids = [l.split()[1] for l in open(OLD / "fig4new/structure/all.fam")]
qo = pd.read_csv(OLD / "fig4new/structure/all.4.Q", sep=r"\s+", header=None); qo.index = oids; qo = qo.loc[ids]
C = np.array([[np.corrcoef(qo[i], qn[j])[0, 1] for j in range(4)] for i in range(4)])
r, c = linear_sum_assignment(-C); qn = qn.iloc[:, list(c)]; qn.columns = [0, 1, 2, 3]
print("matched r", np.round(C[r, c], 4))
prev = W / "groups_k4new_prev"
if not prev.exists(): shutil.copytree(W / "groups_k4new", prev)
a0 = pd.read_csv(prev / "assignments.tsv", sep="\t").set_index("id")
pop = qn.to_numpy().argmax(1) + 1; mx = qn.to_numpy().max(1)
a = pd.DataFrame({"id": ids, "pop": pop, "maxQ": mx}).set_index("id")
a["arch"] = a0.loc[ids, "arch"]; a["k4old"] = a0.loc[ids, "k4old"]
a.reset_index().to_csv(W / "groups_k4new/assignments.tsv", sep="\t", index=False)
for p in range(1, 5): open(W / f"groups_k4new/K4_Pop{p}.txt", "w").write("".join(f"{s}\n" for s in a.index[a["pop"] == p]))
qn.to_csv(W / "struct/admix/all.4.Q.reordered", sep=" ", header=False, index=False, float_format="%.6f")
ch = a0.loc[ids, "pop"] != a["pop"]
print("sizes", a["pop"].value_counts().sort_index().to_dict(), "maxQ<0.7", int((mx < 0.7).sum()))
print("changed vs current groups", int(ch.sum()), [(s, int(a0.loc[s, "pop"]), int(a.loc[s, "pop"]), round(float(a.loc[s, "maxQ"]), 3)) for s in a.index[ch]])
print("max |dQ| vs current", float(np.abs(qn.to_numpy() - qo.to_numpy()).max()), "mean |dQ|", float(np.abs(qn.to_numpy() - qo.to_numpy()).mean()))
