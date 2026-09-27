#!/usr/bin/env python3
"""Pangenome and core-genome accumulation over random assembly-addition orders.

Input: family x assembly presence table (tab-separated, first column family id, 0/1 or gene counts).
For each of --perm random orders (seed 42) the number of families present in >= 1 (pan) and in all (core)
of the first n assemblies is recorded; mean and sample standard deviation are reported for each n.
Families are also classified by occupancy: core (all), soft-core (all but one), shell (>= 2), cloud (1).
"""
import argparse

import numpy as np
import pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument("--families", required=True)
ap.add_argument("--out", required=True)
ap.add_argument("--perm", type=int, default=1000)
ap.add_argument("--seed", type=int, default=42)
a = ap.parse_args()

m = pd.read_csv(a.families, sep="\t", index_col=0)
P = (m.to_numpy() > 0)
n_asm = P.shape[1]
occ = P.sum(1)
cls = np.select([occ == n_asm, occ == n_asm - 1, occ >= 2], ["core", "soft-core", "shell"], "cloud")
print(pd.Series(cls).value_counts().to_string())

rng = np.random.default_rng(a.seed)
pan = np.zeros((a.perm, n_asm)); core = np.zeros((a.perm, n_asm))
for i in range(a.perm):
    order = rng.permutation(n_asm)
    seen_any = np.zeros(P.shape[0], bool); seen_all = np.ones(P.shape[0], bool)
    for j, k in enumerate(order):
        seen_any |= P[:, k]; seen_all &= P[:, k]
        pan[i, j] = seen_any.sum(); core[i, j] = seen_all.sum()
pd.DataFrame({"n_assemblies": np.arange(1, n_asm + 1),
              "pan_mean": pan.mean(0), "pan_sd": pan.std(0, ddof=1),
              "core_mean": core.mean(0), "core_sd": core.std(0, ddof=1)}).to_csv(a.out, sep="\t", index=False)
