#!/bin/bash
# Divergence times with MCMCTree (PAML v4.10.9), approximate-likelihood method
# Inputs: supergene.phy, calibrated_tree.txt (seven uniform calibrations, unit = 100 Ma), jones.dat
# (empirical amino-acid model file from PAML/dat, required by codeml in step 1)
set -euo pipefail

# Step 1 (usedata = 3): gradient and Hessian for the approximate likelihood
mcmctree mcmctree_step1.ctl
if [ -s out.BV ]; then cp out.BV in.BV; else cp rst2 in.BV; fi

# Step 2 (usedata = 2): MCMC; burn-in 50,000, 20,000 samples every 50 iterations
mcmctree mcmctree.ctl

# Posterior means and 95% HPD intervals are in out.txt; the dated tree is FigTree.tre.
# Convergence: run twice and check ESS > 200 for all parameters (e.g. with Tracer on mcmc.txt).
