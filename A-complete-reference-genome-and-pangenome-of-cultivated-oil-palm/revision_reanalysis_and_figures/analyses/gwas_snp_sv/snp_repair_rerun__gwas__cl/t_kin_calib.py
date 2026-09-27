"""Effect of kinship implementation (container emmax-kin aBN -x vs published hBN) on the published genotypes:
Nut weight, chr01B SHELL window (+-500 kb of 3,153,030) and lambda over chr01B."""
import sys, numpy as np, pandas as pd
sys.path.insert(0, "${CLUSTER_WORK}/enh_gwas/cl")
from config import RUN, BED, FAM, KIN_SNP, W
from emmaxpy import GLS, read_bed_rows, read_matrix, read_pheno
from scipy.stats import chi2
ids=[l.split()[1] for l in open(FAM)]; N=len(ids)
y=read_pheno(RUN/"phenotypes/yield/Nut_weight_g.txt", ids); k=np.isfinite(y)
Kh=read_matrix(KIN_SNP); Ka=np.loadtxt("${CLUSTER_WORK}/snp_repair_rerun/gwas/kincal/old_aBN_x.kinf")
rng=pd.read_csv(W/"in/chr_ranges.tsv",sep="\t").set_index("chr").loc["chr01B"]; pos=np.load(W/"in/pos_chr01B.npy")
out={}
for nm,K in [("hBN",Kh),("aBN",Ka)]:
    g=GLS(y[k],np.ones((k.sum(),1)),K[np.ix_(k,k)]); ps=[]
    for s in range(0,int(rng.n),200000):
        e=min(int(rng.n),s+200000); G=read_bed_rows(BED,N,int(rng.start)+s,int(rng.start)+e); ps.append(g.scan(G[:,k])[2])
    p=np.concatenate(ps); out[nm]=p
    j=int(np.argmin(p)); print(nm,"lead",pos[j],p[j],"lambda_chr01B",np.median(chi2.isf(np.clip(p,1e-300,1),1))/chi2.ppf(.5,1), "n_sig", int((p<1.7735e-9).sum()))
d=np.abs(np.log10(out["hBN"])-np.log10(out["aBN"])); print("max |dlog10P|",d.max(),"median",np.median(d))
