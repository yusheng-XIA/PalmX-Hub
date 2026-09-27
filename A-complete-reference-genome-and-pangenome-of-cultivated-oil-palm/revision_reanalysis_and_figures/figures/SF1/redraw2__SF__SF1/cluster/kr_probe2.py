import pickle, numpy as np, sys
H="${ANALYSIS_DIR}/22_answer_reviews/00_ms/03_V3/01_figure1/02_hic_circos_10haps_20260814/RUN-OP41-HICCIRCOS10-20260814-001/results/hic"
M,ns=pickle.load(open(f"{H}/BK_hap1/contact_matrix.pkl","rb"))
print(vars(ns))
A=M.astype(float)
print("sym",np.allclose(A,A.T),"diag sum",np.trace(A),"total",A.sum(),"zero rows",(A.sum(1)==0).sum())
def kr(A,tol=1e-6,delta=0.1,Delta=3):
    n=A.shape[0]; e=np.ones(n); x0=e.copy(); g=0.9; etamax=0.1; eta=etamax; stop_tol=tol*0.5
    x=x0; rt=tol**2; v=x*(A@x); rk=1-v; rho_km1=rk@rk; rout=rho_km1; rold=rout; it=0
    while rout>rt:
        it+=1; k=0; y=e.copy(); innertol=max(eta**2*rout,rt)
        while rho_km1>innertol:
            k+=1
            if k==1:
                Z=rk/v; p=Z.copy(); rho_km1=rk@Z
            else:
                beta=rho_km1/rho_km2; p=Z+beta*p
            w=x*(A@(x*p))+v*p; alpha=rho_km1/(p@w); ap=alpha*p; ynew=y+ap
            if ynew.min()<=delta:
                if delta==0: break
                ind=ap<0; gamma=np.min((delta-y[ind])/ap[ind]); y=y+gamma*ap; break
            if ynew.max()>=Delta:
                ind=ynew>Delta; gamma=np.min((Delta-y[ind])/ap[ind]); y=y+gamma*ap; break
            y=ynew; rk=rk-alpha*w; rho_km2=rho_km1; Z=rk/v; rho_km1=rk@Z
        x=x*y; v=x*(A@x); rk=1-v; rho_km1=rk@rk; rout=rho_km1
        rat=rout/rold; rold=rout; res_norm=np.sqrt(rout); eta_o=eta; eta=g*rat
        if g*eta_o**2>0.1: eta=max(eta,g*eta_o**2)
        eta=max(min(eta,etamax),stop_tol/res_norm)
    return x
keep=A.sum(1)>0
B=A[np.ix_(keep,keep)]
x=kr(B)
N=np.zeros_like(A); N[np.ix_(keep,keep)]=x[:,None]*B*x[None,:]
print("rowsum",N[keep].sum(1).min(),N[keep].sum(1).max())
target=0.006742300978060431/4
import itertools
agp=[l.split() for l in open(f"{H}/BK_hap1/BK_hap1.identity.agp") if l.strip()]
lens=[int(r[2]) for r in agp]; nb=[ -(-L//500000) for L in lens]; print("nbins sum",sum(nb))
edges=np.cumsum([0]+nb)
intra=np.zeros_like(N,bool)
for a,b in zip(edges[:-1],edges[1:]): intra[a:b,a:b]=True
C={}
C["diag median"]=np.median(np.diag(N)); C["diag mean"]=np.diag(N).mean()
C["offdiag1 median"]=np.median(np.diag(N,1))
C["intra median"]=np.median(N[intra]); C["intra nonzero median"]=np.median(N[intra&(N>0)])
C["intra mean"]=N[intra].mean()
for q in (90,95,99): C[f"p{q}"]=np.percentile(N,q)
C["rowmax median"]=np.median(N.max(1))
# raw-scale variants
C["raw median x?"]=np.median(A)
for k,v in C.items(): print(k, v, v/target)
