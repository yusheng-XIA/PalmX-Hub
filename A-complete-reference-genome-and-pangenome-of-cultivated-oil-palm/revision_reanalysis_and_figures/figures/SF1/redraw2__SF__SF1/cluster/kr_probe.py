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
target=0.006742300978060431
for nm,v in [("median all",np.median(N)),("median nonzero",np.median(N[N>0])),("median keep",np.median(N[np.ix_(keep,keep)])),
             ("nanmedian upper",np.median(N[np.triu_indices_from(N)])), ("mean",N.mean())]:
    print(nm, 4*v, 4*v/target)
