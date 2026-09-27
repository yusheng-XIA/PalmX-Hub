import sys, numpy as np, pandas as pd
sys.path.insert(0, "${CLUSTER_WORK}/snp_repair_rerun/gwas/cl")
from config import RUN, BED, BIM, FAM, KIN_SNP
from emmaxpy import GLS, read_bed_rows, read_matrix, read_pheno
O="${CLUSTER_WORK}/snp_repair_rerun/gwas/shell/out"; import os; LEAD=int(os.environ.get("LEAD","3153030"))
ids=[l.split()[1] for l in open(FAM)]; N=len(ids)
bim=pd.read_csv(BIM,sep=r"\s+",header=None,names=["c","id","cm","pos","a1","a2"])
want={"lead":LEAD,"MPOB_L29P":3259207,"AVROS_site":3259200}
rows={}
for k,p in want.items():
    i=int(np.nonzero((bim.c=="chr01B").values&(bim.pos==p).values)[0][0])
    g=read_bed_rows(str(BED),N,i,i+1)[0].astype(float); g[g<0]=np.nan
    rows[k]=(g,bim.iloc[i].a1,bim.iloc[i].a2)
    print(k,p,"A1",bim.iloc[i].a1,"A2",bim.iloc[i].a2,"mean A1 dosage",np.nanmean(g), "miss",np.isnan(g).sum())
K=read_matrix(str(KIN_SNP))
tasks=pd.read_csv(RUN/"manifests/sv_tasks.tsv",sep="\t"); cat=dict(zip(tasks.trait,tasks.category))
import glob
def pheno(tr):
    fs=glob.glob(str(RUN/f"phenotypes/*/{tr}.txt")); return read_pheno(fs[0],ids)
# mutant allele dosages: MPOB mutant = G (alt); AVROS mutant = ref A (FL-Hap2 codes N)
def mut_dose(key,mut):
    g,a1,a2=rows[key]; return g if a1==mut else 2-g
# verified vs VCF: A1 of 3259207 = ALT G (sh-MPOB, L29P); A1 of 3259200 = REF A (sh-AVROS, N31; ALT T restores K)
mp=mut_dose("MPOB_L29P","G"); av=mut_dose("AVROS_site","A"); ld=rows["lead"][0]
sh=mp+av
out=[]
tab=pd.DataFrame({"id":ids,"MPOB":mp,"AVROS":av,"sh_total":sh,"lead_A1":ld})
for tr in ["Nut_weight_g","Shell_thickness_mm","Shell_weight_g","Flesh_thickness_mm","Nut_length_mm","Fruit_weight_g","Mesocarp_to_fruit_ratio"]:
    fs=glob.glob(str(RUN/f"phenotypes/*/{tr}.txt"))
    if not fs: print("no pheno",tr); continue
    y=pheno(tr); tab[tr]=y
    k=np.isfinite(y)&np.isfinite(sh)&np.isfinite(ld)
    X=np.ones((k.sum(),1)); g=GLS(y[k],X,K[np.ix_(k,k)])
    G=np.vstack([mp[k],av[k],sh[k],ld[k],(sh[k]>0).astype(float)])
    Z=[X,(sh[k]>0).astype(float)[:,None]]+([(sh[k]>1).astype(float)[:,None]] if 0<(sh[k]>1).sum()<k.sum() else []); gc=GLS(y[k],np.hstack(Z),K[np.ix_(k,k)]); pl2=gc.scan(ld[k][None,:])[2][0]
    gl=GLS(y[k],np.c_[X,ld[k]],K[np.ix_(k,k)]); pc2=gl.scan((sh[k]>0).astype(float)[None,:])[2][0]
    b,se,p=g.scan(G)
    # conditional: lead | sh ; sh | lead
    g2=GLS(y[k],np.c_[X,sh[k]],K[np.ix_(k,k)]); pl=g2.scan(ld[k][None,:])[2][0]
    g3=GLS(y[k],np.c_[X,ld[k]],K[np.ix_(k,k)]); ps=g3.scan(sh[k][None,:])[2][0]
    r2=np.corrcoef(sh[k],ld[k])[0,1]**2
    med={v:float(np.median(y[k][sh[k]==v])) if (sh[k]==v).sum() else np.nan for v in [0,1,2]}
    nn={v:int((sh[k]==v).sum()) for v in [0,1,2]}
    out.append(dict(trait=tr,n=int(k.sum()),P_MPOB=p[0],P_AVROS=p[1],P_shTotal=p[2],P_lead=p[3],P_shCarrier=p[4],
                    r2_sh_lead=r2,P_lead_given_sh=pl,P_sh_given_lead=ps,P_lead_given_shClass=pl2,P_carrier_given_lead=pc2,n_sh0=nn[0],n_sh1=nn[1],n_sh2=nn[2],
                    med_sh0=med[0],med_sh1=med[1],med_sh2=med[2]))
r=pd.DataFrame(out); pd.set_option("display.width",250); print(r.to_string())
r.to_csv(O+"/g3_shell_known_alleles_v2.tsv",sep="\t",index=False); tab.to_csv(O+"/g3_shell_genotypes_v2.tsv",sep="\t",index=False)
print(pd.crosstab(tab.MPOB,tab.AVROS))
