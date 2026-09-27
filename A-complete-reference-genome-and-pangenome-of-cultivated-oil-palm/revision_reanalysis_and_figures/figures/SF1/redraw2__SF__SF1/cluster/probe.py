import pickle, numpy as np, sys
H="${ANALYSIS_DIR}/22_answer_reviews/00_ms/03_V3/01_figure1/02_hic_circos_10haps_20260814/RUN-OP41-HICCIRCOS10-20260814-001/results/hic"
obj=pickle.load(open(f"{H}/BK_hap1/contact_matrix.pkl","rb"))
def desc(o,d=0):
    print(" "*d, type(o), getattr(o,"shape",""), getattr(o,"dtype",""))
    if isinstance(o,(tuple,list)):
        for x in o[:6]: desc(x,d+2)
    if isinstance(o,dict):
        for k in list(o)[:6]: print(" "*d,k); desc(o[k],d+2)
desc(obj)
