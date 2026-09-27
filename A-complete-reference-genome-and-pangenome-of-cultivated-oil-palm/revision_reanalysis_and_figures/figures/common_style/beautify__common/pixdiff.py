import sys, fitz, numpy as np
def r(p):
    x=fitz.open(p)[0].get_pixmap(dpi=100,alpha=False); return np.frombuffer(x.samples,np.uint8).reshape(x.height,x.width,3).astype(int)
a,b=r(sys.argv[1]),r(sys.argv[2])
if a.shape!=b.shape: print("SHAPE", a.shape,b.shape); sys.exit()
d=np.abs(a-b).sum(2); print(f"mean {d.mean():.4f} px>30: {(d>30).sum()}")
