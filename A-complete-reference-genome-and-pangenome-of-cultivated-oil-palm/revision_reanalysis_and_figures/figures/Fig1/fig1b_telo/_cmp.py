import fitz, numpy as np, sys
def px(p):
    pm=p.get_pixmap(dpi=600); return np.frombuffer(pm.samples,np.uint8).reshape(pm.h,pm.w,pm.n).astype(int)
a=fitz.open(sys.argv[1]); b=fitz.open(sys.argv[2])
print('stream154 identical:', a.xref_stream(154)==b.xref_stream(154))
wa=[(round(w[0],2),round(w[3],2),w[4]) for w in a[0].get_text('words')]; wb=[(round(w[0],2),round(w[3],2),w[4]) for w in b[0].get_text('words')]
print('words identical (text+pos 0.01pt):', wa==wb, len(wa), len(wb))
if wa!=wb: print([ (x,y) for x,y in zip(wa,wb) if x!=y][:20])
A=px(a[0]); B=px(b[0]); D=np.abs(A-B).max(2)
print('pixels differing @600dpi:', int((D>0).sum()), ' >32:', int((D>32).sum()), ' >96:', int((D>96).sum()), 'max', D.max())
if (D>0).any():
    ys,xs=np.nonzero(D>0); print('bbox pt', xs.min()*72/600, ys.min()*72/600, xs.max()*72/600, ys.max()*72/600)
