import fitz, numpy as np, sys
from render_fig3k import X0,X1,Y0,Y1
S='${WORK_DIR}'
O=S+'/deliver/Main_Figures_revised/Figure3.pdf'; N=sys.argv[1]
R=fitz.Rect(X0,Y0,X1,Y1)
po,pn=fitz.open(O)[0],fitz.open(N)[0]
print('page size',po.rect,pn.rect)
def outside(ws): return [w[4] for w in ws if not fitz.Rect(w[:4]).intersects(R)]
wo=po.get_text('words'); wn=pn.get_text('words')
oo,on=outside(wo),outside(wn)
print('words outside rect: orig',len(oo),'new',len(on),'identical',oo==on)
# also compare words with positions
so=sorted((round(w[0],2),round(w[1],2),w[4]) for w in wo if not fitz.Rect(w[:4]).intersects(R))
sn=sorted((round(w[0],2),round(w[1],2),w[4]) for w in wn if not fitz.Rect(w[:4]).intersects(R))
print('outside words+positions identical',so==sn)
print('inside new words:',[w[4] for w in wn if fitz.Rect(w[:4]).intersects(R)])
z=600/72
a=po.get_pixmap(matrix=fitz.Matrix(z,z),alpha=False); b=pn.get_pixmap(matrix=fitz.Matrix(z,z),alpha=False)
A=np.frombuffer(a.samples,np.uint8).reshape(a.h,a.w,3).astype(int); B=np.frombuffer(b.samples,np.uint8).reshape(b.h,b.w,3).astype(int)
d=np.abs(A-B).max(2); m=np.ones_like(d,bool)
m[int(Y0*z):int(np.ceil(Y1*z)),int(X0*z):int(np.ceil(X1*z))]=False
print('600dpi px',A.shape,'max diff outside rect',d[m].max(),'n px >0 outside',(d[m]>0).sum())
b.save(N.replace('.pdf','_600dpi.png'))
fitz.open(N)[0].get_pixmap(matrix=fitz.Matrix(z,z),clip=fitz.Rect(268,366,424.8,512),alpha=False).save(N.replace('.pdf','_3k_zoom600.png'))
fonts=set()
for f in pn.get_fonts(full=True): fonts.add(f[3])
print('fonts',fonts)
