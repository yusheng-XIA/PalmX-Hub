import fitz,re,csv,numpy as np,openpyxl
d=fitz.open('../src/Figure4_submitted_copy.pdf'); s=d.xref_stream(78).decode('latin1')
def parse_rel(seg):
    m=re.search(r'1 0 0 1 ([\d.]+) ([\d.]+) cm',seg); x0,y0=float(m.group(1)),float(m.group(2))
    pts=[(x0+float(a),y0+float(b)) for a,b in re.findall(r'(-?[\d.]+) (-?[\d.]+) [ml]',seg[m.end():])]
    return pts
cols={'Core':'.51 .78 .722','Soft-core':'.902 .365 .427','Shell':'.161 .686 .831','Cloud':'.49 .804 .973'}
lines={}
for c,col in cols.items():
    m=re.search(r'q 1 0 0 1 371\.742 [\d.]+ cm '+re.escape(col)+r' RG \.95 w .*? S Q',s)
    lines[c]=parse_rel(m.group(0))
wb=openpyxl.load_workbook('../../../deliver/Source_Data_split/Source_Data_Fig4.xlsx',read_only=True)
rows=list(wb['Fig.4i_TE_density'].iter_rows(values_only=True))
hdr=rows[0]; sd={}
for r in rows[1:]:
    sd.setdefault(r[0],[]).append(r)
for c in cols:
    v=np.array([r[5] for r in sorted(sd[c],key=lambda r:r[1])])
    pts=np.array(lines[c])
    # x -> bin index
    idx=(pts[:,0]-371.242)/(139.56201/140)-0.5
    ii=np.rint(idx).astype(int)
    pred=283.6333+v[ii]*(364.9493-283.6333)/40
    print(c,len(pts),'idx dev',np.abs(idx-ii).max().round(3),'y dev max',np.abs(pred-pts[:,1]).max().round(3))
from scipy.signal import savgol_filter
def smooth(v):
    r=np.asarray(v,float).copy()
    for a,b,w in ((0,20,7),(20,120,11),(120,140,7)): r[a:b]=savgol_filter(r[a:b],window_length=w,polyorder=2,mode='interp')
    return r
K=(364.9493-283.6333)/40
for c in cols:
    v=smooth([r[5] for r in sorted(sd[c],key=lambda r:r[1])])
    pts=np.array(lines[c]); ii=np.rint((pts[:,0]-371.242)/(139.56201/140)-0.5).astype(int)
    print('smoothed',c,np.abs(283.6333+v[ii]*K-pts[:,1]).max().round(3))
