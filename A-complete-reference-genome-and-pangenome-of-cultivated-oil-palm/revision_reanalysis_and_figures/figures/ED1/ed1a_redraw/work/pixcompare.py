"""Pixel comparison of the new ED1 against the formal deliver version, outside panel a (600 dpi renders of both PDFs)."""
import sys, fitz, numpy as np
SP = "${WORK_DIR}"
OLD = SP + "/deliver/Extended_Data_Figures/Extended_Data_Fig_01.pdf"
NEW = SP + "/fix/ed1a_redraw/out/Extended_Data_Fig_01.pdf"
def rd(p):
    d = fitz.open(p); pg = d[0]
    pm = pg.get_pixmap(dpi=600, alpha=False)
    return pg.rect, np.frombuffer(pm.samples, np.uint8).reshape(pm.height, pm.width, 3).astype(int)
ro, a = rd(OLD); rn, b = rd(NEW)
print("page old", ro, "new", rn, "mm", rn.width / 72 * 25.4, rn.height / 72 * 25.4)
assert a.shape == b.shape
px = 600 / 25.4
# panel-a cell (x 3-95 mm, y 4-84.44 mm) plus the letter 'a' (0-3 mm, 1.5-5 mm)
mask = np.zeros(a.shape[:2], bool)
mask[int(4 * px) - 2:int(84.45 * px) + 3, int(3 * px) - 2:int(95 * px) + 3] = True
diff = (np.abs(a - b).max(2) > 0)
print("changed pixels inside a cell:", int((diff & mask).sum()))
print("changed pixels outside a cell:", int((diff & ~mask).sum()), "of", int((~mask).sum()))
ys, xs = np.nonzero(diff & ~mask)
if len(ys):
    print("outside bbox mm:", xs.min() / px, ys.min() / px, xs.max() / px, ys.max() / px, "max delta", np.abs(a - b).max(2)[diff & ~mask].max())
# where does the new panel a draw? must stay inside the cell
ink = (b.min(2) < 250)
ys, xs = np.nonzero(ink[:int(90 * px), :int(97 * px)])
print("new a ink bbox mm:", xs.min() / px, ys.min() / px, xs.max() / px, ys.max() / px)
