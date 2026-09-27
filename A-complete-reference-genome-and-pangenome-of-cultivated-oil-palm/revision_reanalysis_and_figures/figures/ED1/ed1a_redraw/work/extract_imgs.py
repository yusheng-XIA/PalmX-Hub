"""Extract the original bitmaps of ED1a (ED1a_study_design_arial.pdf) at native resolution, with their soft masks."""
import fitz, numpy as np
from PIL import Image
SRC = "../../beautify/work/ED1/ED1a_study_design_arial.pdf"
d = fitz.open(SRC); p = d[0]
for info in p.get_image_info(xrefs=True):
    x = info["xref"]
    pix = fitz.Pixmap(d, x)
    sm = d.xref_get_key(x, "SMask")
    if sm[0] == "xref":
        mask = fitz.Pixmap(d, int(sm[1].split()[0]))
        pix = fitz.Pixmap(pix, mask)
    if pix.n - pix.alpha != 3:
        pix = fitz.Pixmap(fitz.csRGB, pix)
    pix.save(f"img/x{x}.png")
    print(x, pix.width, pix.height, pix.alpha, [round(v, 3) for v in info["transform"]])
