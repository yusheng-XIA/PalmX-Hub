import fitz
from PIL import Image
Image.MAX_IMAGE_PIXELS=None
from pathlib import Path
H=Path(__file__).resolve().parents[1]
d=fitz.open(H/'Figure2_v3.pdf'); pg=d[0]
print('mm', round(pg.rect.width/72*25.4,1), round(pg.rect.height/72*25.4,1))
pix=pg.get_pixmap(dpi=600, alpha=False)
im=Image.frombytes('RGB',(pix.width,pix.height),pix.samples)
im.save(H/'Figure2_v3.tif', compression='tiff_lzw', dpi=(600,600)); print(im.size)
pg.get_pixmap(dpi=150).save(H/'Figure2_v3_preview.png')
