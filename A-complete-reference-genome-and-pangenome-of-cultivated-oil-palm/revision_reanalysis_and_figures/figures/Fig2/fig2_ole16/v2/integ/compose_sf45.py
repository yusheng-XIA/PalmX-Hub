#!/usr/bin/env python3
"""SF4 (+ d, former Fig. 2c heatmap) and SF5 (+ j, OLE16a upstream classes): vector composition from the pre-integration
backups; writes PDF, 600-dpi PNG and a 2400-px-wide docx PNG for each."""
from pathlib import Path
import fitz
from PIL import Image
Image.MAX_IMAGE_PIXELS = None
H = Path(__file__).resolve().parent; V = H.parent; B = V / 'backup_integrate'
OUT = H / 'sf_out'; OUT.mkdir(exist_ok=True)


def compose(src, pan, name, gap=6):
    s, p = fitz.open(src), fitz.open(pan)
    w, h = s[0].rect.width, s[0].rect.height
    pr = p[0].rect
    out = fitz.open(); pg = out.new_page(width=w, height=h + gap + pr.height)
    pg.show_pdf_page(fitz.Rect(0, 0, w, h), s, 0)
    pg.show_pdf_page(fitz.Rect(0, h + gap, pr.width, h + gap + pr.height), p, 0)
    out.save(OUT / f'{name}.pdf', garbage=4, deflate=True)
    pix = pg.get_pixmap(dpi=600, alpha=False)
    im = Image.frombytes('RGB', (pix.width, pix.height), pix.samples)
    im.save(OUT / f'{name}.png', dpi=(600, 600))
    im.resize((2400, round(im.height * 2400 / im.width)), Image.LANCZOS).save(OUT / f'{name}_docx.png')
    print(name, f'{w/72*25.4:.1f} x {(h+gap+pr.height)/72*25.4:.1f} mm', im.size)


compose(B / 'Supplementary_Fig_04.pdf', V.parent / 'panels/sf4_metab.pdf', 'Supplementary_Fig_04')
compose(B / 'Supplementary_Fig_05.pdf', H / 'sf5j/SF5j.pdf', 'Supplementary_Fig_05', gap=4)
