"""Before/after images: whole ED1 side by side (300 dpi) and panel a enlarged side by side (800 dpi)."""
import fitz
from PIL import Image, ImageDraw, ImageFont
SP = "${WORK_DIR}"
OLD = SP + "/deliver/Extended_Data_Figures/Extended_Data_Fig_01.pdf"
NEW = SP + "/fix/ed1a_redraw/out/Extended_Data_Fig_01.pdf"
OUT = SP + "/fix/ed1a_redraw/"
k = 72 / 25.4
font = ImageFont.truetype("/System/Library/Fonts/Supplemental/Arial Bold.ttf", 44)
def render(p, dpi, clip=None):
    pm = fitz.open(p)[0].get_pixmap(dpi=dpi, alpha=False, clip=clip)
    return Image.frombytes("RGB", (pm.width, pm.height), pm.samples)
def side(a, b, labels, out):
    pad, top = 40, 80
    c = Image.new("RGB", (a.width + b.width + 3 * pad, max(a.height, b.height) + top + pad), "white")
    c.paste(a, (pad, top)); c.paste(b, (2 * pad + a.width, top))
    d = ImageDraw.Draw(c)
    d.text((pad, 18), labels[0], fill="black", font=font); d.text((2 * pad + a.width, 18), labels[1], fill="black", font=font)
    d.rectangle([pad - 1, top - 1, pad + a.width, top + a.height], outline="#999999", width=2)
    d.rectangle([2 * pad + a.width - 1, top - 1, 2 * pad + a.width + b.width, top + b.height], outline="#999999", width=2)
    c.save(out, optimize=True)
side(render(OLD, 300), render(NEW, 300), ["Before (deliver, current)", "After (ed1a_redraw)"], OUT + "compare_ED1_before_after.png")
clip = fitz.Rect(0, 0, 96 * k, 85.5 * k)
side(render(OLD, 800, clip), render(NEW, 800, clip), ["Before: panel a (text 2.2-4.2 pt)", "After: panel a (text >= 5 pt)"],
     OUT + "compare_ED1a_zoom.png")
