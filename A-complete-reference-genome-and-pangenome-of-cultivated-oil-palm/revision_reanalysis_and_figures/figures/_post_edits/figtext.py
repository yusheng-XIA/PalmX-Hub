"""Vector in-figure text edits for the round-1 fixes (text-only redaction; graphics untouched).

Each edit selects one text line (all its spans), removes the text and writes the replacement segments on the same
baseline and direction with the same size and colour (Arial / Arial Bold / Arial Italic embedded).  After saving, the
page text is compared with the original: only the targeted strings may differ.
"""
import fitz, io, math, re, shutil, sys
from pathlib import Path
from PIL import Image

AR = "/System/Library/Fonts/Supplemental/Arial.ttf"
ARB = "/System/Library/Fonts/Supplemental/Arial Bold.ttf"
ARI = "/System/Library/Fonts/Supplemental/Arial Italic.ttf"
FONTS = {"r": AR, "b": ARB, "i": ARI}
norm = lambda t: t.replace("\xa0", " ").replace("\xad", "-").replace("\u037e", ";")


def lines(p):
    for b in p.get_text("dict")["blocks"]:
        for l in b.get("lines", []):
            yield l


def rgb(c):
    return ((c >> 16) & 255) / 255, ((c >> 8) & 255) / 255, (c & 255) / 255


def style_of(span):
    f = span["font"]
    return "b" if "Bold" in f else ("i" if "Italic" in f else "r")


def edit_line(p, select, segments, align="center", tag="x", size=None):
    """select(line_text, line) -> bool; segments: list of (text, style) with style r/b/i, or None to keep styles
    ('same' -> single segment in the style of the first span)."""
    hits = [l for l in lines(p) if select(norm("".join(s["text"] for s in l["spans"])), l)]
    assert len(hits) == 1, (tag, len(hits), [norm("".join(s["text"] for s in l["spans"])) for l in hits])
    l = hits[0]
    sp = l["spans"]
    c, s_ = l["dir"]                                   # baseline direction, y down
    size, col = (size or sp[0]["size"]), rgb(sp[0]["color"])
    ox, oy = sp[0]["origin"]
    old_len = sum(fitz.Font(fontfile=FONTS[style_of(s)]).text_length(s["text"], fontsize=s["size"]) for s in sp)
    assert size is None or True
    if segments and segments[0][1] == "same":
        segments = [(segments[0][0], style_of(sp[0]))]
    new_len = sum(fitz.Font(fontfile=FONTS[st]).text_length(t, fontsize=size) for t, st in segments)
    shift = {"left": 0.0, "center": (old_len - new_len) / 2, "right": old_len - new_len}[align]
    # remove text of the line (shrunken span boxes; text only)
    for s in sp:
        x0, y0, x1, y1 = s["bbox"]
        if abs(c) > 0.99:                              # horizontal
            h = y1 - y0
            p.add_redact_annot(fitz.Rect(x0 + 0.1, y0 + 0.3 * h, x1 - 0.1, y1 - 0.3 * h), fill=False)
        else:                                          # rotated: redact character boxes only (shrunk)
            for rb in p.get_text("rawdict")["blocks"]:
                for rl in rb.get("lines", []):
                    if [round(v, 1) for v in rl["bbox"]] != [round(v, 1) for v in l["bbox"]]:
                        continue
                    for rs in rl["spans"]:
                        for ch in rs["chars"]:
                            cx0, cy0, cx1, cy1 = ch["bbox"]
                            cx, cy = (cx0 + cx1) / 2, (cy0 + cy1) / 2
                            p.add_redact_annot(fitz.Rect(cx - 0.6, cy - 0.6, cx + 0.6, cy + 0.6), fill=False)
            break
    p.apply_redactions(images=fitz.PDF_REDACT_IMAGE_NONE, graphics=fitz.PDF_REDACT_LINE_ART_NONE,
                       text=fitz.PDF_REDACT_TEXT_REMOVE)
    ang = math.degrees(math.atan2(-s_, c))              # counter-clockwise angle as seen on the page
    x, y = ox + c * shift, oy + s_ * shift
    for k, (t, st) in enumerate(segments):
        pt = fitz.Point(x, y)
        if abs(ang) < 0.01:
            p.insert_text(pt, t, fontname=f"F{tag}{k}{st}", fontfile=FONTS[st], fontsize=size, color=col)
        else:
            p.insert_text(pt, t, fontname=f"F{tag}{k}{st}", fontfile=FONTS[st], fontsize=size, color=col,
                          morph=(pt, fitz.Matrix(1, 0, 0, 1, 0, 0).prerotate(ang)))
        w = fitz.Font(fontfile=FONTS[st]).text_length(t, fontsize=size)
        x, y = x + c * w, y + s_ * w
    return "".join(s["text"] for s in sp), "".join(t for t, _ in segments)


def export(pdf_in, pdf_out, png_out=None, docx_out=None, dpi=600, docx_width=2400):
    pg = fitz.open(pdf_out)[0]
    pm = pg.get_pixmap(matrix=fitz.Matrix(dpi / 72, dpi / 72), alpha=False)
    im = Image.open(io.BytesIO(pm.tobytes("png"))).convert("RGB")
    if png_out:
        if str(png_out).endswith((".tif", ".tiff")):
            im.save(png_out, compression="tiff_lzw", dpi=(dpi, dpi))
        else:
            im.save(png_out, dpi=(dpi, dpi))
    if docx_out:
        h = round(im.height * docx_width / im.width)
        im.resize((docx_width, h), Image.LANCZOS).save(docx_out, dpi=(round(docx_width / (pg.rect.width / 72)),) * 2)
    return im


def before_after(pdf_before, pdf_after, rect, out, dpi=400):
    a = fitz.open(pdf_before)[0].get_pixmap(clip=fitz.Rect(rect), dpi=dpi)
    b = fitz.open(pdf_after)[0].get_pixmap(clip=fitz.Rect(rect), dpi=dpi)
    A = Image.open(io.BytesIO(a.tobytes("png"))); B = Image.open(io.BytesIO(b.tobytes("png")))
    C = Image.new("RGB", (A.width, A.height * 2 + 20), "white"); C.paste(A, (0, 0)); C.paste(B, (0, A.height + 20))
    C.save(out)


def text_words(pdf):
    return sorted(norm(fitz.open(pdf)[0].get_text()).split())
