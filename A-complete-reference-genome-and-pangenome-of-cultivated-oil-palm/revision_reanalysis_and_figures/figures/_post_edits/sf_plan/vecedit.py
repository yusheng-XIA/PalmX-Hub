"""Minimal vector text replacement for matplotlib-made PDFs (horizontal or 90-degree text), preserving font style, size,
colour and the text's anchored end. Only the matched spans are redacted (graphics and images untouched)."""
import fitz
FONTS = {"ArialMT": "/System/Library/Fonts/Supplemental/Arial.ttf",
         "Arial-BoldMT": "/System/Library/Fonts/Supplemental/Arial Bold.ttf",
         "Arial-ItalicMT": "/System/Library/Fonts/Supplemental/Arial Italic.ttf"}

def spans(page):
    for b in page.get_text("dict")["blocks"]:
        for l in b.get("lines", []):
            for s in l["spans"]:
                yield s, l["dir"]

def replace(page, old, new, keep="start", where=None, expect=None):
    """Replace every span whose text == old (optionally filtered by where(bbox)). keep: 'start' keeps the reading-start
    position, 'end' keeps the reading-end position (e.g. right-/top-aligned tick labels)."""
    hits = [(s, d) for s, d in spans(page) if s["text"] == old and (where is None or where(fitz.Rect(s["bbox"])))]
    if expect is not None:
        assert len(hits) == expect, (old, len(hits))
    for s, d in hits:
        page.add_redact_annot(fitz.Rect(s["bbox"]))
    page.apply_redactions(images=fitz.PDF_REDACT_IMAGE_NONE, graphics=fitz.PDF_REDACT_LINE_ART_NONE,
                          text=fitz.PDF_REDACT_TEXT_REMOVE)
    for s, d in hits:
        fn = s["font"]; ff = FONTS[fn]; name = {"ArialMT": "AR", "Arial-BoldMT": "AB", "Arial-ItalicMT": "AI"}[fn]
        page.insert_font(fontname=name, fontfile=ff)
        f = fitz.Font(fontfile=ff); size = s["size"]
        lo, ln = f.text_length(old, fontsize=size), f.text_length(new, fontsize=size)
        ox, oy = s["origin"]; c = s["color"]; col = ((c >> 16 & 255) / 255, (c >> 8 & 255) / 255, (c & 255) / 255)
        if abs(d[0] - 1) < 1e-3:          # horizontal
            x = ox if keep == "start" else ox + lo - ln
            page.insert_text((x, oy), new, fontname=name, fontsize=size, color=col)
        elif abs(d[1] + 1) < 1e-3:        # reading upwards (rotated 90)
            y = oy if keep == "start" else oy - lo + ln
            page.insert_text((ox, y), new, fontname=name, fontsize=size, color=col, rotate=90)
        else:
            raise ValueError(("unsupported direction", d))
    return len(hits)

def words(page):
    return sorted(w[4] for w in page.get_text("words"))
