"""Shared helpers for the Fig. 4 candidate edits (PyMuPDF, content-stream level)."""
import re
from pathlib import Path
import fitz

FIX = Path(__file__).resolve().parents[1]
SCR = FIX.parents[1]
ORIG = SCR / "deliver/Main_Figures_revised/Figure4.pdf"
SD_XLSX = SCR / "deliver/Source_Data_split/Source_Data_Fig4.xlsx"
MAIN_FORM = 78            # Illustrator form XObject that holds panels a-i, k
PAGE_H = 620.787          # form y = PAGE_H - page y


def num(v, nd=3):
    s = f"{v:.{nd}f}".rstrip("0").rstrip(".")
    if s.startswith("0."):
        s = s[1:]
    elif s.startswith("-0."):
        s = "-" + s[2:]
    return "0" if s in ("", "-", "-0") else s


def get_stream(doc, xref=MAIN_FORM):
    return doc.xref_stream(xref).decode("latin1")


def put_stream(doc, text, xref=MAIN_FORM):
    doc.update_stream(xref, text.encode("latin1"))


def replace_once(text, old, new):
    n = text.count(old)
    if n != 1:
        raise SystemExit(f"expected exactly one occurrence, found {n}: {old[:80]!r}")
    return text.replace(old, new)


def tt0_widths(doc, xref=10):
    obj = doc.xref_object(xref)
    first = int(doc.xref_get_key(xref, "FirstChar")[1])
    arr = [int(x) for x in re.search(r"/Widths\s*\[([^\]]*)\]", obj).group(1).split()]
    return lambda s: sum(arr[ord(c) - first] for c in s) / 1000.0


def words(page, clip=None):
    return [w[4] for w in page.get_text("words", clip=clip, sort=True)]
